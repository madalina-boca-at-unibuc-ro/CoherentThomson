# Lattice stack: what was implemented

Record of the finished work. The proposal is
[CoherentThomson-lattice-stack-implementation.md](CoherentThomson-lattice-stack-implementation.md); discussion and open
questions are in [lattice_stack-work.md](lattice_stack-work.md).

| Part | Status |
|---|---|
| Preparation: trajectory memory | Done, commit `3052872` |
| Phase 1: uniform stack | Done, commits `2fd28d0`, `ec1517d`, `8ccd0e6`; run by the user on `M` = 1, 2, 4 |
| Phase 2: density model | Done; the user chose moved sheets |
| Phase 3a: `b₁` | Done, reviewed by the user (2026-11-09) |
| Beam distribution rework | Done (2026-11-09), awaiting the user's review; see the last section |
| Phase 3b: `b₂₀` | Not started |

## Preparation: trajectory memory (2026-11-06)

**Problem.** Every electron kept its full trajectory until the end of the run, about 2 MB per electron at 20,000 points.
A 280-sheet stack (301,280 electrons) would have needed about 630 GB.

**Change** (commit `3052872`, branch `lattice-stack`):

- `Electron::compute_trajectory` now reserves the storage, instead of the constructor.
- New `Electron::release_trajectory()` frees it. `run_simulation` calls it after each electron's `compute_radiation`.
- Memory now scales with threads, not electrons. Nothing may read a trajectory after `run_simulation`. `electron.dat` is
  written before it, so it is unaffected.

**Check.** The old and new binaries ran the same configs on 4 threads with GCC 15.2. The configs are based on
`config/config_lattice.cfg`, with a 6 × 6 circular screen and 2 harmonics.

| Case | Electrons | `.dat` files identical |
|---|---|---|
| `simplified`, random beam | 60 | 5 of 5 |
| `direct`, random beam | 60 | 5 of 5 |
| `long_distance_simplified`, random beam | 60 | 5 of 5 |
| `long_distance_direct`, random beam | 60 | 5 of 5 |
| `long_distance_simplified_approx`, random beam | 60 | 5 of 5 |
| `long_distance_direct_approx`, random beam | 60 | 5 of 5 |
| `long_distance_simplified_approx`, polar lattice | 3,319 | 5 of 5 |

Peak memory of the lattice case, from `/usr/bin/time`:

| | Peak memory | Wall time |
|---|---|---|
| Before | 3.4 GB | 7.5 s |
| After | 22 MB | 7.0 s |

The outputs are in `~/output_tests/lattice_stack/memory_fix/{baseline,fixed}/`. The fix has not yet been run on the
cluster or at the size of a real stack.

## Phase 1: uniform stack (2026-11-06)

**Config.** `beam_stack_sheets M` (optional, default 1) and `beam_stack_spacing d` (needed when `M > 1`), with
`beam_distribution polar_lattice`. Model config: `config/config_lattice_stack.cfg`.

**Geometry.** Sheets sit at `z_j = (j − (M−1)/2) d`, before `beam_center` and the beam rotation. A stack covers `H = M d`,
with a first-to-last span of `(M−1) d`.

**Rejected configs.** These stop with an error: `beam_stack_spacing` without `beam_stack_sheets`; `M` below 1 or not an
integer; `M > 1` without a spacing; a stack with a random beam; non-zero `beam_cylinder_height` or `beam_sigma_z`.

**Charge.** Each electron has `charge_weight = 1/M`, so the total charge equals the flat disk's.

| Where | What the weight does |
|---|---|
| `radiation.cpp` | Scales the field once per (frequency, screen point), for `F_l`, `F_s` and `F_b`, in all six modes |
| `form_factor.cpp` | `b = (1/W) Σ w_j e^{iq·r_j}` and a weighted mean momentum |
| `run_log.txt` | Lists sheets, sites per sheet, spacing, `H`, span, represented charge and trajectory count, for `M > 1` |
| Particle motion | Not affected |

`plot_beam_density_z.py` draws a stack as the uniform density it samples, over `H`.

### Checkpoint results

Run A's laser and lattice (`a0 = 0.1`, 1,076 sites per sheet, pure Gaussian pulse), 16 × 8 screen, harmonics N = 1, 2,
`long_distance_direct_approx` unless stated. Full output: `theory/lattice_stack_tests/phase1_results.txt`.

`b` in these tables used `q = k_s − k_L`, as the code did then; since the next commit it uses the exact `q` (below).
`c` is the best complex factor in `F_M = c F_1`. `B` is the discrete `B_M` at the measured slope `2 N k_1/(1+β₀)`.

**`M = 1` reproduces the flat lattice.** Against the pre-change binary, all `.dat` files are byte-identical for the six
radiation formulas (random beam, 60 electrons) and for the lattice (3,319 electrons), with the key omitted and with
`beam_stack_sheets 1`.

**Stack scan, `d = 1e-4 λ_L`, fixed charge.**

| `M` | `H/Λ₁` | `b` | N = 1: `B` | N = 1: `c` | N = 1: `\|c−B\|` | N = 2: `B` | N = 2: `c` | N = 2: `\|c−B\|` |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.145 | +0.974064 | +0.974065 | +0.974065 | 3e-8 | +0.897604 | +0.897604 | 1e-7 |
| 7 | 0.509 | +0.631124 | +0.631129 | +0.631129 | 6e-8 | −0.017446 | −0.017446 | 5e-8 |
| 14 | 1.017 | −0.017000 | −0.016994 | −0.016994 | 5e-9 | +0.017421 | +0.017421 | 5e-8 |
| 21 | 1.526 | −0.209764 | −0.209766 | −0.209766 | 2e-8 | −0.017379 | −0.017379 | 5e-8 |
| 28 | 2.034 | +0.016976 | +0.016969 | +0.016969 | 8e-9 | +0.017320 | +0.017320 | 5e-8 |

Across all 11 stacks of the scan, `|c−B|` is at most 6e-8 at N = 1 and 1e-7 at N = 2. The residual over the screen,
`|F_M − B F_1|/|F_1|`, is at most 7e-8 and 3e-7. `radiation_field.dat` has 7 significant digits, so these may be close to
the output's own rounding.

**The exact minimum separates the two `q`.** With `H = Λ₁ = 2π/|k_s−k_L|`, `|b|` is 6e-9 (limited by the 7 digits of `q` read from
the log), but the field is `c = +7e-6`, which is what `B` predicts at the measured slope. The absolute error of `c` is 6e-9 at
N = 1 and 5e-8 at N = 2, the same at `M` = 14, 28 and 56.

**Spacing refinement at fixed `H` and charge.** The change from `d = 1e-4` to `0.5e-4 λ_L` is the discrete prediction:

| `H` [`λ_L`] | N | `\|F_fine−F_coarse\|/\|F_1\|` | Prediction |
|---|---|---|---|
| 7e-4 | 1 | 4.106e-3 | 4.105e-3 |
| 7e-4 | 2 | 4.525e-4 | 4.525e-4 |
| 14e-4 | 1 | 1.105e-4 | 1.105e-4 |
| 14e-4 | 2 | 4.518e-4 | 4.518e-4 |
| 21e-4 | 1 | 1.365e-3 | 1.365e-3 |
| 21e-4 | 2 | 4.507e-4 | 4.507e-4 |

**Off-centre stack** (`M = 7`, `beam_center_z = 1e-3 λ_L`): the complex `c` matches, to 6e-8 (N = 1) and 5e-8 (N = 2).

**Control, `q d = 2π`:** `|b| = 1.000000000` at `M` = 2, 3, 5, diagnostic only.

**Sheet translation:** see `theory/lattice_stack-work.md`. Relative residual power is at most 1e-12 (N = 1) and 2e-10
(N = 2) across ±140 nm, and `|c| − 1` is at most 2e-6.

**Radiation formulas**, one on-axis point, N = 1, totals:

| Formula | `H = Λ₁`: `c` (`B` = 6.91e-6) | `H = Λ₁/2`: `c` (`B` = 0.642) | Residual |
|---|---|---|---|
| `direct` | 6.908133e-6 | 0.6419986 | 3e-9, 1e-7 |
| `simplified` | 6.908133e-6 | 0.6419986 | 3e-9, 1e-7 |
| `long_distance_direct_approx` | 6.908131e-6 | 0.6419985 | 3e-9, 5e-8 |
| `long_distance_simplified_approx` | 6.908131e-6 | 0.6419985 | 3e-9, 5e-8 |

The endpoint and pulse-edge effect feared in the discussion does not appear at this level: the residual at the minimum
is 2e-8 at N = 1 and 1e-7 at N = 2, about the output's rounding.

### Exact `q` in the form factor (2026-11-06)

`b` now uses `q = k_s − (K_s·P / K_L·P) k_L`, undressed mean `P`, and the log also lists `k_s − k_L`. No `.dat` output changes
(checked on the lattice case). Checks, in `theory/lattice_stack_tests/form_factor_exact_q.py`:

- At `H = 2π/|k_s−k_L|`, `b = 6.900e-6` against the field's `c = 6.909e-6`. The old `q` gave `|b|` about 6e-9.
- At θ = 0.9π, `a0 = 0.2`, small plane-wave disk: `b`'s phase is the conjugate of the field's, to 2e-8 rad for an x shift and
  7e-8 rad for a z shift of 1e-3 λ_L.

### Not tested

- A full 128 × 64 screen, the cluster, GCC 13.3, or an LG laser with `l ≠ 0`.
- Momentum spread in a stack.
- Stacks longer than 56 sheets, or more than `H = 2.03 Λ₁` (the translation test covers ±140 nm).
- `plot_beam_density_z.py` ran on a stack run with the analytic density only, not with the scatter histogram.

### To reproduce

```
cmake -B build/ && cmake --build build/ --clean-first
.venv/bin/python theory/lattice_stack_tests/phase1_translation_slope.py
.venv/bin/python theory/lattice_stack_tests/phase1_stack_scan.py      # about 30 min on 24 threads
.venv/bin/python theory/lattice_stack_tests/q_formula_head_on.py
.venv/bin/python theory/lattice_stack_tests/q_formula_tilted_small_disk.py
.venv/bin/python theory/lattice_stack_tests/q_formula_tilted_wide_disk.py
```

Run them from the repository root. Finished runs in `~/output_tests/lattice_stack/` are reused, so delete them after
changing the solver. These runs came from a dirty build on top of `db70d2b`, and each folder holds its `source.patch`.

## Phase 2: density model (2026-11-06)

**Tool.** `theory/lattice_stack_tests/phase2_density_model.py`, standalone. With no arguments it runs the checkpoint in about 1 min
and writes `phase2_results.txt`; `--zeta` gives one case at any `q`, `N`, sheet count and origin; `--run` compares with a run's log.

**Two options, both at fixed total charge.** With `φ = K z + origin`, `K = 2π/Λ_mod`, uniform sheets `z_j = (j − (SP−1)/2) d`, `d = Λ_mod/S`:

- *Charges:* sheets stay at `z_j`, each with the charge of its bin, from the map's exact crossings of the bin edges on all branches.
- *Moved:* equal charges, `K z_j → K z_j + ζ sin(K z_j + origin)`.

**Checkpoint results.**

| Check | Result |
|---|---|
| Charge | Sums to 1 within 2e-16, no negative bin. At `ζ = 1` the lightest bin is about 0.5 of the mean, the heaviest 23 at `S = 280` |
| Density | Bin charges match the brief's formula away from folds: 3e-8 at `ζ = 0.5` and 2e-8 at `ζ = 1.5`, `S = 4480` (median) |
| Direct integral against Bessel | At most 4e-16, for `n` = 1, 2, 20, 40 and `ζ` = 0, 1, 1.111, 1.5 |
| Uniform control, `ζ = 0` | `|b|` at most 4e-16 |
| Log, uniform stacks | Tool against `run_log.txt`: 4e-8 for `M` = 2, 4, 14, about the 7 digits of the log's `q` |
| Charges by a second route | Nearest-sheet binning of 4e7 phases against exact crossings: at most 3e-7 |

**Errors of the sheet sum, relative to `|b_n|`, origin 0:**

| N (`n`, `Λ_mod = 20Λ`) | `ζ` | `S` | Charges, sheet spacing | Charges, sinc divided out | Charges, 4e7-point binning | Moved |
|---|---|---|---|---|---|---|
| 1 (20) | 1 | 120 | 10% | 5.9% | 3e-8 | 2e-16 |
| 1 (20) | 1 | 280 | 1.6% | 0.80% | 9e-8 | 1e-15 |
| 1 (20) | 1 | 560 | 0.37% | 0.16% | 2e-7 | 2e-15 |
| 1 (20) | 1.111 | 280 | 2.1% | 1.3% | 6e-8 | 2e-16 |
| 2 (40) | 1 | 120 | 44% | 33% | 8e-7 | 3e-15 |
| 2 (40) | 1 | 280 | 7.3% | 4.1% | 1e-6 | 1e-15 |
| 2 (40) | 1 | 560 | 1.6% | 0.80% | 1e-6 | 1e-15 |
| 2 (40) | 1.111 | 280 | 10% | 7.2% | 9e-7 | 2e-16 |

The spacing error dominates; computing the charges adds almost nothing. Near `ζ = 1.5`, `J₂₀(30)` is close to zero (`b = 0.0048`), so
relative errors there are not a useful measure.

**Findings beyond the brief.**

- *Phase origin:* the charges' error depends on where the bins fall relative to the density's cusp. At `n = 20`, `ζ = 1`, `S = 280` it
  is 1.6%, 5.9% and 0.37% at origins 0, 0.7 and −2.3. At origin 0.7 it is mostly a phase: `|b|` is 0.75% low, with a phase error of about
  0.06 rad. Moved sheets stay below 1e-14.
- *Sheet count:* moved sheets fail abruptly. `S = 120` is exact at N = 1, 2 and has 100% error at N = 3 (`n = S/2`). `S = 180` is good to N = 3
  at `ζ` = 1 and 1.111, and `S = 120` reaches only 1.6e-5 at N = 2 for `ζ = 1.5`. Choose `S` from the highest harmonic wanted.
- *Envelope:* on the peak the two agree. Off it they can differ, because moved sheets leave `H` for `ζ > 1` or a non-zero origin (0.10 `Λ_mod` at
  each end for `ζ = 1`, origin 0.7, and 0.16 for `ζ = 1.5`, origin 0.7). There, `|b|` half a line width from the peak is 0.127 moved and 0.112 charges
  (`P = 3`). Folding the moved sheets back into `H` recovers the charges' value, to about 3e-5. At `ζ = 1`, origin 0, nothing leaves `H`
  and the charges converge to the moved value as `S` grows.
- *Replica:* at `q d = 2π`, `|b| = 1` for charges and 0.11 to 0.18 for moved sheets (`ζ = 1`). Moved stacks have no complete replica.

**What this does not show.** The moved sheets' accuracy at `q = n K` is partly by construction: their sum is the rule for the same integral that
defines `b_n`. None of this has gone through the solver. Modulated stacks have not been compared with a log or a field.

**To reproduce:** `.venv/bin/python theory/lattice_stack_tests/phase2_density_model.py`, from the repository root.

## Phase 3a: bunched stack in the solver (2026-11-06)

**Config.** `beam_stack_modulation` (optional) picks the stack: `true` uses the keys below and ignores `beam_stack_sheets` and `beam_stack_spacing`, noting that in the
log; `false` does the opposite. Both sets can stay in a config. A `beam_modulation_*` key without the switch is an error, so it can't go silently unused.

| Key | Meaning |
|---|---|
| `beam_modulation_zeta` | `ζ ≥ 0`; `0` is the uniform control |
| `beam_modulation_period_ratio` | `r = Λ_mod/Λ₁`, any positive real number |
| `beam_modulation_sheets_per_period` | `S`, a whole number |
| `beam_modulation_periods` | `P`; `M = round(S P)`, with a warning if that isn't a whole number of periods |
| `beam_modulation_phase` | Optional, default 0: `φ = K z + phase`, `z` from the stack centre |

**Geometry.** `Λ₁ = 2π/|q·axis|`, from the exact `q` at harmonic 1 along the detector normal, with the configured `average_p`.
`d = Λ_mod/S`, `H = M d`. Sheets start at `(j − (M−1)/2) d` and move by `(ζ/K) sin(K z + phase)`, with charge `1/M` each.
`main.cpp` now creates the detector before the beam. Folding sheets back into `H` is left out, as agreed.

**Log.** `Λ₁`, `Λ_mod`, `q_transverse/|q_axial|`, `S`, `H/Λ_mod`, `d`, `H`, how far sheets leave `±H/2`, and for each harmonic the sheet sum
against the map's phase integral. It warns on a non-whole `P` and on `S < 3 N_max r`.

**Plot.** `plot_beam_density_z.py` draws a bunched stack's sheets against the drift map, using `Λ_mod`, `H` and the phase from the log. Bins cover only the sheets and average about
four of them, at least 5 bins, since finer ones hold 0, 1 or 2 discrete sheets and show gaps that aren't in the density. With few sheets the bars are coarse (14 sheets: 5, 1, 2, 1, 5 per bin).

**Unchanged beams.** Against the pre-change binary, all `.dat` files are byte-identical for a flat lattice, an `M = 14` stack, and a random
beam with `simplified`, `direct` and `long_distance_direct_approx`. Thirteen bad configs stop with an error: ten listed in
`phase3a_noop_and_rejects.py` and three with a required key missing. The switch came after the checkpoint runs below, which set only `beam_modulation_*`. Adding it was checked
to leave a bunched stack's `.dat` files byte-identical, so those results stand; `phase3a_checkpoint.py` now sets the switch and will recompute instead of reusing them.

### Checkpoint results

Run A (`a0 = 0.1`, 1,076 sites per sheet, 16 × 8 screen, `long_distance_direct_approx`), fixed total charge. `c` is the best complex factor in
`F = c F_flat`; the prediction is `(-1)^n J_n(nζ) e^{+inφ₀}`. Residual: `|F − c F_flat|/|F_flat|`. Full output: `lattice_stack_tests/phase3a_results.txt`.

| Stack | N | Prediction | `\|c−pred\|` | `\|c\|²` | Residual |
|---|---|---|---|---|---|
| Uniform, `ζ = 0`, S = 14, P = 3 | 1 | 0 | 1e-8 | below 1e-15 | 2e-8 |
| | 2 | 0 | 6e-8 | below 1e-14 | 1e-7 |
| `ζ = 1`, S = 14, P = 3 | 1 | −0.440051 | 2e-8 | 0.19364 | 2e-8 |
| | 2 | +0.352834 | 7e-8 | 0.12449 | 2e-7 |
| `ζ = 1`, S = 28, P = 3 | 1 | −0.440051 | 2e-8 | 0.19364 | 2e-8 |
| | 2 | +0.352834 | 6e-8 | 0.12449 | 1e-7 |
| `ζ = 1.5`, S = 14, P = 1 | 1 | −0.557937 | 2e-8 | 0.31129 | 2e-8 |
| | 2 | +0.486091 | 2e-7 | 0.23628 | 1e-7 |
| `ζ = 1`, S = 14, P = 1, phase 0.7 | 1 | −0.336569 − 0.283488i | 3e-8 | 0.19364 | 2e-8 |
| | 2 | +0.059970 + 0.347700i | 6e-8 | 0.12449 | 2e-7 |

- **Grid:** `|F(S=28) − F(S=14)|/|F_flat|` is 2e-10 at N = 1 and 3e-8 at N = 2 (`ζ = 1`, P = 3). This replaces the brief's 1 nm against 0.5 nm.
- **Power against the flat disk:** `|J₁(1)|² = 0.1936` and `|J₂(2)|² = 0.1245`, against about zero for the uniform stack of the same length.
- **Beam file:** `b` recomputed from `electron_beam_scatter.dat` and the log's `q` agrees with the log's `b` to at most 3e-7, the limit of the 7 digits of `q`.
- **Radiation formulas** (one on-axis point, N = 1, `ζ = 1`): `direct`, `simplified` and both approximate long-distance modes give `|c−pred|` at most 4e-8.
- **`a0 = 0.2`** with the same keys: `Λ₁` becomes 1.3969e-3 `λ_L` (13.97 nm), and at N = 1 `|c−pred|` is 2e-8.
- **Tilted, θ = 0.9π** (small plane-wave disk, `a0 = 0.2`, `q_transverse/|q_axial|` = 4.3e-4): `|c−pred|` is 9e-9 and 6e-8.

**The switch** (`phase3a_switch.py`, against the solver just before it): `true` with both key sets, `false` with both key sets, and no switch with no modulation keys each give
all 5 `.dat` files byte-identical to a run with only the keys in use.

**Preview of `b₂₀`** (`r = 20`, `S` = 120 and 240, 64 sites per sheet, 16 × 8 screen, residual 1e-7 to 8e-7):

| `ζ` | N | Prediction | `c` | `\|c−pred\|` |
|---|---|---|---|---|
| 1 | 1 | +0.164748 | +0.164748 | 5e-8 |
| 1 | 2 | +0.130781 | +0.130780 | 4e-7 |
| 1.111 | 1 | +0.243302 | +0.243302 | 1e-7 |
| 1.111 | 2 | +0.157086 | +0.157087 | 8e-7 |

**The user's run on a 64 × 64 screen** (`~/output/lattice_stack/20261106_164759`, full 1,076-site lattice, `ζ = 1`, ratio 1, S = 14, P = 1, 15,064 trajectories, 1,705 s on 24 threads), against the
flat run `20261106_140716`, whose config differs only in the stack keys:

| | N = 1 | N = 2 |
|---|---|---|
| Predicted `c` | −0.440051 | +0.352834 |
| `\|c−pred\|` | 2e-8 | 7e-8 |
| Residual after fitting `c` | 4e-8 | 8e-8 |
| `P_z` bunched / flat, from the logs | 0.19364 | 0.12449 |
| Predicted, `J_N(N)²` | 0.19364 | 0.12449 |

`w_eff` equals the flat run's at both harmonics (76.88 and 110.82 `λ_L`), and `Λ_zz/P_z` at N = 2 stays −0.1512612, against −0.1512611 in theory. Bunching only rescales the flat disk's field, as expected in this cold, equal-momentum
model. These are ratios; absolute intensities still wait for the single-electron benchmark in `CLAUDE.md`. The starburst at the centre of the flux plots is probably the polar grid cells being drawn, which was not confirmed. It isn't caused by
bunching, given the residual above.

**On-axis N = 2 is not a useful test.** At one on-axis point, `|c−pred|` at N = 2 is 5e-4 to 1e-2, and the four formulas disagree by about 1e-3. There, the flat
field is about 2.6e4 times below the screen's rms (8.6e-23 against 2.3e-18), so the deviation of about 1.6e-24 is about 7e-7 of that rms. It doesn't change with `S`. Its origin
is not established, and on the screen it has no visible effect.

### Not tested

- The full 1,076-site `b₂₀` stack and a 128 × 64 screen (`lg_decomposition.py` on a bunched stack, in particular), the cluster, GCC 13.3, an LG laser with `l ≠ 0`.
- Momentum spread or energy modulation: every sheet keeps `average_p`.
- A detector with `q·axis < 0`. The code handles the sign, but no run covers it.
- Off-peak field scans of bunched stacks, in particular non-integer `r` and a non-whole `P`.

### To reproduce

```
.venv/bin/python theory/lattice_stack_tests/phase3a_noop_and_rejects.py <solver built at 8ccd0e6>
.venv/bin/python theory/lattice_stack_tests/phase3a_switch.py <solver built before the switch; that build no longer exists>
.venv/bin/python theory/lattice_stack_tests/phase3a_checkpoint.py     # about 20 min on 24 threads
.venv/bin/python theory/lattice_stack_tests/phase3a_on_axis_n2.py
```

## Beam distribution rework (2026-11-09)

Plan and decisions: "Beam distribution rework" in [lattice_stack-work.md](lattice_stack-work.md). The Phase 1 and 3a
sections above describe the keys as they were; this section supersedes them.

**Layout.**

| | Random (`beam_distribution random`) | Lattice stack (`polar_lattice`) |
|---|---|---|
| Along the axis | Uniform in `H` (`beam_cylinder_height`), plus Gaussian wings `beam_sigma_z` | `M` sheets (`beam_stack_sheets`) at `z_m = (m − (M−1)/2) H/M`, charge `1/M` each |
| `H = 0` | Flat disk | One sheet; `M > 1` is reduced to 1, with a warning |
| Wings | Allowed | Refused (`beam_sigma_z` must be 0) |
| `M = 1`, `H > 0` | (no `M`) | Warning, `H` set to 0 |
| Other type's keys | `beam_lattice_spacing`, `beam_stack_sheets` ignored | `beam_particle_count` ignored |

**Modulation**, for either: `beam_modulation true` (optional, default false), `beam_modulation_zeta`,
`beam_modulation_period` (a length, unit required), optional `beam_modulation_phase`. Each electron moves by its own
position, `z → z + (ζ/K) sin(K z + phase)`, `K = 2π/Λ_mod`; for a random beam after the wings. The map draws no random
numbers, so the RNG stream is the unmodulated beam's. A beam with no extent along its axis gets a warning that the
modulation only translates it.

**Unit `lambda_1`** `= 2π/|q_axial|`, the axial part of `b(q)`'s `q` (harmonic 1, detector normal, undressed `P`):
accepted by `beam_cylinder_height`, `beam_sigma_z`, `beam_modulation_period`; refused elsewhere. Parsing keeps such a
length in units of `Λ₁`; `Particle::resolve_beam_lengths` converts it in the beam factory and the run log, the two places
with the laser and the detector. It throws when a config uses the unit and `q` is perpendicular to the beam axis.

**Removed keys**, refused with the replacement named: `beam_stack_spacing`, `beam_stack_modulation`,
`beam_modulation_sheets_per_period`, `beam_modulation_periods`, `beam_modulation_period_ratio`. The old `S`, `P`, ratio
`r` are now `M = round(S P)`, `H = (M r / S) lambda_1`, `beam_modulation_period r lambda_1`.

**Run log.** `Λ₁` always (or "undefined"); `H` and `sigma_z` resolved, with the `lambda_1` value given; `d = H/M`; the
notes as `WARNING` lines. With the modulation: `Λ_mod`, the ratio, the periods `H/Λ_mod`, and per harmonic the sum
`(1/count) Σ exp(i q_N z_j)` over the sheets (lattice) or every electron (random, with `1/√N_e` printed) against the
cold-beam value. The cold value needs integer `n = N Λ_mod/Λ₁` and whole periods (each within a relative `1e-9`) and no
wings. `write_run_log` now takes the beam, to read the electrons' axial positions.

**Code.** `IoUtils::parse_cylinder_beam_params` and `CylinderBeamParams` (`io_utils.hpp`); `Particle::axial_translation`,
`resolve_beam_lengths`, `modulated_z`, `sheet_positions`, `axial_positions`, and the reworked `modulation_geometry`,
`modulation_harmonics`, `modulation_warnings` (`electron_factory.cpp`); `run_log.cpp`.
`py_scripts/particle/plot_beam_density_z.py` reads `H` and `sigma_z` from `run_log.txt` and plots a modulated random
beam too; runs from before the rework still plot. `config/config_lattice_stack.cfg` is migrated (same beam).
`phase1_stack_scan.py`, `form_factor_exact_q.py`, `phase3a_checkpoint.py` use the new keys (not rerun);
`phase3a_noop_and_rejects.py` and `phase3a_switch.py` are marked historical (old keys only).

**Checks** (`lattice_stack_tests/rework_checks.py`, results in `rework_results.txt`; old and new solvers built on the
same compute node):

1. Old solver against new, same beam: all `.dat` files byte-identical for the flat lattice, a uniform 14-sheet stack
   (`d = 1e-4 λ` against `H = 1.4e-3 λ`), bunched stacks `S = 14, P = 1` and `S = 14, P = 2, phase π, ζ = 1.111`
   (against `H` and the period in `lambda_1`), and random beams under three formulas and with `sigma_z`.
2. Modulated random beam (`ζ = 1`, `H = Λ_mod = Λ₁`): `|sum − cold|` at N = 1, 2 is 0.5 to 1.6 times `1/√N_e`
   at `N_e` = 2000 and 8000, two seeds each: sampling noise only. The density plot follows the drift map.
3. All 15 refused configs stop with the intended message, including `q ⊥` axis with `lambda_1`.
4. The corrections and warnings fire as intended, with the electron counts they imply.

### Sheets per period against aliasing (2026-11-09)

**Rule.** A lattice stack needs at least

  `S ≥ n (1 + ζ) + 5 n^(1/3)` sheets per modulation period, `n = N_max Λ_mod/Λ₁`

for an error of `b_n` below about `1e-4` (relative, or `1e-7` absolute for harmonics below `1e-3` of full bunching).
In sheets per `Λ₁` at ratio 1: `N_max (1 + ζ) + 5 N_max^(1/3)`. For `1e-6` the margin is about 7 instead of 5.
`Particle::modulation_warnings` now warns below it, in the log and on the terminal. It replaces the earlier
`S ≥ 3 N_max Λ_mod/Λ₁`, which allowed 3 sheets at `n = 1`, `ζ = 1` (error about 26%) and fell short at `ζ = 1.5` from
`n ≈ 20` on.

**Why.** One period of `S` sheets at equally spaced `φ₀`, moved to `φ₀ + ζ sin φ₀`, makes `b_n` an `S`-point sum of
`exp(i n (φ₀ + ζ sin φ₀))` = `Σ_k J_k(nζ) exp(i (n + k) φ₀)`. The sum keeps every term with `n + k ≡ 0 (mod S)`: the
exact `(−1)ⁿ Jₙ(nζ)` plus the aliases `J_{mS−n}(nζ)`, `m ≠ 0`, the largest being `J_{S−n}(nζ)`. `J_ν(x)` collapses
once `ν` passes `x` by the transition width `~x^(1/3)`, hence `S − n > nζ + c n^(1/3)`. At `ζ = 0` this is the
grating condition `d < Λ_N` (no Bragg lobe from the spacing). A random beam has no aliasing, only `1/√N_e` noise.

**Check** (`lattice_stack_tests/sheets_per_period_aliasing.py`, `sheets_per_period_aliasing_results.txt`): exact
minimum `S` (worst of phase 0, π/2, π) for `n` = 1 to 60 and `ζ` = 0.5, 1, 1.111, 1.5. The rule is within one sheet of
the `1e-4` column everywhere, except `ζ = 1.5`, `n = 20` (`|b|` ≈ 5e-3): 64 against 68, an error about `1e-2` there.
Examples at `1e-4`: `ζ = 1`, n = 1, 2, 5, 20, 40: 7, 10, 18, 52, 95 sheets. Phase 3b (`ratio` 20, `ζ` = 1.111):
N = 1, 2, 3 need 55, 100, 146 per period.

