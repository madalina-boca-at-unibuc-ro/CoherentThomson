# Lattice stack: what was implemented

Record of the finished work. The proposal is
[CoherentThomson-lattice-stack-implementation.md](CoherentThomson-lattice-stack-implementation.md); discussion and open
questions are in [lattice_stack-work.md](lattice_stack-work.md).

| Part | Status |
|---|---|
| Preparation: trajectory memory | Done, commit `3052872` |
| Phase 1: uniform stack | Done, awaiting the user's review; not committed yet |
| Phase 2: density model | Not started |
| Phase 3a: `b₁` | Not started |
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
