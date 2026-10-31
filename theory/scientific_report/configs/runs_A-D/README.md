# a0 scan on the server: LG (0,0) and (2,2), a0 = 0.1 and 0.2

Moved here on 2026-10-31 from `config/expansion_a0_scan_server/` (with `run_batch.sh`), next to the four configs
it describes; the configs in `config/` were removed then.

Prepared 2026-10-29 on `ad1` from a first pass with 8192 electrons on a 10 w0 screen (the "R10 runs", tables
below). These four runs are for the progress report. a0 = 0.3 is left out: on the R10 screen it was clipped by 11–21%
and its N = 2 had visibly larger sampling noise. The report just says a0 = 0.3 is expected to be worse (stronger
ponderomotive lens).

## Configs

| Config | LG (p, l) | a0 | Screen | Electrons |
|---|---|---|---|---|
| `config_expansion_p0l0_a0_01_R15.cfg` | (0,0) | 0.1 | R_max 15 w0, N_R × N_phi = 128 × 64 | 65536 |
| `config_expansion_p2l2_a0_01_R15.cfg` | (2,2) | 0.1 | R_max 15 w0, 128 × 64 | 65536 |
| `config_expansion_p0l0_a0_02_R25.cfg` | (0,0) | 0.2 | R_max 25 w0, 128 × 64 | 65536 |
| `config_expansion_p2l2_a0_02_R25.cfg` | (2,2) | 0.2 | R_max 25 w0, 128 × 64 | 65536 |

Common to all four:
- **Laser:** pure Gaussian pulse (`laser_flat_duration 0`, `laser_wing_sigma 25 cycles_adim`, `cutoff 4`), circular
  polarization, `w0 = 25 λ`, along −Oz.
- **Beam:** γ ≈ 13.5 head-on, radius from the LG table (`beam_cylinder_radius < 0`).
- **Numerics:** `long_distance_direct_approx` with `trajectory_NT 100`.
- **Screen:** at `Z = 160000 w0 = 4×10⁶ λ`, radius sized from the R10 fits (99.5% containment angle × ~1.5).
- Why the Gaussian pulse and the direct formula: `theory/future_dev/pulse_envelope_and_lg_normalization.md` and
  `theory/future_dev/numerical_parameters_and_cost.md`. Short version: the old flat-top pulse leaked N = 1
  radiation into N = 2, breaking its azimuthal symmetry.

## Running

Inside a compute-node job (`theory/srun.txt`: `srun` inside `tmux` on the head node, `--exclude=hpc-ctrl01,hpc-db01`):

```
git pull
bash theory/scientific_report/configs/runs_A-D/run_batch.sh
```

It builds on that node (`--clean-first`, `-march=native`), then runs the four cases one after another, each on all
CPUs of the allocation. Progress and the config → run-folder mapping go to `~/output/a0_scan_server_logs/index.txt`,
with one solver log per config (the logs of the report's runs A–D are now in
`~/output/progress_report_20261029/logs/a0_scan_server/`). Expected roughly 25–45 min per case on 256 cores (an
extrapolation from `ad1`: each case is 16× an R10 run's electrons × screen points, and an R10 run took 13 min on 32
cores).

## Processing each run

For each run folder `R` from `index.txt`:

```
python3 py_scripts/radiation/plot_all_components.py R
python3 py_scripts/radiation/plot_observables.py R
Z0=$(python3 -c "print(','.join(str(v) for v in range(-1500000, 1500001, 100000)))")
python3 py_scripts/radiation/lg_decomposition.py R --z0=$Z0 --z0-unit lambda
python3 py_scripts/radiation/plot_processed.py radial_phase total 0 1 R
```

- **Use the bounded z0 list** (−1.5×10⁶ … 1.5×10⁶ λ in 10⁵ steps). Write it as `--z0=…` because it starts with a minus
  sign. With the default `--scan-z0`, the N = 1 scan reaches z0 ≈ Z (the screen). There the trial modes are a few
  rings wide (flagged UNDER-RESOLVED), their projections exceed 1, and the "largest overlap" rule picks that
  artefact.
- **For looking at the field, use the `radial_phase` maps** (`png_folder/radiation/processed/radial_phase/`). The raw
  Re/Im/arg maps carry the full `e^{ikR}` (~20 fringes over the screen at N = 2) and look aliased. The physics
  doesn't depend on that.
- **Expected azimuthal index in `E_−`** (spin −1 about lab +z; `J = N(m+σ)` with `m`, `σ` taken about lab +z):
  (0,0): N = 1 `l = 0`, N = 2 `l = −1`; (2,2): N = 1 `l = −2`, N = 2 `l = −5`. `E_+` carries ~10⁻¹² of the power.

Report per case and harmonic, for `E_−`:
- `l` purity (the "azimuthal power" line, ring FFT);
- dominant mode, w0, z0, its share `|a|²/‖E‖²`, and Σ|a|²/‖E‖²;
- containment (power beyond 0.9 R_max; flagged above 1%);
- `Λ_zz/P_z` ÷ theory (`plot_observables.py`'s OAM-to-energy table);
- a look at the `radial_phase` maps and the `lg_w0_z0_map` / `lg_coefficients` figures.

## Comparison: R10 runs on `ad1` (8192 electrons, R_max 10 w0, 64 × 64), 2026-10-29

`E_−`, waist position from the bounded z0 scan. "Beyond" = power beyond 0.9 R_max (above ~1% the screen clips the
field and the projections are underestimated).

**LG (0,0)** (runs `20261029_131618`, `_132913`, `_134538` on `ad1`):

| a0 | N | `l` purity | dominant mode | share | w0 [λ] | z0 [λ] | Σ\|a\|² | beyond | Λ_zz/P_z ÷ theory |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 1 | 0.998 | `LG_0^0` | 0.991 | 21.6 | 6.0×10⁵ | 1.007 | 0.1% | — (theory 0) |
| 0.1 | 2 | 0.943 | `LG_0^{−1}` | 0.887 | 10.8 | 5.0×10⁵ | 0.921 | 1.9% | 1.061 |
| 0.2 | 1 | 0.998 | `LG_0^0` | 0.946 | 9.5 | 3.0×10⁵ | 0.973 | 4.4% | — |
| 0.2 | 2 | 0.989 | `LG_0^{−1}` | 0.817 | 5.0 | 3.0×10⁵ | 0.900 | 14% | 0.994 |
| 0.3 | 1 | 0.990 | `LG_0^0` | 0.871 | 8.2 | 3.0×10⁵ | 0.921 | 14% | — |
| 0.3 | 2 | 0.901 | `LG_0^{−1}` | 0.718 | 4.8 | 2.0×10⁵ | 0.813 | 21% | 0.966 |

**LG (2,2)** (runs `20261029_140004`, `_141400`, `_142701` on `ad1`):

| a0 | N | `l` purity | dominant mode | share | w0 [λ] | z0 [λ] | Σ\|a\|² | beyond | Λ_zz/P_z ÷ theory |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 1 | 0.987 | `LG_2^{−2}` | 0.975 | 25.4 | 1.0×10⁵ | 0.985 | 0.4% | 0.989 |
| 0.1 | 2 | 0.951 | `LG_0^{−5}` | 0.662 | 9.9 | 0 | 0.718 | 1.8% | 0.982 |
| 0.2 | 1 | 0.988 | `LG_2^{−2}` | 0.899 | 24.3 | 1.0×10⁵ | 0.973 | 1.7% | 0.988 |
| 0.2 | 2 | 0.950 | `LG_0^{−5}` | 0.791 | 9.6 | 1.0×10⁵ | 0.882 | 3.6% | 0.971 |
| 0.3 | 1 | 0.985 | `LG_0^{−2}` | 0.791 | 12.3 | 0 | 0.907 | 11% | 0.985 |
| 0.3 | 2 | 0.919 | `LG_0^{−5}` | 0.808 | 9.4 | 1.0×10⁵ | 0.879 | 3.9% | 0.940 |

**What the R10 runs showed, and what the server runs should confirm:**
- **(0,0):** each harmonic is close to one LG mode, `LG_0^0` (N = 1) and `LG_0^{−1}` (N = 2), with the waist not at
  the interaction point. The ponderomotive dressing acts as a lens: the waist shrinks and moves downstream as a0
  grows. The share falls with a0 partly because the R10 screen clipped the wider beams; the larger screens here
  should show how much of that is real.
- **(2,2):** N = 1 copies the laser mode (`LG_2^{−2}` at ~25 λ, little lens, because the (2,2) laser's peak field is
  lower). N = 2 has clean `l = −5` but is a mixture of radial modes, as the coherent-source model predicted at rest
  (`theory/lg_decomposition_implementation.md`, Verification 7). So "almost pure LG" applies to N = 1 of both lasers
  and N = 2 of the Gaussian laser, not to N = 2 of (2,2).
- **Sampling noise:** at 8192 electrons the off-`l` power was mostly an incoherent background in the outer rings,
  not a coherent admixture. One electron on axis is `l`-pure to 10⁻⁹ at a0 = 0.3, NT 100 and 200 alike. 65536
  electrons should cut that background ~8×; at a0 = 0.1, N = 2 of (0,0) should rise above 0.95 in `l` purity.
- **Caveats for the report:** results at one frequency per harmonic; an idealized, fully coherent zero-thickness
  disk; for the LG laser, "a0" means an on-axis peak of √2·a0 (`theory/future_dev/pulse_envelope_and_lg_normalization.md` §2).
