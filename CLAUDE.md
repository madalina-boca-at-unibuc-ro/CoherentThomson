# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

Pipeline (`src/app/main.cpp`): build a config-driven laser and detector, generate the beam, integrate up to 10
electrons' trajectories for `electron.dat`, compute the form factor `b(q)` (diagnostic), run
`Simulation::run_simulation` over the whole beam, and write `.dat` files plus `run_log.txt` to
`~/<output_folder>/YYYYMMDD_HHMMSS/` (`output_folder` is relative to `$HOME`). `py_scripts/` plots the outputs.

**Open issues** (detail in `theory/dev_notes.md` unless noted):
- The single-electron Thomson-dipole benchmark (electron at rest, full-4π spherical detector, weak field →
  `1+cos²θ` for circular polarization) is **not** reproduced yet. Don't trust absolute intensities/angular patterns
  for complex setups until it is.
- `plot_observables.py`'s `Lambda_zz/P_z` sign is fixed (`lz_sign`); theory column `L_N/omega_N`,
  `L_N = N(m+σ)−σ`, `σ = 2 Im(ζ₁*ζ₂)`. It was ~5% low at N = 2, likely the old centered-difference `∂φ` (replaced
  by an exact FFT derivative on 2026-10-29): re-run `plot_observables.py` on older runs before reading OAM numbers.
- `b(q)` (`Simulation::compute_form_factor`, spec `theory/thomson_form_factor_diagnostic.md`): one scalar at
  harmonic 1 along the detector normal, using the beam's actual mean momentum. Possibly inconsistent away from
  head-on collisions: don't rely on it in tilted geometries.
- `LaguerreGaussLaser` matches the standard LG mode at the waist and at `Z = 10 z_R` (to `1e-7`, after the
  2026-10-27 curvature-sign fix), but isn't physics-reviewed against the literature (Allen et al. 1992; Siegman
  ch. 17). `PhysUtils::non_linear_Thomson_formula` is flagged in-code as needing review.
- **The flat-top envelope leaks N = 1 radiation into N = 2** (N = 2's off-axis crescent and ~9% `l = 0`
  admixture). A pure Gaussian (`laser_flat_duration 0`) removes it. The LG on-axis field is √2 × the plane wave's at
  the same `a0`. Both: `theory/future_dev/pulse_envelope_and_lg_normalization.md`. Unacted plans:
  `theory/future_dev/README.md`.
- **Ponderomotive dressing acts as a lens** (electrons at rest): the drift ∝ local `⟨a²⟩` dents the disk, adding a
  radial phase that doesn't vanish as `Z → ∞`, so `N ≥ 2` isn't a single LG mode. `lg_decomposition.py`'s
  prediction with the dressing factor `d_N(r)` reaches fidelity 0.9999/0.998 (N = 1/2) at `a0 = 0.1`; tests 2–4
  open. `theory/future_dev/ponderomotive_lens.md`; its test configs write to `~/output_tests/ponderomotive_lens/`.

## Machines

Two laptops share a checkout through Dropbox; two servers sync through git. Python is `.venv` on the laptops and
Anaconda on the servers (editor type-checking config points at `.venv` only). Compilers differ between machines.

**Slurm cluster:** head node `hpc-ctrl01`, compute nodes `hpc-node[01-20]` (256 CPUs, GCC 13.3), shared home.
Commands, partition limits and the long-run recipe: `theory/srun.txt`.
- **Always `--exclude=hpc-ctrl01,hpc-db01`**: every CPU partition contains them, and a job on the head node looks
  like the login shell (check `squeue -u $USER`, not `hostname`).
- Compute nodes have only Anaconda's CMake **3.21**: don't raise `cmake_minimum_required` or the presets' minimum.
- Build on the compute node that runs the job (`-march=native`). Claude's shell may be on the head node while the
  user's terminal is in a job: check `hostname` before building.
- `num_threads 0` follows the job's allocation, not the node's 256 CPUs.
- Long runs: `srun ... --pty bash` inside `tmux` on the head node; never `exit` in that tmux (it ends the job).
  `radiation_field.dat` and `run_log.txt` appear only when the run ends.

## Commands

```
cmake -B build/ && cmake --build build/                    # build (out-of-source enforced)
cmake --preset release && cmake --build --preset release  # same via preset; `debug` preset opts out of Release
python3 py_scripts/compile_and_run.py <config_file>        # configure + build + run
./bin/coherent_thomson_solver config/config.cfg            # run (config argument required)
```

**Always built as `Release`** (`-O3 -march=native`, C++20): `CMakeLists.txt` overrides any other
`CMAKE_BUILD_TYPE` unless `-DCOHERENT_THOMSON_ALLOW_NON_RELEASE=ON` (set by the `debug` preset; stays cached until
`=OFF`). Debug is ~12x slower: check the configure step's `Build type:` line before trusting a timing.

**Build traps, check these first:**
- **Incremental rebuilds have missed header dependencies**, possibly because Dropbox syncs `build/` and `bin/`
  between laptops. After editing a widely-included header, or if a build acts strangely, use `--clean-first`.
- **Every build writes the same `bin/coherent_thomson_solver`**, and a no-op build doesn't relink. To compare builds,
  `--clean-first` and copy each binary out. Both presets share `build/`.
- **GCC 15.2's SLP vectorizer miscompiles the long-distance radiation functions.** Both carry
  `__attribute__((optimize("no-tree-slp-vectorize")))` (GCC only). Keep it on every machine and compiler version.
  `theory/implementation_details.md` §6.

**Build provenance:** the binary embeds git commit, dirty state, compiler, flags and host (`cmake/build_info.cmake`,
`Core::BuildInfo`), printed in `run_log.txt`. A dirty build writes `<run>/source.patch` (apply to the recorded
commit). `require_clean_build true` refuses to run a dirty build. Git state is read at build time: commit, then build.

**Iterating on code changes:** copy a config (don't edit `config/config.cfg` unless the task is about it) and drop
`beam_particle_count` to e.g. `50`. Configs:
- `config/config.cfg`: default, for current experiments; parameters change freely.
- `config/config_cross_check.cfg`: matched to the Python reference's `main.py` `INPUTS`. Use it for every
  C++-vs-Python comparison; update it when `INPUTS` change.
- `config/config_spectrum.cfg`: `dense_frequency_spectrum true` on one screen point (`plot_point_spectrum*.py`).
- `config/config_expansion*.cfg`: LG laser, large fine circular screen, for LG-mode decomposition
  (`theory/LG_mode_decomposition.md`). `_rest*` = electrons at rest (65536 electrons, ~27 min on 256 CPUs). Their
  header comments lag behind the keys: read the keys (or the run's `run_log.txt`).
- `config/expansion_a0_scan_server/`: `a0` scan for LG-mode purity; its `README.md` covers running (`run_batch.sh`),
  processing and results. `expansion_a0_scan_wide/` holds one config. The removed R10/R20 set is in commit `aa12078`.
- `theory/scientific_report/`: the 2026-10-29 group progress report; its `README.md` maps every quoted result to a
  run and a config (`configs/`, frozen copies) and rebuilds figures/PDF. Its runs live in
  `~/output/progress_report_20261029/` (manifest, logs, `SHA256SUMS`). If a quoted result changes, update the report.
- `theory/future_dev/tests/`: test inputs for the future-dev plans, all writing to `~/output_tests/...`:
  `one_electron_on_axis{,_spectrum}.cfg` (seconds per run), `nt_scan.py` (`trajectory_NT` convergence, runs the
  binary itself, ~3 min), and `bunch_length/`, `ponderomotive_lens/` (each with configs and a `README.md`).

Formatting: `clang-format -i` on touched files (Google style, 120 cols). **No test suite**, no lint command.

**Verification** (stands in for tests; use a small copied config):
- **`simplified` vs `direct`** (`radiation_formula`), comparing totals: the main check. They agree to `O(dτ²)`
  (4x per doubling of `trajectory_NT`) down to a 1 λ screen; a `1/(kR)` mismatch means a short-range bug.
- **Rotation covariance:** rotating laser, beam momentum and forward detector by `R` gives `F → Λ F Λᵀ` (6 digits).
- **Spherical screen:** `B` exactly ⊥ `n0`, `E` not (`|E_r|/|E_θ|` ~ 0.1), via `plot_spherical_components.py`.
- **Byte-identity** for changes that should be no-ops: diff the `.dat` outputs with the same config and `random_seed`.
- **Python reference** with `config/config_cross_check.cfg` (Physics notes).

### Plotting

Every script plots the latest run by default; the optional trailing argument is a **run folder**, not a `.dat`.
"Latest" is looked up under `config/config.cfg`'s `output_folder` (`run_output_utils.py`), whatever config made the
run: for a run under `~/output_tests/` or another folder, pass the run folder explicitly.
PNGs go to `<run>/png_folder/<module>/` (`MODULE_NAME` comes from the script's folder; don't hardcode it).
`plot_observables.py`, `plot_point_spectrum_observables.py` and `lg_decomposition.py` also **write** their own
sections into `run_log.txt`. Conventions and internals: `py_scripts/PLOTTING_NOTES.md`.

```
python3 py_scripts/{laser,detector,particle}/<script>.py   # some need plot_*=true flags in the config
python3 py_scripts/radiation/plot_field.py <long|short|boundary|total> <mu> <nu> [--incident] [run_folder]
python3 py_scripts/radiation/plot_all_components.py [--incident]
python3 py_scripts/radiation/plot_processed.py radial_phase <range|all> <mu> <nu>   # flat detectors only
python3 py_scripts/radiation/plot_point_spectrum{,_all,_observables}.py           # dense_frequency_spectrum runs
python3 py_scripts/radiation/plot_spherical_components.py <long|short|boundary> <E|B> <r|theta|phi>
python3 py_scripts/radiation/plot_observables.py [--incident]   # rect./circular facing ±Oz (--incident: any)
python3 py_scripts/radiation/lg_decomposition.py [--incident] [--w0 ... --w0-unit lambda] [--scan-z0] [--self-test]
                                                 # options: theory/lg_decomposition_reference.md
```
`--incident` uses `incident_field.dat` (analytic laser field; rectangular/circular detectors only). `mu == nu` is
rejected; `mu > nu` uses `F^{nu mu} = -F^{mu nu}`. `E_i = c·F^{i0}`, `B = (F^{32}, F^{13}, F^{21})`.

## Architecture

### Layout

- `main.cpp` only orchestrates; all logic is in the static library `coherent_thomson_core` (`src/core/`, headers in
  `src/core/include/<module>/`). Namespaces are `Core::<CamelCase of folder>`. `PhysUtils::AtomicUnits`:
  `c = 137.035999084`, `q_0 = -1` (signed charge), `e_0 = +1`, `epsilon_0 = 1/(4π)`.
- **Class + factory** (`laser/`, `detector/`): `create_<thing>(const ConfigMap&, ...)` in `<thing>_factory.cpp`
  holds *all* config parsing and units. `particle/`'s factory is `generate_cylinder_beam`. New component: class +
  factory, add the `.cpp` to `src/core/CMakeLists.txt`, wire it into `main.cpp`.
- `py_scripts/<module>/` mirrors the C++ modules; scripts run as plain files and import `utils.*` via `sys.path`.
- **C++/Python copies, keep in sync:** constants (from `phys_utils.hpp`), `convert_unit_to_number` (Python in
  `faraday_frame_utils.py` *raises* on unknown units where C++ warns), and `rotation_matrix_from_direction`.

### Frames and geometry

- **Baseline** (`config/config.cfg`): electrons and screen along +Oz, laser along −Oz (`laser_direction_theta 1.0 pi`).
- **Only the laser moves** with `laser_direction_theta/phi`; beam and detector keep their own lab geometry, so a
  tilt is a real change of collision geometry. `epsilon_1/2 = R x̂, R ŷ` with `R =
  MathUtils::rotation_matrix_from_direction(n)` (proper, so helicity and LG `l` sign are relative to `n`; for
  `θ=π`: `epsilon_2 = −ŷ`). That function snaps to ±Oz only within ~1e-12 rad: don't loosen it, small tilts matter.
- `laser_field*.dat` and `incident_field.dat` are in the **laser's own frame** `(epsilon_1, epsilon_2, n)`.
  `incident_field.dat` is at `z = Z` downstream of the waist on a beam-sized grid (`Radiation::incident_screen`);
  the run writes `incident_config.cfg` for the scripts. Pre-2026-10-27 runs differ (dev notes).
  `radiation_field.dat` is lab-frame components on the detector's local grid (mirrored `(x, −y)` facing −Oz).
- The beam axis follows the mean momentum; `beam_center_x/y/z` apply **before** that rotation.
- `beam_distribution polar_lattice` (with `beam_lattice_spacing`) places deterministic equal-area rings, ignoring
  `beam_particle_count`: no speckle, far fewer electrons for coherent tests.
- `beam_sigma_z` smooths the cylinder ends with a Gaussian z displacement (`theory/gaussian_smoothed_electron_cylinder.md`);
  `0` is byte-identical to the plain cylinder. A fixed `random_seed` with fewer electrons gives the first electrons
  of a larger beam: use different seeds for independent samples.
- Negative `beam_cylinder_radius` = `IoUtils::lg_beam_radius_in_w0(p, l) · w0` (LG only, `p ≤ 4`, `|l| ≤ 6`); the
  table is copied from `theory/table of beam electron radius in units of w_0 for different p,l.md`; keep in sync.
- `CircularDetector`'s radial grid is **equal-area**: `r_i = sqrt(R_min² + i*(R_max²-R_min²)/(N_R-1))`.
- A detector with grid counts `1` is one exact point.

### Radiation calculation (the hot path)

Modes, formulas, optimizations, precision and timings: `theory/implementation_details.md`. **Keep it in sync when
changing `radiation.cpp`.** Formulas: `theory/FT_Faraday_tensor-direct_and_simplified_forms-v2.md`; long-distance
modes: `theory/long_distance_direct_simplified_coding_guide.md`; OAM: `theory/numerical_calculation_of_angular_momentum.md`.

- **Three terms** (`LR_`/`SR_`/`BR_` columns): long-range, short-range, boundary. **Only the sum is physical** (e.g.
  on-axis `F_l^{03}` and `F_b^{03}` cancel). In the LD modes the split differs, so compare only totals.
- `radiation_formula` (required): `simplified`, `direct` (`F_b ≡ 0`), `long_distance_{simplified,direct}` (no
  `F_s`, error ~`1/D`), `long_distance_{simplified,direct}_approx` (stepped phase, fastest, ≤1e-7 from unstepped).
- `general_factor` uses the **signed** `q_0`; direct `F_s` has an explicit `c²` to compensate.
- Trapezoid weights apply to `F_l`/`F_s` only; `F_b` is an unweighted endpoint evaluation.
- **Load-bearing:** prefactors stay in `radiation.cpp`'s anonymous namespace; `radiation_phase_argument` stays
  **frequency-independent** (frequency recurrence over evenly spaced `frequencies_list`); screen point outer, tau
  inner; in LD modes the constant phase is applied at reconstruction, not in the tau loop; stepped loops must stay
  vectorizable (`-fopt-info-vec`).
- **Phases carry only `e^{ik(R - |X|)}`** (computed without cancellation); `run_simulation` applies `e^{ik|X|}`.
  Never put `R` or `|x_0|` back into a phase: rounding `kR ~ 1e9–1e17` per electron destroys coherence
  (`implementation_details.md` §4.5; don't implement `stabilize_phase_axis_aligned_detector.md`).
- **Threading:** contiguous electron chunks, one `std::thread` each, per-thread `PackedRadiationField` reduced after
  `join`; `compute_radiation` must not write shared state. `num_threads 0` = affinity mask, capped at that and the
  electron count (check "Using N threads"). "Time per electron" is CPU time, not wall-clock.

### Frequencies

- `frequencies_list[i] = non_linear_Thomson_formula(k1, q, n2, N_harmonics_min + i)`, `q` dressed with
  `<a²> = a0²/2` for every polarization. **Don't change it to `a0²`** (a past bug).
- `dense_frequency_spectrum=true`: linear scan of `N_omega` points at one screen point (bounds in `omega_laser` or
  `first_harmonic_frequency` only).
- Internally frequencies are `omega/c`; `radiation_field.dat`'s `omega` column is `omega/omega_1`.

### Laser

- `get_faraday_tensor` builds on the pure-virtual `complex_amplitude` (`{amplitude, d/dx, d/dy}`). `Ez`/`Bz` are
  first-order paraxial (`div E` residual ~1e-6 to ~1e-2 is expected). Carrier sign `exp(-iφ)`, the `+i/k`
  coefficient, Gouy sign and curvature sign (`+kρ²/2R`) must flip together.
- `zeta_1/zeta_2` are normalized in `create_laser`, so `a0` alone sets the strength. LG `Npn` normalization is
  deliberately custom.
- **Known bugs, revisit together:** `laser_delay` has no effect (overwritten with `wing_sigma_cutoff * wing_sigma`),
  `tau_0_traj` is hardcoded `0.0`, and the heatmap snapshot time depends on both.
- `export_field_heatmap_z0`'s window is `±beam_cylinder_radius` and ignores `beam_center_x/y`.

## Config file format

Flat `key value [unit]` lines; units via `IoUtils::convert_unit_to_number` (check it before adding one).
- **An unrecognized unit only warns** (assumes `1.0`). **A blank value drops the key** (`std::out_of_range`).
  Unknown keys are ignored. To see what a run actually used, read its `run_log.txt` (it also gives SI values).
- All keys are required except `require_clean_build` and `beam_distribution`. `w0` units need `laser_type=laguerre_gauss`.
- `run_log.txt` re-reads keys the way the factories do (`read_scaled`), not via getters on the physics classes.

## Physics notes worth knowing before comparing results

- **Rectangular and spherical screens differ for real** (exact `R`): they agree only when `N_F = a²/(Dλ) << 1`.
  Compare `Re/Im(F01)`, not `|F01|`.
- The Python reference is `~/Dropbox/work/bin/python/Superradiant_Thomson` on the laptops (not on the servers;
  ask the user). It has caught several real bugs.

Full history and rationale: `theory/dev_notes.md`.
