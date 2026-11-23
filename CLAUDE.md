# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

Pipeline (`src/app/main.cpp`): build a config-driven laser and detector, generate the beam, integrate up to 10
electrons' trajectories for `electron.dat`, compute the form factor `b(q)` (diagnostic), run
`Simulation::run_simulation` over the whole beam, and write `.dat` files plus `run_log.txt` to
`~/<output_folder>/YYYYMMDD_HHMMSS/` (`output_folder` is relative to `$HOME`). `py_scripts/` plots the outputs.

**Open issues** (detail in `theory/dev_notes.md` unless noted):
- **Old runs need their observables re-run:** OAM numbers before 2026-10-29 (centered-difference `∂φ`, ~5% low at
  N = 2) and absolute spectral densities / photon numbers before 2026-11-21 (`2π` too small; ratios unaffected,
  `theory/spectral_density_2pi_fix.md`). Re-run `plot_observables.py` / `plot_point_spectrum_observables.py`
  before quoting them.
- `b(q)` (`Simulation::compute_form_factor`, `theory/thomson_form_factor_diagnostic.md`): harmonic 1 along the
  detector normal, `q = k_s − (K_s·P / K_L·P) k_L` with undressed mean `P`; a translated electron's field gets
  `e^{-i q·r}` and `q_N = N q_1` (not `N k_1 + k_L`).
- `LaguerreGaussLaser` matches the standard LG mode numerically but isn't physics-reviewed (Allen 1992; Siegman
  ch. 17); `PhysUtils::non_linear_Thomson_formula` is flagged in-code for review.
- The flat-top envelope leaks N = 1 into N = 2 (`laser_flat_duration 0` removes it); LG on-axis field is √2 × the
  plane wave's at the same `a0`: `theory/future_dev/pulse_envelope_and_lg_normalization.md`. Unacted plans:
  `theory/future_dev/README.md`.
- Ponderomotive dressing acts as a lens, so `N ≥ 2` isn't a single LG mode: `theory/future_dev/ponderomotive_lens.md`.
- **Lattice stack:** Phases 1–3a done; the beam rework (2026-11-09) and Gaussian stack (2026-11-10) await the user's
  review; Phase 3b next. State and "Resume here": `theory/lattice_stack-work.md`.

## Machines

Two laptops share a checkout through Dropbox; two servers sync through git. Python is `.venv` on the laptops and
Anaconda on the servers. Compilers differ between machines.

**Slurm cluster** (`theory/srun.txt` for commands, limits, long-run recipe): head node `hpc-ctrl01`, compute nodes
`hpc-node[01-20]` (256 CPUs, GCC 13.3), shared home.
- **Always `--exclude=hpc-ctrl01,hpc-db01`**: a job on the head node looks like the login shell (check
  `squeue -u $USER`, not `hostname`).
- Compute nodes have only CMake **3.21**: don't raise `cmake_minimum_required` or the presets' minimum.
- Build on the compute node that runs the job (`-march=native`); Claude's shell may be on the head node while the
  user's terminal is in a job, so check `hostname` before building. `num_threads 0` follows the job's allocation.
- Long runs go in `tmux` on the head node; never `exit` that tmux. Unattended: `srun ... bash -c "cd ... && <build>
  && <run>"` so the job ends with the run (a `--pty bash` job holds the node).
- `radiation_field.dat` appears only at the end. `run_log.txt` is written first (`STATUS: IN PROGRESS`) and rewritten
  at the end, so a stopped run still supports the non-radiation scripts.

## Commands

```
cmake -B build/ && cmake --build build/                    # build (out-of-source enforced)
cmake --preset release && cmake --build --preset release  # same via preset; `debug` preset opts out of Release
python3 py_scripts/compile_and_run.py <config_file>        # configure + build + run (config needed in practice)
./bin/coherent_thomson_solver config/config.cfg            # run (config argument required)
```

**Always built as `Release`** (`-O3 -march=native`, C++20) unless `-DCOHERENT_THOMSON_ALLOW_NON_RELEASE=ON` (set by
the `debug` preset; stays cached until `=OFF`). Debug is ~12x slower: check `Build type:` before trusting a timing.
The `debug` preset uses `build-debug/`; non-Release binaries go to `bin/<type>/`, so `bin/coherent_thomson_solver`
is always Release.

**Build traps, check these first:**
- **Incremental rebuilds have missed header dependencies** (Dropbox may sync `build/`, `bin/`). After editing a
  widely-included header, or if a build acts strangely, use `--clean-first`.
- **Every Release build writes the same `bin/coherent_thomson_solver`** and a no-op build doesn't relink: to compare
  builds (e.g. two commits), `--clean-first` and copy each binary out.
- **GCC 15.2's SLP vectorizer miscompiles the long-distance radiation functions.** Keep their
  `__attribute__((optimize("no-tree-slp-vectorize")))` everywhere (`theory/implementation_details.md` §6).

**Build provenance:** the binary embeds commit, dirty state, compiler, flags, host (`cmake/build_info.cmake`), printed
in `run_log.txt`; a dirty build writes `<run>/source.patch`. `require_clean_build true` refuses dirty builds. Git
state is read at build time: commit, then build.

Formatting: `clang-format -i` on touched files (Google style, 120 cols); not installed on the servers, so say so in
the commit message and format on a laptop. **No test suite**, no lint command.

### Configs and test inputs

**Iterating on code changes:** copy a config (don't edit `config/config.cfg` unless the task is about it) and drop
`beam_particle_count` to e.g. `50`. `config/` holds only generic configs:
- `config.cfg`: default random beam, current experiments. `config_lattice.cfg`: the same plus `polar_lattice`; keep
  its other keys in step with `config.cfg`.
- `config_lattice_stack.cfg` (run A, Gaussian pulse, `a0 = 0.1`) and `config_modulation_stack.cfg` (Gaussian stack,
  64 × 64 screen): current stack experiments, not kept in step. Estimate cost (sheets × sites) before running.
- `config_spectrum.cfg`: one-point dense spectrum. `config_cross_check.cfg`: matched to the Python reference's
  `main.py` `INPUTS`; use it for every C++-vs-Python comparison and update it when `INPUTS` change.

Frozen experiments:
- `theory/scientific_report/`: the 2026-10-29 progress report; its `README.md` maps every quoted result to a run and a
  frozen config. Runs in `~/output/progress_report_20261029/`. If a quoted result changes, update the report.
- `theory/lattice_stack_tests/`: lattice-stack checks (results in `*_results.txt`). They **reuse finished runs in
  `~/output_tests/lattice_stack/`: delete those after changing the solver.** `rework_checks.py` and
  `gaussian_stack_checks.py` take `<old_solver> <new_solver>`; `phase3a_noop_and_rejects.py`/`phase3a_switch.py` need
  pre-rework binaries.
- `theory/future_dev/tests/`: one-electron configs, the Thomson-dipole benchmark, `nt_scan.py` (`trajectory_NT`
  convergence, ~3 min), `bunch_length/`, `ponderomotive_lens/`. All write to `~/output_tests/...`.

### Verification (stands in for tests; use a small copied config)

- **`simplified` vs `direct`** (`radiation_formula`), comparing totals: the main check. They agree to `O(dτ²)` (4x per
  doubling of `trajectory_NT`) down to a 1 λ screen; a `1/(kR)` mismatch means a short-range bug.
- **Thomson dipole** (the only absolute check, under a second): run `theory/future_dev/tests/thomson_dipole.cfg`, then
  `python3 theory/future_dev/tests/thomson_dipole_check.py <run>`; agrees to ~1e-4.
- **Rotation covariance:** rotating laser, beam momentum and forward detector by `R` gives `F → Λ F Λᵀ` (6 digits).
- **Spherical screen:** `B` exactly ⊥ `n0`, `E` not (`|E_r|/|E_θ|` ~ 0.1), via `plot_spherical_components.py`.
- **Byte-identity** for intended no-ops: diff the `.dat` outputs with the same config and `random_seed`.
- **Python reference** `~/Dropbox/work/bin/python/Superradiant_Thomson` (laptops only; on servers ask the user) with
  `config_cross_check.cfg`. It has caught several real bugs.

### Plotting

Every script plots the latest run under `config/config.cfg`'s `output_folder`, whatever config made it; for runs
elsewhere (e.g. `~/output_tests/`) pass the **run folder** (not a `.dat`) as the trailing argument. PNGs go to
`<run>/png_folder/<module>/` (`MODULE_NAME` from the script's folder; don't hardcode it). `plot_observables.py`,
`plot_point_spectrum_observables.py` and `lg_decomposition.py` also **write** sections into `run_log.txt`.
Conventions: `py_scripts/PLOTTING_NOTES.md`.

```
python3 py_scripts/{laser,detector,particle}/<script>.py   # some need plot_*=true flags in the config
python3 py_scripts/radiation/plot_field.py <long|short|boundary|total> <mu> <nu> [--incident] [run_folder]
python3 py_scripts/radiation/plot_all_components.py [--incident]
python3 py_scripts/radiation/plot_processed.py radial_phase <long|short|boundary|total> <mu> <nu>   # or `radial_phase all`; flat detectors only
python3 py_scripts/radiation/plot_point_spectrum{,_all,_observables}.py           # dense_frequency_spectrum runs
python3 py_scripts/radiation/plot_spherical_components.py <long|short|boundary> <E|B> <r|theta|phi>
python3 py_scripts/radiation/plot_observables.py [--incident]   # rect./circular facing ±Oz (--incident: any)
python3 py_scripts/radiation/lg_decomposition.py [--incident] ...   # options: theory/lg_decomposition_reference.md
```
`--incident` uses `incident_field.dat` (analytic laser field; rectangular/circular detectors only). `mu == nu` is
rejected; `mu > nu` uses `F^{nu mu} = -F^{mu nu}`. `E_i = c·F^{i0}`, `B = (F^{32}, F^{13}, F^{21})`.

## Architecture

### Layout

- `main.cpp` only orchestrates; all logic is in the static library `coherent_thomson_core` (`src/core/`, headers in
  `src/core/include/<module>/`), namespaces `Core::<CamelCase of folder>`. `PhysUtils::AtomicUnits`:
  `c = 137.035999084`, `q_0 = -1` (signed charge), `e_0 = +1`, `epsilon_0 = 1/(4π)`.
- **Class + factory** (`laser/`, `detector/`): `create_<thing>(const ConfigMap&, ...)` in `<thing>_factory.cpp` holds
  *all* config parsing and units (`particle/`: `generate_cylinder_beam`). New component: class + factory, add the
  `.cpp` to `src/core/CMakeLists.txt`, wire it into `main.cpp`.
- `py_scripts/<module>/` mirrors the C++ modules; scripts import `utils.*` via `sys.path`.
- **C++/Python copies, keep in sync:** constants (`phys_utils.hpp`), `convert_unit_to_number` (Python in
  `faraday_frame_utils.py` *raises* on unknown units where C++ warns), `rotation_matrix_from_direction`.

### Frames and geometry

- **Baseline** (`config/config.cfg`): electrons and screen along +Oz, laser along −Oz (`laser_direction_theta 1.0 pi`).
- **Only the laser moves** with `laser_direction_theta/phi`; beam and detector keep their lab geometry, so a tilt
  changes the collision geometry. `epsilon_1/2 = R x̂, R ŷ`, `R = MathUtils::rotation_matrix_from_direction(n)`
  (proper: helicity and LG `l` sign are relative to `n`; for `θ=π`, `epsilon_2 = −ŷ`). It snaps to ±Oz only within
  ~1e-12 rad: don't loosen it.
- `laser_field*.dat` and `incident_field.dat` are in the **laser's frame** `(epsilon_1, epsilon_2, n)`;
  `incident_field.dat` is at `z = Z` past the waist (`Radiation::incident_screen`, plus `incident_config.cfg`).
  `radiation_field.dat` is lab-frame components on the detector's local grid (mirrored `(x, −y)` facing −Oz).
- The beam axis follows the mean momentum; `beam_center_x/y/z` apply **before** that rotation.
- **Beam along its axis** (full rules: `theory/lattice-stack-implemented.md`, last two sections): `random` is uniform
  in `H = beam_cylinder_height` with Gaussian end smoothing `beam_sigma_z` (`0` is byte-identical to the plain
  cylinder); `polar_lattice` puts `beam_stack_sheets` sheets over `L`. With `beam_stack_sigma_z > 0` the Gaussian is in
  the sheet **charges**, not positions (no `1/√M` noise). Each type ignores the other's keys, so switching is one edit.
  Weights scale the field and form factor, never the motion. `beam_modulation true` displaces every electron
  deterministically. The solver warns on too-coarse spacing and too few sheets per period (aliasing).
  Unit `lambda_1 = 2π/|q_axial|` (axial part of `b(q)`'s `q`; needs laser, detector normal and beam momentum) is
  accepted only by `beam_cylinder_height`, `beam_sigma_z`, `beam_stack_sigma_z`, `beam_modulation_period` and
  **throws** on any other key. Parsing keeps it in units of `lambda_1`; `Particle::resolve_beam_lengths` converts it
  (so the detector is created before the beam).
  Old stack/modulation keys are refused. Don't judge an `N = 2` stack at one on-axis point (~2.6e4 below the screen).
- A fixed `random_seed` with fewer electrons gives the first electrons of a larger beam: use different seeds for
  independent samples.
- Negative `beam_cylinder_radius` = `IoUtils::lg_beam_radius_in_w0(p, l) · w0` (LG only); its table is copied from
  `theory/table of beam electron radius in units of w_0 for different p,l.md`: keep in sync.
- `CircularDetector`'s radial grid is **equal-area**. A detector with grid counts `1` is one exact point.

### Radiation calculation (the hot path)

Modes, optimizations, precision, timings: `theory/implementation_details.md` (**keep in sync with `radiation.cpp`**).
Formulas: `theory/FT_Faraday_tensor-direct_and_simplified_forms-v2.md`; long-distance modes:
`theory/long_distance_direct_simplified_coding_guide.md`; OAM: `theory/numerical_calculation_of_angular_momentum.md`.

- **Three terms** (`LR_`/`SR_`/`BR_` columns): long-range, short-range, boundary. **Only the sum is physical**; the
  split differs in LD modes, so compare only totals.
- `radiation_formula`: `simplified`, `direct` (`F_b ≡ 0`), `long_distance_{simplified,direct}` (no `F_s`, error
  ~`1/D`), `long_distance_{simplified,direct}_approx` (stepped phase, fastest, ≤1e-7 from unstepped).
- `general_factor` uses the **signed** `q_0`; direct `F_s` has an explicit `c²` to compensate. Trapezoid weights apply
  to `F_l`/`F_s` only; `F_b` is an unweighted endpoint evaluation.
- **Load-bearing:** prefactors stay in `radiation.cpp`'s anonymous namespace; `radiation_phase_argument` stays
  **frequency-independent** (frequency recurrence over evenly spaced `frequencies_list`); screen point outer, tau
  inner; in LD modes the constant phase is applied at reconstruction; stepped loops must stay vectorizable
  (`-fopt-info-vec`).
- **Phases carry only `e^{ik(R - |X|)}`**; `run_simulation` applies `e^{ik|X|}`. Never put `R` or `|x_0|` back into a
  phase: rounding `kR ~ 1e9–1e17` per electron destroys coherence (`implementation_details.md` §4.5).
- **Threading:** contiguous electron chunks, one `std::thread` each, per-thread `PackedRadiationField` reduced after
  `join`; `compute_radiation` must not write shared state. "Time per electron" is CPU time.
- **Trajectories are freed per electron** (`Electron::release_trajectory()`): nothing may read a trajectory after
  `run_simulation`; `electron.dat` is written before it.

### Frequencies

- `frequencies_list[i] = non_linear_Thomson_formula(k1, q, n2, N_harmonics_min + i)`, `q` dressed with
  `<a²> = a0²/2` for every polarization. **Don't change it to `a0²`** (a past bug).
- `dense_frequency_spectrum=true`: `N_omega` points at one screen point (bounds in `omega_laser` or
  `first_harmonic_frequency` only). Internally frequencies are `omega/c`; the output `omega` column is `omega/omega_1`.

### Laser

- `get_faraday_tensor` builds on the pure-virtual `complex_amplitude` (`{amplitude, d/dx, d/dy}`); `Ez`/`Bz` are
  first-order paraxial (`div E` residual up to ~1e-2 is expected). Carrier sign `exp(-iφ)`, the `+i/k` coefficient,
  Gouy sign and curvature sign (`+kρ²/2R`) must flip together.
- `zeta_1/zeta_2` are normalized in `create_laser`, so `a0` alone sets the strength. LG `Npn` normalization is
  deliberately custom.
- **Known bugs, revisit together:** `laser_delay` has no effect (overwritten with `wing_sigma_cutoff * wing_sigma`),
  `tau_0_traj` is hardcoded `0.0`, and the heatmap snapshot time depends on both. The heatmap window ignores
  `beam_center_x/y`.

## Config file format

Flat `key value [unit]` lines; units via `IoUtils::convert_unit_to_number` (check it before adding one).
- **An unrecognized unit only warns** (assumes `1.0`). **A blank value drops the key.** Unknown keys are ignored. What
  a run actually used (with SI values) is in its `run_log.txt`, which re-reads keys the way the factories do.
- All keys are required except `require_clean_build`, `beam_distribution`, `beam_lattice_spacing`, `beam_stack_*` and
  `beam_modulation*`. `w0` units need `laser_type=laguerre_gauss`.

## Physics notes

- **Rectangular and spherical screens differ for real** (exact `R`): they agree only when `N_F = a²/(Dλ) << 1`.
  Compare `Re/Im(F01)`, not `|F01|`.

Full history and rationale: `theory/dev_notes.md`.
