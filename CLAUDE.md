# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

Pipeline (`src/app/main.cpp`): build a config-driven laser and detector, generate the beam, integrate up to 10
electrons' trajectories for `electron.dat`, compute the form factor `b(q)` (diagnostic), run
`Simulation::run_simulation` over the whole beam, and write `.dat` files plus `run_log.txt` to
`~/<output_folder>/YYYYMMDD_HHMMSS/` (`output_folder` is relative to `$HOME`). `py_scripts/` plots the outputs.

**Open issues and the status of every document:** `theory/README.md`. Read it before starting work; it also
points to the lattice-stack "Resume here".

## Spec-driven workflow

Work arrives as a numbered task spec, `theory/specs/NN_<name>.md` (next free number in `theory/README.md`). When it is
done, `theory/implemented/NN_<name>.md`, with the same number, records the status, the commit, where the code is,
deviations from the spec, judgment calls (e.g. tolerances) and how it was checked.

- **A spec with a `Status: DRAFT, NOT APPROVED` notice is for discussion:** read it and the code it names, reply with
  your analysis, and wait. Remove the notice only when the user approves.
- Before implementing, read the real code and the earlier `implemented/` notes; reuse existing classes and
  conventions rather than creating competing copies. Implement only the current spec.
- If a needed API in earlier code is missing or broken, fix it there and record it in that spec's `implemented/` note.
- Afterwards write the `implemented/` note and update `theory/README.md` (status, open issues). Living references
  (`theory/reference/`) are kept in sync with the code; new log entries go to `theory/work/dev_notes.md`.
- Rejected specs are not numbered (archive or delete them).

## Machines

Two laptops share a checkout through Dropbox; two servers sync through git. Python is `.venv` on the laptops and
Anaconda on the servers. Compilers differ between machines.

**Slurm cluster** (`theory/cluster.md` for commands, limits, long-run recipe): head node `hpc-ctrl01`, compute nodes
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
  `__attribute__((optimize("no-tree-slp-vectorize")))` everywhere (`theory/reference/implementation_details.md` §6).

**Build provenance:** the binary embeds commit, dirty state, compiler, flags, host (`cmake/build_info.cmake`), printed
in `run_log.txt`; a dirty build writes `<run>/source.patch`. `require_clean_build true` refuses dirty builds. Git
state is read at build time: commit, then build.

Formatting: `clang-format -i` on touched files (Google style, 120 cols); not installed on the servers, so say so in
the commit message and format on a laptop. **No test suite**, no lint command.

### Configs and test inputs

**Iterating on code changes:** copy a config to your scratch space (don't edit `config/config.cfg` or a reference
unless the task is about it) and shrink it (`beam_particle_count`, screen points; a lattice stack shrinks only through
`beam_stack_sheets` or the spacing). `config/` holds:
- `reference/`: one config per beam situation (`random_2d`, `lattice_2d`, `random_3d`, `random_3d_modulated`,
  `lattice_stack_uniform`, `lattice_stack_modulated`, `spectrum_lattice_2d`, `spectrum_line_lattice_2d`), identical
  except for the beam section (and the spectra's screens): **a change to the shared base goes into all eight.**
  Production size; costs in `reference/README.md`. The user copies one to `config/work/` (ignored by git) to work.
- `config.cfg`: the default working config (a copy of `reference/random_2d.cfg`); scripts read its `output_folder`.
- `config_cross_check.cfg`: matched to the Python reference's `main.py` `INPUTS`; use it for every C++-vs-Python
  comparison and update it when `INPUTS` change.

Elsewhere:
- `reports/<name>/`: a report with every experiment that fed it (configs, batch and comparison scripts). Frozen configs
  are moved, never edited. `reports/2026-10-29_progress/`: the progress report; its `README.md` maps every quoted
  result to a run and a frozen config. Runs in `~/output/progress_report_20261029/`. If a quoted result changes,
  update the report.
- `checks/` (index: `checks/README.md`): reusable checks, all writing to `~/output_tests/...`: the Thomson dipole,
  one-electron configs, `nt_scan` (`trajectory_NT` convergence, ~3 min) and the lattice-stack checks. The latter
  **reuse finished runs in `~/output_tests/lattice_stack/`: delete those after changing the solver.**

### Verification (stands in for tests; use a small copied config)

- **`simplified` vs `direct`** (`radiation_formula`), comparing totals: the main check. They agree to `O(dτ²)` (4x per
  doubling of `trajectory_NT`) down to a 1 λ screen; a `1/(kR)` mismatch means a short-range bug.
- **Thomson dipole** (the only absolute check, under a second): run `checks/thomson_dipole/thomson_dipole.cfg`, then
  `python3 checks/thomson_dipole/thomson_dipole_check.py <run>`; agrees to ~1e-4.
- **Rotation covariance:** rotating laser, beam momentum and forward detector by `R` gives `F → Λ F Λᵀ` (6 digits).
- **Spherical screen:** `B` exactly ⊥ `n0`, `E` not (`|E_r|/|E_θ|` ~ 0.1), via `plot_spherical_components.py`.
- **Byte-identity** for intended no-ops: diff the `.dat` outputs with the same config and `random_seed`.
- **Python reference** `~/Dropbox/work/bin/python/Superradiant_Thomson` (laptops only; on servers ask the user) with
  `config_cross_check.cfg`. It has caught several real bugs. It keeps the pre-spec-14 arrangement, so compare only in
  the old geometry (`laser_delay = cutoff·σ`, as that config has).

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
python3 py_scripts/radiation/plot_field_carrier_removed.py [<range> [<mu> <nu>]] [run_folder]   # F e^{-ik_N(R-Z)}; all components by default; flat detectors only
python3 py_scripts/radiation/plot_point_spectrum{,_all,_observables}.py           # dense_frequency_spectrum, one point
python3 py_scripts/radiation/plot_spectrum_map.py [run_folder] [--range R]         # dense spectrum on a line: (theta, omega) maps
python3 py_scripts/radiation/plot_spherical_components.py <long|short|boundary> <E|B> <r|theta|phi>
python3 py_scripts/radiation/plot_observables.py [--incident]   # rect./circular facing ±Oz (--incident: any)
python3 py_scripts/radiation/lg_decomposition.py [--incident] ...   # options: theory/reference/lg_decomposition_reference.md
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
- **Beam along its axis** (full rules: `theory/implemented/09_lattice_stack.md`, last two sections): `random` is uniform
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
  `theory/specs/03_lg_beam_radius_table.md`: keep in sync.
- `CircularDetector`'s radial grid is **equal-area**. A detector with grid counts `1` is one exact point.

### Radiation calculation (the hot path)

Modes, optimizations, precision, timings: `theory/reference/implementation_details.md` (**keep in sync with
`radiation.cpp`**).
Formulas: `theory/physics/FT_faraday_tensor_v2.md`; long-distance modes:
`theory/specs/01_long_distance_formulas.md`; OAM: `theory/physics/angular_momentum_numerics.md`.

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
- `dense_frequency_spectrum=true`: `N_omega` points at one screen point or on a **line** (rectangular `Nx` or `Ny` 1,
  circular `N_phi` 1; spec 15); a 2D grid keeps only its first point. Memory ~ points × `N_omega` × threads: estimated
  on stdout and in `run_log.txt`, refused above 8 GB of buffers. Bounds in `omega_laser` or
  `first_harmonic_frequency` only. Internally frequencies are `omega/c`; the output `omega` column is `omega/omega_1`.
  Map only φ-independent quantities on a line (dP_z/dω, |E±|), not Cartesian components.

### Laser

- `get_faraday_tensor` builds on the pure-virtual `complex_amplitude` (`{amplitude, d/dx, d/dy}`); `Ez`/`Bz` are
  first-order paraxial (`div E` residual up to ~1e-2 is expected). Carrier sign `exp(-iφ)`, the `+i/k` coefficient,
  Gouy sign and curvature sign (`+kρ²/2R`) must flip together.
- `zeta_1/zeta_2` are normalized in `create_laser`, so `a0` alone sets the strength. LG `Npn` normalization is
  deliberately custom.
- **Collision timing** (spec 14): `laser_delay` = the time the plateau start (the peak, with no plateau) crosses the
  origin, where the LG focus is; beam positions are free-flight positions at t = 0, so `laser_delay 0` +
  `beam_center_z 0` collide at the focus. Each electron starts at the pulse front on its own free-flight line
  (`Particle::free_flight_crossing`); `b(q)`, `electron_beam_scatter.dat` and the axial diagnostics use
  `get_reference_position()` (t = 0), never `get_initial_position()` (the start). `laser_delay = cutoff·σ` reproduces
  the old off-focus geometry. `run_log.txt` prints the collision point.
- **⚠ Collision geometry:** the start time, pull-back and collision point are derived for φ = k(ct − n̂·x), focus at
  the origin, free flight before the pulse, head-on baseline. **Changing the collision geometry** (co-propagating
  or crossed beams, a moved focus, a non-φ envelope) **means re-deriving spec 14 §3 and re-running every check**;
  1 − n̂·β → 0 for co-propagation makes the start run away.
- The heatmap window ignores `beam_center_x/y`.

## Config file format

Flat `key value [unit]` lines; units via `IoUtils::convert_unit_to_number` (check it before adding one).
- **An unrecognized unit only warns** (assumes `1.0`). **A blank value drops the key.** Unknown keys are ignored. What
  a run actually used (with SI values) is in its `run_log.txt`, which re-reads keys the way the factories do.
- All keys are required except `require_clean_build`, `beam_distribution`, `beam_lattice_spacing`, `beam_stack_*` and
  `beam_modulation*`. `w0` units need `laser_type=laguerre_gauss`.

## Physics notes

- **Rectangular and spherical screens differ for real** (exact `R`): they agree only when `N_F = a²/(Dλ) << 1`.
  Compare `Re/Im(F01)`, not `|F01|`.

History and rationale: `theory/work/dev_notes_2026-07-28_to_2026-11-23.md` (frozen); new entries in
`theory/work/dev_notes.md`.
