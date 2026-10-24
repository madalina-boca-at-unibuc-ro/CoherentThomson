# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

Pipeline (`src/app/main.cpp`): build a config-driven laser pulse and detector screen, generate an electron beam,
integrate the first up-to-10 electrons' trajectories for `electron.dat`, compute the initial-position form factor
`b(q)` (diagnostic only; see Open issues), run the coherent radiation calculation
(`Simulation::run_simulation`) over the whole beam, and export `.dat` files plus `run_log.txt` into a per-run
directory `~/<output_folder>/YYYYMMDD_HHMMSS/` (`output_folder` is relative to `$HOME`, so one config works on
every machine). Python scripts in `py_scripts/` plot those outputs.

**Open issues** (detail in `theory/dev_notes.md`):
- The single-electron Thomson-dipole benchmark (electron at rest, full-4π spherical detector, weak field →
  `1+cos²θ` for circular polarization) has **not** been reproduced yet. Don't trust absolute intensities/angular
  patterns for more complex setups until it is.
- `plot_observables.py`'s `Lambda_zz/P_z` now has the right sign everywhere (the old backward-screen sign flip was
  a local-grid vs. lab-component mix, fixed with `lz_sign`). Theory column: `L_N/omega_N` with
  `L_N = N(m+σ)−σ` for the emitted harmonics (`σ = 2 Im(ζ₁*ζ₂)`, 0 for linear; `m`, `σ` times `sign(n_laser,z)`);
  `m/omega` with `m = laser_lg_l` for `--incident`. Still ~5% low at N = 2 and not settled; measurements, and why
  flipping only the screen at `γ ≠ 0` isn't a clean test, are in dev notes.
- The initial form factor `b(q)` (`Simulation::compute_form_factor`; spec `theory/thomson_form_factor_diagnostic.md`)
  is one scalar at harmonic 1 along the detector normal, using the generated beam's mean momentum (not the
  configured `average_p`), printed to stdout and `run_log.txt`. It uses the dressed `omega_s` with
  `q = k_s − k_L`, which may be inconsistent away from head-on collisions (dev notes). Don't rely on `b` in tilted
  geometries until that's settled.
- `LaguerreGaussLaser`'s mode formula isn't physics-reviewed against the literature (Allen et al. PRA 45, 8185
  (1992); Siegman ch. 17), though its derivatives match finite differences. `PhysUtils::non_linear_Thomson_formula`
  is flagged in-code as still needing review.

## Machines

Two laptops share a checkout through Dropbox (under `~/Dropbox`); two servers are synced through git. Python is
`.venv` on the laptops and Anaconda on the servers (`pyrefly.toml`/`pyrightconfig.json` point at `.venv`, so editor
type-checking won't find the interpreter on the servers). Compilers differ between machines.

**Slurm cluster** (one of the servers): head node `hpc-ctrl01`, compute nodes `hpc-node[01-20]` (256 CPUs each,
GCC 13.3), with the home directory shared between them. Runs go on a compute node, e.g. interactively with
`srun --partition=debug --time=01:00:00 --cpus-per-task=64 --pty bash` (`debug` caps a session at 1 h; other
`srun` variants, e.g. `--exclusive`, are in `theory/srun.txt`).
- Compute nodes have no system CMake, only Anaconda's **3.21**, which is why `cmake_minimum_required` (and the
  presets' `cmakeMinimumRequired`) is 3.21. Don't raise it.
- Because of `-march=native`, build on the compute node that runs the job: `compile_and_run.py` inside the job does
  that. Claude's shell may be on the head node even while the user's terminal is in a job (check `hostname`), so
  don't build there for a compute-node run.
- `num_threads 0` follows the job's allocation (the affinity mask), not the node's 256 CPUs.

## Commands

```
cmake -B build/ && cmake --build build/                    # build (out-of-source enforced)
cmake --preset release && cmake --build --preset release  # same via preset; `debug` preset opts out of Release
python3 py_scripts/compile_and_run.py <config_file>        # configure + build + run (arg is optional to the
                                                           # script, but without it the solver just prints usage)
./bin/coherent_thomson_solver config/config.cfg            # run (config argument required)
```

**Always built as `Release`** (`-O3 -march=native -mtune=native`, C++20, `-Wall -Wextra -Wpedantic`).
`CMakeLists.txt` forces it over any cached or IDE-supplied `CMAKE_BUILD_TYPE`. Another type needs
`-DCOHERENT_THOMSON_ALLOW_NON_RELEASE=ON` (the `debug` preset sets it), which stays cached until you pass `=OFF`.
Debug is ~12–13x slower, so check the configure step's `Build type: ...` line before trusting any timing.
`-march=native` means a binary must be built on the machine that runs it (on the cluster: see Machines).

**Build traps, check these first:**
- **Incremental rebuilds have missed header dependencies** (not root-caused). After editing a widely-included
  header (e.g. `phys_utils.hpp`), rebuild with `cmake --build build/ --clean-first`.
- **Every build writes the same `bin/coherent_thomson_solver`,** and a no-op build doesn't relink. To compare
  builds, use `--clean-first` and copy each binary out right after building. The `release` and `debug` presets
  share one `build/` directory and CMake cache. Both export `build/compile_commands.json`.
- **On the laptops, Dropbox syncs `build/` and `bin/`** (gitignored, not Dropbox-ignored), so the other laptop's
  binaries/objects and changed timestamps can arrive. This may explain the missed rebuilds. If a build acts
  strangely, `--clean-first` on this machine before debugging anything else.
- **GCC 15.2's SLP vectorizer miscompiles the long-distance radiation functions** (wrong sums at every frequency
  after the first). `compute_radiation_long_distance` and `compute_radiation_long_distance_stepped` (in
  `radiation.cpp`) both carry `__attribute__((optimize("no-tree-slp-vectorize")))`, guarded only by GCC-not-Clang.
  GCC 13.3 doesn't have the bug, but the workaround is **kept on purpose on every machine**: don't remove it or
  make it depend on the GCC version.
  Details: `theory/implementation_details.md` §6.

**Build provenance (reproducibility):** every build runs `cmake/build_info.cmake`, which embeds the git commit,
dirty state (untracked files count), compiler, flags and build host into the binary (`Core::BuildInfo`,
`build_info/`; the generated `build/src/core/build_info_data.cpp` is rewritten only when it changes, so no-op
builds still don't relink). `run_log.txt`'s `Build` section prints it, with the run host. A dirty build also
embeds the uncommitted changes and writes them to `<run>/source.patch`: `git apply` it on a checkout of the
recorded commit to get the built tree back. Config key `require_clean_build true` (optional, default `false`)
refuses to run a dirty build before any output is written. Use it for publication runs. The git state is read when
you build, not when you run, so commit first, then build.

**Iterating on code changes:** copy a config (don't edit `config/config.cfg` in place unless the task is about the
default config) and drop `beam_particle_count` to e.g. `50`; run time scales with electron count. Configs:
- `config/config.cfg`: default, for current experiments; its parameters change freely.
- `config/config_cross_check.cfg`: parameter-matched to the Python reference's `main.py` `INPUTS` (mapping and
  date in its header). Use it for every C++-vs-Python comparison, and update it when those `INPUTS` change.
- `config/config_initial_momentum.cfg`.

Formatting: `clang-format -i` on touched files (Google style, 120 cols, 2-space indent). **There is no test suite**
and no enforced lint command.

**Verification** (what stands in for tests; run on a small copied config):
- **`simplified` vs `direct`** (`radiation_formula`), comparing totals: the main correctness check. They agree to the
  trapezoid's `O(dτ²)` (4x per doubling of `trajectory_NT`) down to a 1 λ screen. A mismatch of order `1/(kR)`
  means a short-range bug. Simplified `F03` near the axis and weak harmonics need a finer `trajectory_NT`.
- **Rotation covariance:** rotating laser, beam momentum and a forward detector by the same `R` must give
  `F → Λ F Λᵀ` (holds to the 6-digit output precision).
- **Spherical screen:** per contribution, `B` is exactly ⊥ `n0` but `E` is not (`|B_r|/|B_θ|` ~ 1e-5,
  `|E_r|/|E_θ|` ~ 0.1), via `plot_spherical_components.py`.
- **Byte-identity** for a change that should be a no-op in some limit (e.g. a refactor, or `beam_sigma_z 0`): diff
  the `.dat` outputs before and after with the same config and `random_seed`.
- **Python reference** (Physics notes below) with `config/config_cross_check.cfg`.

### Plotting

With no arguments, every script plots the latest run. The optional trailing argument is a **run folder**, not a
`.dat` path. PNGs go to `<run>/png_folder/<module>/`. Each script sets `MODULE_NAME` from its own folder name;
don't hardcode it. Rendering conventions and `plot_observables.py` internals: `py_scripts/PLOTTING_NOTES.md`.

```
python3 py_scripts/laser/plot_field.py
python3 py_scripts/laser/plot_heatmap_z0.py        # needs plot_field_heatmap=true
python3 py_scripts/detector/plot_stereographic.py
python3 py_scripts/detector/plot_scatter.py        # needs plot_detector_scatter=true
python3 py_scripts/particle/plot_trajectory.py
python3 py_scripts/particle/plot_beam_scatter.py   # needs plot_beam_scatter=true
python3 py_scripts/particle/plot_beam_density_z.py # analytic z-density from config; + histogram if plot_beam_scatter=true
python3 py_scripts/radiation/plot_field.py <long|short|boundary|total> <mu> <nu> [--incident] [run_folder]
python3 py_scripts/radiation/plot_all_components.py [--incident]   # plot_field over all 4 ranges x 6 (mu,nu)
python3 py_scripts/radiation/plot_processed.py radial_phase <long|short|boundary|total> <mu> <nu> [run_folder]
                                                   # or `radial_phase all`; F·e^{-ik_N(R-Z)}, flat detectors only
                                                   # (theory/post_processing.md); -> png_folder/radiation/processed/<processing>/
python3 py_scripts/radiation/plot_point_spectrum.py <long|short|boundary|total> <mu> <nu> [run_folder]
                                                   # dense_frequency_spectrum=true; -> png_folder/radiation/point_spectrum/
python3 py_scripts/radiation/plot_point_spectrum_all.py [run_folder]          # 4 ranges x 6 components + |F| overviews
python3 py_scripts/radiation/plot_point_spectrum_observables.py [run_folder]  # S_z, Sigma_zz, u, P_z + spin/photon
                                                   # ratios; local only (no L_z/Lambda_zz: they need a screen grid)
python3 py_scripts/radiation/plot_spherical_components.py <long|short|boundary> <E|B> <r|theta|phi>  # spherical only
python3 py_scripts/radiation/plot_observables.py [--incident] [run_folder]   # rectangular/circular facing ±Oz only
                                                                             # (--incident: any laser direction)
```
`--incident` analyzes `incident_field.dat` (the incident laser's analytic field, written only for
rectangular/circular detectors) instead of `radiation_field.dat`. For `mu > nu` the scripts use
`F^{nu mu} = -F^{mu nu}`; `mu == nu` is rejected. Field extraction: `E_i = c·F^{i0}`, `B = (F^{32}, F^{13}, F^{21})`.

## Architecture

### Layout

- `src/app/main.cpp` only orchestrates. All logic is in the static library `coherent_thomson_core` (`src/core/`),
  with public headers in `src/core/include/<module>/`. `math_utils/`, `phys_utils/`, `io_utils/` are header-only.
  `logging/` (`Core::Logging::write_run_log`) writes `run_log.txt`.
- Namespaces are `Core::<CamelCase of folder>` (`Core::Radiation`, `Core::PhysUtils`, ...). `PhysUtils::AtomicUnits`
  holds the constants: `c = 137.035999084`, `q_0 = -1` (signed charge), `e_0 = +1` (magnitude), `epsilon_0 = 1/(4π)`.
- **Class + factory pattern** (`laser/`, `detector/`): `<thing>_factory.cpp` defines `create_<thing>(const
  ConfigMap&, ...)`, and *all* config-key parsing and unit handling lives there. `particle/`'s factory is
  `generate_cylinder_beam`. `simulation/` is free functions. Each module with output has a `<module>_plotter`.
- **Adding a physics component:** class + factory, add the `.cpp` to `add_library(...)` in
  `src/core/CMakeLists.txt`, wire the factory into `main.cpp`.
- `py_scripts/<module>/` mirrors the C++ modules. Scripts run as plain files (not `python -m`): each inserts
  `py_scripts/` onto `sys.path` and imports `utils.*` (a namespace package, no `__init__.py`).
- **C++/Python duplication, keep in sync:** constants (`C_LIGHT`, `EPSILON_0`, ...) are copied in the scripts from
  `phys_utils.hpp`. `radiation/faraday_frame_utils.py::convert_unit_to_number` copies the C++
  `IoUtils::convert_unit_to_number`, but **raises** on an unknown unit where C++ warns. Add any new unit to both.

### Frames and geometry

- **Baseline orientation** (`config/config.cfg`): electrons and screen along +Oz, laser along −Oz (head-on,
  `laser_direction_theta 1.0 pi`). The forward screen's local frame is the lab frame, so its plots need no mirror.
- The laser propagates along `laser_direction_theta/phi` (required; `(0,0)` = +Oz), in the single lab frame used
  everywhere. **Only the laser moves**: the beam and detector keep their own lab-frame geometry, so a tilt is a real
  change of collision geometry.
  `epsilon_1/epsilon_2` (what `zeta_1/zeta_2` refer to) are `R x̂, R ŷ` with `R =
  MathUtils::rotation_matrix_from_direction(n)`, the detector's own Oz→n rotation. It's proper, so helicity and
  the sign of LG `l` are relative to `n` (for `θ=π`: `epsilon_1 = x̂`, `epsilon_2 = −ŷ`). `run_log.txt` prints all
  three axes. `k1` and the `d_tau_traj` Doppler scaling (`(p⁰ − n·p)/mc`) follow `n`.
- `laser_field.dat` and `laser_field_heatmap_z0.dat` are in the **laser's own frame**: the heatmap is the laser's
  transverse plane through the waist, and both files' `Ex..Bz` are components along `(epsilon_1, epsilon_2, n)`.
  `incident_field.dat` is in the laser's own frame too (grid on its `(epsilon_1, epsilon_2)` plane at the waist,
  `F^{mu nu}` rotated to `(epsilon_1, epsilon_2, n)`), so it always shows the laser unrotated, as if along +Oz
  (older runs: lab frame; dev notes).
  `radiation_field.dat` stays lab-frame components on the detector's local grid (mirrored `(x, −y)` on a screen
  facing −Oz).
- The detector has its own direction (`detector_direction_theta/phi`, default `(0,0)` = forward).
- The beam cylinder's axis follows the mean momentum (`average_px/py/pz`), not Oz. `beam_center_x/y/z` are applied
  **before** that rotation, so `beam_center_z` moves along the direction of motion.
- Electron positions are independent random draws (no lattice). `beam_sigma_z` (required) adds an unclipped Gaussian
  displacement to the uniform z, smoothing the cylinder's sharp ends whose `sinc` form factor otherwise leaves a
  large coherent term at `N ~ 1e5–1e6` (`theory/gaussian_smoothed_electron_cylinder.md`). At `0` the Gaussian draw is
  skipped, so the beam (and every output) is byte-identical to the plain cylinder. With a fixed `random_seed`, a
  smaller `beam_particle_count` gives the first electrons of a larger one, so use different seeds for independent
  samples.
- A negative `beam_cylinder_radius` (value and unit otherwise ignored) sets the radius to
  `IoUtils::lg_beam_radius_in_w0(p, l) · w0`, just larger than the LG mode; it needs `laser_type=laguerre_gauss` and
  throws outside the tabulated `p = 0..4`, `|l| = 0..6`. The table is copied from
  `theory/table of beam electron radius in units of w_0 for different p,l.md` (same values as
  `table-with-R_beam.png`); keep them in sync. `run_log.txt` notes when it was used.
- `MathUtils::rotation_matrix_from_direction` (detector, laser `epsilon_1/2`, beam axis; Python copy in
  `faraday_frame_utils.py`, keep in sync) snaps to the exact ±Oz matrix only within ~1e-12 rad. Don't loosen it:
  small tilts are physical here (`θ_spot ~ laser tilt · ω/ω₁`).
- `CircularDetector`'s radial grid is **equal-area**: `r_i = sqrt(R_min² + i*(R_max²-R_min²)/(N_R-1))`. Anything
  reconstructing `r` from `i` must use this.
- Every detector collapses to one exact point when its grid counts are `1`.

### Radiation calculation (the hot path)

All six `radiation_formula` modes, their formulas, the optimizations, precision limits and timings are in
`theory/implementation_details.md`. **Keep it in sync when changing `radiation.cpp`.** Current formulas:
`theory/FT_Faraday_tensor-direct_and_simplified_forms-v2.md` (supersedes the unsuffixed v1). The long-distance
modes follow `theory/long_distance_direct_simplified_coding_guide.md` (reviewed against the exact code in
`long_distance_implementation_review.md`). `plot_observables.py`'s angular-momentum quantities are derived in
`theory/numerical_calculation_of_angular_momentum.md`.

- **Three terms**, exported as `LR_`/`SR_`/`BR_F<mu><nu>` columns: long-range `F_l`, short-range `F_s`, boundary
  `F_b`. **Only the sum is physical**; e.g. on-axis `F_l^{03}` and `F_b^{03}` are individually large and cancel.
- `radiation_formula` is **required** (no default; unknown names throw): `simplified`, `direct` (needs the stored
  4-acceleration; `F_b ≡ 0`), `long_distance_{simplified,direct}` (far screens only; no `F_s`; error ~`1/D`), and
  `long_distance_{simplified,direct}_approx` (stepped phase, fastest, matches unstepped LD to ≤1e-7). In the LD
  modes the `LR_`/`BR_` split differs from the exact one, so compare only totals.
- `general_factor = q_0 / (2π · 4πε₀c²)` uses the **signed** charge `q_0`. The direct `F_s` carries an explicit
  `c²` in `short_range_prefactor_direct` to compensate.
- **Quadrature:** trapezoidal weights apply to `F_l`/`F_s` only; `F_b` is an exact endpoint evaluation with **no**
  weight.
- **Keep these, they're load-bearing:** prefactors stay isolated in `radiation.cpp`'s anonymous namespace, and
  `radiation_phase_argument` must stay **frequency-independent** (the phase uses a frequency recurrence, valid
  because `frequencies_list` is evenly spaced). Loop order is screen point outer, tau inner. In the LD modes the
  constant phase `e^{ik(|x_0| - |X| + r^0(tau_m))}` is applied at reconstruction, **not** inside the tau loop
  (precision, §4.3). The stepped loops must stay vectorizable (check `-fopt-info-vec`, §6).
- **Every mode carries only `e^{ik(R - |X|)}`** (`distance_minus_reference`, computed as `(|r|² − 2X·r)/(R + |X|)`
  without cancellation); `run_simulation` applies `e^{ik|X|}` per (frequency, screen point), next to
  `general_factor`. Never put `R` or `|x_0|` itself back into a phase: per-electron rounding of `kR ~ 1e9–1e17`
  destroys coherence at far screens (grows like `γ⁴` at fixed Fresnel number). The exported field still has the
  full `e^{ikR}`. Details and measurements: `implementation_details.md` §4.5, `stabilize_phase_split_proposal.md`
  (supersedes `stabilize_phase_axis_aligned_detector.md`, an Oz-only plan: don't implement that one).
- The whole beam is generated upfront and held in memory (deliberate, for now).
- **Threading** (`simulation.cpp`): the beam is split into `num_threads` contiguous electron chunks, one
  `std::thread` each. Each thread accumulates into its own `PackedRadiationField` (6 packed complex components per
  point), and the per-thread fields are reduced to the full Faraday tensor after `join`. `compute_radiation` must not
  write any shared state. `num_threads` is a required key (`0` = every CPU in the process's affinity mask, which
  under Slurm is the job's allocation, not the whole node; it is silently capped at that count and at the electron
  count, so check the printed "Using N threads" line when timing). The reported
  "Time per electron" sums each thread's own elapsed time, so it is CPU time, not wall-clock time.

### Frequencies

- Default: `frequencies_list[i] = non_linear_Thomson_formula(k1, q, n2, N_harmonics_min + i)`, with `n2` the
  detector direction and `q` the **ponderomotively dressed** momentum with `<a²> = a0²/2` for every polarization.
  **Do not change it to `a0²`**; that was a real bug (dev notes).
- `dense_frequency_spectrum=true`: a linear scan of `N_omega` points for a single screen point; bounds accept only
  `omega_laser` or `first_harmonic_frequency` units. With more than one detector point, `main.cpp` warns and keeps
  point 0.
- Internally frequencies are `omega/c`. `radiation_field.dat`'s `omega` column is **`omega/omega_1`**;
  `run_log.txt` records `fundamental_frequency` in atomic units.

### Laser

- `get_faraday_tensor` is built on the pure-virtual `complex_amplitude` (`{amplitude, d/dx, d/dy}`). `Ez`/`Bz` are
  first-order paraxial: a `div E` residual of ~1e-6 (fundamental) to ~1e-2 (higher `p`/`|l|`) is expected. The
  carrier sign `exp(-iφ)` and the `+i/k` coefficient must flip together.
- Polarization `zeta_1`/`zeta_2` are complex and **normalized in `create_laser`**, so `a0` alone sets the strength.
- LG: radial profile via `hypergeometric_1F1_neg_int_a`; `Npn` is a deliberately custom normalization; azimuthal
  factor `(x + i·sign(l)·y)^|l|`; waist at the origin.
- Envelope wings are `exp(-Δφ²/wing_sigma²)` (not `/(2σ²)`, a fixed bug).
- **Known bugs/hardcodes, revisit together:** `laser_delay` has no effect (overwritten with
  `wing_sigma_cutoff * wing_sigma`), `tau_0_traj` is hardcoded to `0.0`, and the heatmap snapshot time is keyed off
  `wing_sigma_cutoff*wing_sigma` because of them.
- `export_field_heatmap_z0`'s window is `±beam_cylinder_radius` and ignores `beam_center_x/y`.

## Config file format

Flat `key value [unit]` lines. `IoUtils::convert_unit_to_number` resolves units (`lambda`, `pi`, `mc`,
`cycles_adim`, `omega_laser`, `w0`, `a.u.`, ...); check it before adding a unit.
- **An unrecognized unit only warns** and assumes `1.0`, so a typo fails silently.
- Keys are required (`get_required`) except `require_clean_build`.
- **A key with a blank value is dropped**, so reading it throws `std::out_of_range` (`IoUtils::get_required`
  names the missing key; a bare `config.at(...)` doesn't).
- `w0` is valid only with `laser_type=laguerre_gauss`.
- Unknown keys are ignored, so stale keys from old configs don't error.

`run_log.txt` re-reads config keys the way the factories do (`read_scaled`) rather than adding getters to the
physics classes. Beside each atomic-unit value it writes SI columns (lengths in mm, wavelength in nm, period and
duration in fs) from `PhysUtils::AtomicUnits::bohr_radius_m`/`atomic_time_unit_s` (2018 CODATA), so a run can be
read in physical units. To confirm what a run actually used (unit typos only warn), read its `run_log.txt` rather
than the config.

## Physics notes worth knowing before comparing results

- **Rectangular vs. spherical screens differ for real** (exact `R`, no far-field approximation). They agree only
  when the Fresnel number `N_F = a²/(Dλ) << 1`; growing `D` at a fixed *angular* window makes `N_F` larger. Compare
  `Re/Im(F01)`, not `|F01|`.
- The independent Python reference is `~/Dropbox/work/bin/python/Superradiant_Thomson` on the laptops (not on the
  servers; ask the user). Comparing against it has caught several real bugs.

Full investigation history, measurements, and rationale: `theory/dev_notes.md`.
