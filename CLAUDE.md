# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

Pipeline (`src/app/main.cpp`): build a config-driven laser pulse and detector screen, generate an electron beam,
integrate the first up-to-10 electrons' trajectories for `electron.dat`, run the coherent radiation calculation
(`Simulation::run_simulation`) over the whole beam, and export `.dat` files plus `run_log.txt` into a per-run
directory `~/<output_folder>/YYYYMMDD_HHMMSS/`. Python scripts in `py_scripts/` plot those outputs.

**Open issues** (detail in `theory/dev_notes.md`):
- The single-electron Thomson-dipole benchmark (electron at rest, full-4π spherical detector, weak field →
  `1+cos²θ` for circular polarization) has **not** been reproduced yet. Don't trust absolute intensities/angular
  patterns for more complex setups until it is.
- `plot_observables.py`'s `Lambda_zz/P_z` matches the theoretical `m/omega_N` in magnitude (~0.6%) but has the
  **opposite sign** on the backward-detector config. Planned check: flip only `detector_direction_theta` 0↔1 π
  and the sign of `average_pz`, keeping `|l|` and the electron count.
- `LaguerreGaussLaser`'s mode formula isn't physics-reviewed against the literature (Allen et al. PRA 45, 8185
  (1992); Siegman ch. 17), though its derivatives match finite differences. `PhysUtils::non_linear_Thomson_formula`
  is flagged in-code as still needing review.

## Machines

Two laptops share a checkout through Dropbox (under `~/Dropbox`); two servers are synced through git. Python is
`.venv` on the laptops and Anaconda on the servers (`pyrefly.toml`/`pyrightconfig.json` point at `.venv`, so editor
type-checking won't find the interpreter on the servers). Compilers differ between machines.

## Commands

```
cmake -B build/ && cmake --build build/                    # build (out-of-source enforced)
cmake --preset release && cmake --build --preset release  # same via preset; `debug` preset opts out of Release
python3 py_scripts/compile_and_run.py [config_file]        # configure + build + run
./bin/coherent_thomson_solver config/config.cfg            # run (config argument required)
```

**Always built as `Release`** (`-O3 -march=native -mtune=native`, C++20, `-Wall -Wextra -Wpedantic`).
`CMakeLists.txt` forces it over any cached or IDE-supplied `CMAKE_BUILD_TYPE`. Another type needs
`-DCOHERENT_THOMSON_ALLOW_NON_RELEASE=ON` (the `debug` preset sets it), which stays cached until you pass `=OFF`.
Debug is ~12–13x slower, so check the configure step's `Build type: ...` line before trusting any timing.
`-march=native` means a binary must be built on the machine that runs it.

**Build traps, check these first:**
- **Incremental rebuilds have missed header dependencies** (not root-caused). After editing a widely-included
  header (e.g. `phys_utils.hpp`), rebuild with `cmake --build build/ --clean-first`.
- **Every build writes the same `bin/coherent_thomson_solver`,** and a no-op build doesn't relink. To compare
  builds, use `--clean-first` and copy each binary out right after building. The `release` and `debug` presets
  share one `build/` directory and CMake cache.
- **On the laptops, Dropbox syncs `build/` and `bin/`** (gitignored, not Dropbox-ignored), so the other laptop's
  binaries/objects and changed timestamps can arrive. This may explain the missed rebuilds. If a build acts
  strangely, `--clean-first` on this machine before debugging anything else.
- **GCC 15.2's SLP vectorizer miscompiles the long-distance radiation functions** (wrong sums at every frequency
  after the first). Both carry `__attribute__((optimize("no-tree-slp-vectorize")))`. GCC 13.3 doesn't have the
  bug, but the workaround is **kept on purpose on every machine**: don't remove it or make it compiler-conditional.
  Details: `theory/implementation_details.md` §6.

**Iterating on code changes:** copy a config (don't edit `config/config.cfg` in place unless the task is about the
default config) and drop `beam_particle_count` to e.g. `50`; run time scales with electron count. Configs:
- `config/config.cfg`: default, for current experiments; its parameters change freely.
- `config/config_cross_check.cfg`: parameter-matched to the Python reference's `main.py` `INPUTS` (mapping and
  date in its header). Use it for every C++-vs-Python comparison, and update it when those `INPUTS` change.
- `config/config_initial_momentum.cfg`.

Formatting: `clang-format -i` on touched files (Google style, 120 cols, 2-space indent). **There is no test suite**
and no enforced lint command.

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
python3 py_scripts/radiation/plot_field.py <long|short|boundary|total> <mu> <nu> [--incident] [run_folder]
python3 py_scripts/radiation/plot_all_components.py [--incident]   # plot_field over all 4 ranges x 6 (mu,nu)
python3 py_scripts/radiation/plot_point_spectrum.py <long|short|boundary> <mu> <nu>   # dense_frequency_spectrum=true
python3 py_scripts/radiation/plot_spherical_components.py <long|short|boundary> <E|B> <r|theta|phi>  # spherical only
python3 py_scripts/radiation/plot_observables.py [--incident] [run_folder]   # rectangular/circular facing ±Oz only
```
`--incident` analyzes `incident_field.dat` (the incident laser's analytic field, written only for
rectangular/circular detectors) instead of `radiation_field.dat`. For `mu > nu` the scripts use
`F^{nu mu} = -F^{mu nu}`; `mu == nu` is rejected. Field extraction: `E_i = c·F^{i0}`, `B = (F^{32}, F^{13}, F^{21})`.

## Architecture

### Layout

- `src/app/main.cpp` only orchestrates. All logic is in the static library `coherent_thomson_core` (`src/core/`),
  with public headers in `src/core/include/<module>/`. `math_utils/`, `phys_utils/`, `io_utils/` are header-only.
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

- The laser **always propagates along Oz**; there is a single frame everywhere. (Arbitrary laser direction was
  removed in `d04bb89`; the old code is on the `general-laser-direction-legacy` branch.)
- The detector has its own direction (`detector_direction_theta/phi`, default `(0,0)` = forward).
- The beam cylinder's axis follows the mean momentum (`average_px/py/pz`), not Oz. `beam_center_x/y/z` are applied
  **before** that rotation, so `beam_center_z` moves along the direction of motion.
- `CircularDetector`'s radial grid is **equal-area**: `r_i = sqrt(R_min² + i*(R_max²-R_min²)/(N_R-1))`. Anything
  reconstructing `r` from `i` must use this.
- Every detector collapses to one exact point when its grid counts are `1`.

### Radiation calculation (the hot path)

All six `radiation_formula` modes, their formulas, the optimizations, precision limits and timings are in
`theory/implementation_details.md`. **Keep it in sync when changing `radiation.cpp`.** Current formulas:
`theory/FT_Faraday_tensor-direct_and_simplified_forms-v2.md` (supersedes the unsuffixed v1).

- **Three terms**, exported as `LR_`/`SR_`/`BR_F<mu><nu>` columns: long-range `F_l`, short-range `F_s`, boundary
  `F_b`. **Only the sum is physical**; e.g. on-axis `F_l^{03}` and `F_b^{03}` are individually large and cancel.
- `radiation_formula` is **required** (no default; unknown names throw): `simplified`, `direct` (needs the stored
  4-acceleration; `F_b ≡ 0`), `long_distance_{simplified,direct}` (far screens only; no `F_s`; error ~`1/D`), and
  `long_distance_{simplified,direct}_approx` (stepped phase, fastest, matches unstepped LD to ≤1e-7). In the LD
  modes the `LR_`/`BR_` split differs from the exact one, so compare only totals.
- **Cross-checking `simplified` against `direct` is the main correctness tool.** They agree to the trapezoid's
  `O(dτ²)` (4x per doubling of `trajectory_NT`) down to a 1 λ screen. A mismatch of order `1/(kR)` means a
  short-range bug. Simplified `F03` near the axis and weak harmonics need a finer `trajectory_NT`.
- `general_factor = q_0 / (2π · 4πε₀c²)` uses the **signed** charge `q_0`. The direct `F_s` carries an explicit
  `c²` in `short_range_prefactor_direct` to compensate.
- **Quadrature:** trapezoidal weights apply to `F_l`/`F_s` only; `F_b` is an exact endpoint evaluation with **no**
  weight.
- **Keep these, they're load-bearing:** prefactors stay isolated in `radiation.cpp`'s anonymous namespace, and
  `radiation_phase_argument` must stay **frequency-independent** (the phase uses a frequency recurrence, valid
  because `frequencies_list` is evenly spaced). Loop order is screen point outer, tau inner. In the LD modes the
  constant phase `e^{ik(|x_0| + r^0(tau_m))}` is applied at reconstruction, **not** inside the tau loop
  (precision, §4.3). The stepped loops must stay vectorizable (check `-fopt-info-vec`, §6).
- The whole beam is generated upfront and held in memory (deliberate, for now).

### Frequencies

- Default: `frequencies_list[i] = non_linear_Thomson_formula(k1, q, n2, N_harmonics_min + i)`, with `n2` the
  detector direction and `q` the **ponderomotively dressed** momentum with `<a²> = a0²/2` for every polarization.
  **Do not change it to `a0²`**; that was a real bug, confirmed three ways (dev notes).
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
- **A key with a blank value is dropped**, and `config.at(...)` then throws `std::out_of_range`.
- `w0` is valid only with `laser_type=laguerre_gauss`.
- Unknown keys are ignored, so stale keys from old configs don't error.

`run_log.txt` re-reads config keys the way the factories do (`read_scaled`) rather than adding getters to the
physics classes.

## Physics notes worth knowing before comparing results

- **Rectangular vs. spherical screens differ for real** (exact `R`, no far-field approximation). They agree only
  when the Fresnel number `N_F = a²/(Dλ) << 1`; growing `D` at a fixed *angular* window makes `N_F` larger. Compare
  `Re/Im(F01)`, not `|F01|`.
- Per contribution, `B` is exactly ⊥ `n0` but `E` is not (`|B_r|/|B_θ|` ~ 1e-5, `|E_r|/|E_θ|` ~ 0.1): a useful
  regression check for `plot_spherical_components.py`.
- The independent Python reference is `~/Dropbox/work/bin/python/Superradiant_Thomson` on the laptops (not on the
  servers; ask the user). Cross-checking against it found several past bugs (missing trapezoid weight, `e_0` vs
  `q_0` sign, envelope `σ` convention, `<a²>`).

Full investigation history, measurements, and rationale: `theory/dev_notes.md`.
