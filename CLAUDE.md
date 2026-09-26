# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

Pipeline (`src/app/main.cpp`): build a config-driven laser pulse and detector screen, generate an electron beam,
integrate the first up-to-10 electrons' trajectories for `electron.dat`, run the coherent radiation calculation
(`Simulation::run_simulation`) over the whole beam, and export `.dat` files plus `run_log.txt` into a per-run
directory `~/<output_folder>/YYYYMMDD_HHMMSS/`. Python scripts in `py_scripts/` plot those outputs.

**Open issues** (see "Physics notes" below and `theory/dev_notes.md` for detail):
- The single-electron Thomson-dipole benchmark (electron at rest, full-4π spherical detector, weak field →
  `1+cos²θ` for circular polarization) has **not** been reproduced yet. Don't trust absolute intensities/angular
  patterns for more complex setups until it is.
- `plot_observables.py`'s `Lambda_zz/P_z` matches the theoretical `m/omega_N` in magnitude (~0.6%) but has the
  **opposite sign** on the backward-detector config. Planned check: a controlled forward-vs-backward comparison
  (flip only `detector_direction_theta` 0↔1 π and the sign of `average_pz`, same `|l|` and electron count).
- `LaguerreGaussLaser`'s mode formula isn't physics-reviewed against the literature (Allen et al. PRA 45, 8185
  (1992); Siegman ch. 17), though its derivatives match finite differences. `PhysUtils::non_linear_Thomson_formula`
  is flagged in-code as still needing review.

## Commands

```
cmake -B build/ && cmake --build build/                    # build (out-of-source enforced)
cmake --preset release && cmake --build --preset release  # same, but forces CMAKE_BUILD_TYPE=Release
python3 py_scripts/compile_and_run.py [config_file]        # configure + build + run
./bin/coherent_thomson_solver config/config.cfg            # run (config argument required)
```

Build types: `Release` (default, `-O3 -march=native -mtune=native`), `Debug` (`-g`, no `-O`), `RelWithDebInfo`,
`MinSizeRel`. C++20; `-Wall -Wextra -Wpedantic` except on `RelWithDebInfo`/`MinSizeRel`.

**Build traps, check these first:**
- **The build type sticks in the cache.** `cmake -B build/` and `compile_and_run.py` (which never passes
  `-DCMAKE_BUILD_TYPE`) keep whatever `build/CMakeCache.txt` already holds, and an IDE can silently set it to
  `Debug`. A Debug build is ~12–13x slower. Before trusting any timing, run `grep CMAKE_BUILD_TYPE
  build/CMakeCache.txt`, or use the `release` preset.
- **Incremental rebuilds have missed header dependencies.** After editing a widely-`#include`d header (e.g.
  `phys_utils.hpp`), rebuild with `cmake --build build/ --clean-first`. Not root-caused.

**Iterating on code changes:** copy a config (don't edit `config/config.cfg` in place unless the task is about the
default config) and drop `beam_particle_count` to e.g. `50` (default `16384`); run time scales with electron
count. Configs: `config/config.cfg` (default; parameter-matched to the Python reference, so use it for
cross-checks), `config/debug.cfg` (debug mode, see below), `config/config_initial_momentum.cfg`.

Formatting: `clang-format -i` on touched files (`.clang-format`: Google style, 120 cols, 2-space indent).
`py_scripts/` is type-checked against `.venv` via `pyrefly.toml`/`pyrightconfig.json` (editor config only, no
enforced lint command). **There is no test suite.**

### Plotting

With no arguments, every script plots the latest run (`utils/run_output_utils.find_latest_output_file`). The
optional trailing argument is a **run folder**, not a `.dat` path: each script appends its own filename and fails
clearly if it's missing. PNGs go to `<run>/png_folder/<module>/`. Each script sets `MODULE_NAME` from its own
folder name; don't hardcode it.

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
python3 py_scripts/radiation/plot_observables.py [--incident] [run_folder]   # rectangular/circular only
python3 py_scripts/debug/plot_integrand.py <long|short|boundary> <mu> <nu>   # debug=true
python3 py_scripts/debug/plot_exponent.py                                    # debug=true
```
`--incident` analyzes `incident_field.dat` (the incident laser's own analytic field, from
`Radiation::export_incident_field_fourier`, written only for rectangular/circular detectors) instead of
`radiation_field.dat`. For `mu > nu` the scripts use `F^{nu mu} = -F^{mu nu}`; `mu == nu` is rejected.

## Architecture

### Layout

- `src/app/main.cpp`: the executable, which only orchestrates. All logic is in the static library
  `coherent_thomson_core` (`src/core/`). Public headers live in `src/core/include/<module>/`, with `.cpp` files in
  `src/core/<module>/`.
- `math_utils/`, `phys_utils/`, `io_utils/` are header-only.
- **Class + factory pattern** (`laser/`, `detector/`): `<thing>.hpp/.cpp` holds the class; `<thing>_factory.hpp`
  declares `create_<thing>(const ConfigMap&, ...)`, defined in `<thing>_factory.cpp`, where *all* config-key
  parsing and unit handling lives. `particle/`'s factory is `generate_cylinder_beam(config, tau_0, d_tau, N_tau)`,
  which returns `std::vector<Electron>`. `simulation/` is plain free functions (`init_simulation_parameters`,
  `run_simulation`).
- Each module with output has a `<module>_plotter.hpp/.cpp` exporter, called from `main.cpp`.
- **Adding a physics component:** class + factory as above, add the `.cpp` to `add_library(...)` in
  `src/core/CMakeLists.txt`, and wire the factory into `main.cpp`.
- `py_scripts/<module>/` mirrors the C++ modules, with one plotting script per exported `.dat`. Scripts run as
  plain files (not `python -m`): each inserts `py_scripts/` onto `sys.path` and imports `utils.*`, an implicit
  namespace package with no `__init__.py`. Scripts in the same folder import each other by filename. There is no
  shared constants module between C++ and Python: `C_LIGHT = 137.035999084`, `EPSILON_0 = 1/(4π)` etc. are
  duplicated in the scripts, so keep them in sync with `phys_utils.hpp`.

### Namespaces (`Core::<CamelCase of folder>`)

| Module | Key contents |
|---|---|
| `MathUtils` | `FourVector<T>`/`FourTensor<T>` (`Real*`/`Complex*` aliases; `Complex = std::complex<double>`), Minkowski contractions, `unpack_bivector`, 3D rotation helpers |
| `PhysUtils` | `AtomicUnits::` constants (`c = 137.035999084`, `q_0 = -1` signed charge, `e_0 = +1` magnitude, `epsilon_0 = 1/(4π)`, SI conversions); `dressed_momentum`, `non_linear_Thomson_formula` |
| `IoUtils` (`ConfigMap` is in `Core`) | config parsing, `convert_unit_to_number`, per-key accessors, `make_run_output_directory`, `copy_config_to_run_directory` |
| `Particle` | `Electron` (RK4 Lorentz-force integrator; every `State` records position, momentum, and 4-acceleration), beam generation, trajectory export |
| `Laser` | `LaserField` base + `PlaneWaveLaser`/`LaguerreGaussLaser`; `create_laser` dispatches on `laser_type` |
| `Detector` | `Detector_2D` base + `Rectangular`/`Circular`/`SphericalDetector`; `create_detector` |
| `Simulation` | `init_simulation_parameters` (frequency list, fundamental, dressed `q`, trajectory grid), `run_simulation`, `PackedFaraday`/`Faraday` field containers |
| `Radiation` | `compute_radiation` (per-electron hot path), `plot_radiation_field`, `export_incident_field_fourier` |
| `Debug` | per-tau integrand diagnostics (see "Debug mode") |
| `Logging` | `write_run_log` → `run_log.txt` |

### Frames and geometry

- The laser **always propagates along Oz**; `unity_n`/`epsilon_1`/`epsilon_2` are constants. (Arbitrary laser
  direction was removed in `d04bb89`; the old implementation is on the `general-laser-direction-legacy` branch.)
  There is a single frame everywhere, with no lab/canonical distinction.
- The detector has its own direction (`detector_direction_theta/phi`, default `(0,0)` = along the laser, i.e.
  forward). `Detector_2D::local_rotation` orients the screen orthogonal to it.
- The beam cylinder's axis follows the beam's mean momentum (`average_px/py/pz`), not Oz. `beam_center_x/y/z`
  offsets are applied **before** that rotation, so `beam_center_z` translates along the direction of motion.
- `CircularDetector` spaces its radial grid by **equal area**: `r_i = sqrt(R_min² + i*(R_max²-R_min²)/(N_R-1))`.
  Anything reconstructing `r` from `i` must use this formula.
- Every detector type collapses to one exact point when its grid counts are `1`.
  `Detector_2D::restrict_to_first_point` exists for the dense-spectrum case below.

### Radiation calculation (the hot path)

- `compute_radiation` sums, per `(screen point, tau, frequency)`, prefactor × bivector
  `n0^α u^β − n0^β u^α` × `exp(i·freq·(x⁰+R))`. Values accumulate in the packed 6-component form
  (`PackedFaraday`) and are unpacked to a 4x4 tensor after the per-thread accumulators are summed.
  `run_simulation` splits electrons across `num_threads` (`0` = all hardware threads) and scales the result by
  `general_factor = q_0 / (2π · 4πε₀c²)`, where `q_0` is the signed charge, not `e_0`.
- **Three terms**, exported separately as `LR_`/`SR_`/`BR_F<mu><nu>` columns in `radiation_field.dat`: long-range
  `F_l`, short-range `F_s`, and boundary `F_b`. **Only the sum `F_l+F_s+F_b` is physical.** Individual terms are
  artifacts of the integration-by-parts derivation; e.g. on-axis `F_l^{03}` and `F_b^{03}` are individually large
  and cancel. The plotting scripts' `total` range and `plot_observables.py` always use the sum.
- **Current formulas: `theory/FT_Faraday_tensor-direct_and_simplified_forms-v2.md`**, which supersedes the
  unsuffixed v1 doc for the simplified `F_s`. Simplified form: `F_l ∝ (-ik/R)(nᵅuᵝ−nᵝuᵅ)`,
  `F_s ∝ (sᵅuᵝ−sᵝuᵅ)/R²` with `s = (0, 𝐧)` (its own bivector, `short_range_bivector_element`; not the `(n, u)`
  one), `F_b ∝ [(nᵅuᵝ−nᵝuᵅ)/(R(n·u))]` evaluated between the two endpoints. The direct `F_s` carries an explicit
  `c²` in `short_range_prefactor_direct`, because the doc's direct `F_s` has a bare `q` while `general_factor`
  applies `q/c²` to everything.
- `radiation_formula` selects `simplified` (default: velocity-only integrand plus the boundary term `F_b`) or
  `direct` (Liénard–Wiechert FT that needs the stored 4-acceleration; `F_b ≡ 0`). Cross-checking the two
  formulas against each other is the main correctness tool. With the v2 `F_s` they agree to the trapezoid's
  `O(dτ²)` error (4x smaller per doubling of `trajectory_NT`) down to a 1 λ screen distance. A mismatch of order
  `1/(kR)` means a short-range bug. Simplified `F03` near the axis is a ~1000x cancellation of `F_l` and `F_b`,
  so it needs a finer `trajectory_NT` than the other components to reach the same accuracy.
- **Quadrature:** trapezoidal weights (`d_tau`, halved at the endpoints) apply to `F_l`/`F_s` only. `F_b` is an
  exact endpoint evaluation (nonzero only at the first/last tau, opposite signs) and gets **no** weight.
- The prefactors are isolated in `radiation.cpp`'s anonymous namespace (`long_range_prefactor`,
  `short_range_prefactor`, `*_direct`, `boundary_prefactor`, `radiation_phase_argument`) so they can be swapped
  independently. `radiation_phase_argument` must stay **frequency-independent**, because of the recurrence below.
- **Performance tuning (keep it):** the loops run screen point outer, tau inner, so each screen point accumulates
  in local tensors and `field` is written once per `(i_d, i_freq)`. Because `frequencies_list` is always evenly
  spaced, the phase factor is computed by recurrence (`cexp *= cexp_step`) instead of a `std::polar` call per
  frequency. This is guarded by a per-electron evenly-spaced check. Loops are deliberately left uninstrumented.
  On a matched config the Release build is ~4x faster than the Python reference. An older "C++ is 3x slower"
  result came from a Debug build.
- `run_simulation` prints per-thread timings and aggregate CPU-seconds per electron (per screen point), computed
  from the **sum of per-thread times**, not wall-clock time.
- The whole beam is generated upfront and held in memory (a deliberate simplification for now).

### Frequencies

- Default: `frequencies_list[i] = non_linear_Thomson_formula(k1, q, n2, N_harmonics_min + i)` for `N_harmonics`
  harmonics. `k1` is along Oz, `n2` is the detector direction (not the electron's), and `q` is the
  **ponderomotively dressed** momentum `PhysUtils::dressed_momentum`, with cycle-averaged `<a²> = a0²/2` for
  every polarization. **Do not change it to `a0²`.** That was a real bug, and the `/2` has been confirmed three
  independent ways (see dev notes).
- `dense_frequency_spectrum=true`: a linear scan of `N_omega` points from `omega_min` to `omega_max`, meant for a
  single screen point (resolving a line shape). The bounds accept only `omega_laser` (raw multiple of the incident
  frequency) or `first_harmonic_frequency` (multiple of the emitted, Doppler-shifted fundamental). If the detector
  has more than one point, `main.cpp` warns and restricts it to point 0.
- Internally frequencies are `omega/c` (`k`). `radiation_field.dat`'s `omega` column is **`omega/omega_1`**
  (normalized to `fundamental_frequency`, the N=1 harmonic of `q`), not raw omega. `run_log.txt` records
  `fundamental_frequency` in raw atomic units, and `plot_observables.py` reads it from there to recover true
  `omega` for its `1/omega` factors.

### Laser

- `get_faraday_tensor` is a single non-virtual implementation built on the pure-virtual `complex_amplitude`,
  which returns `{amplitude, d/dx, d/dy}` already scaled by `E0_c` (`{0, 0}` derivatives for `PlaneWaveLaser`).
  `Ez`/`Bz` are built from those derivatives: this is first-order paraxial, with a `div E` residual of ~1e-6 for
  the fundamental Gaussian and up to ~1e-2 for higher `p`/`|l|`. That is expected. The carrier sign `exp(-iφ)` and
  the `+i/k` coefficient must flip together.
- Polarization `zeta_1`/`zeta_2` are complex (`laser_zeta_{1,2}_{re,im}` keys) and **normalized in
  `create_laser`**, so `a0` alone sets the field strength. Linear: `(1,0),(0,0)`. Circular: `(1,0),(0,1)`.
- LG: the radial profile uses `hypergeometric_1F1_neg_int_a` (`generalized_laguerre` exists but is unused). `Npn`
  is a deliberately custom normalization. The azimuthal factor is the polynomial `(x + i·sign(l)·y)^|l|` to stay
  smooth on-axis. The beam waist is at the origin.
- Envelope: Gaussian wings `exp(-Δφ²/wing_sigma²)` (not `/(2σ²)`, which was a fixed bug) around a flat top.
- **Known bugs/hardcodes:** `laser_delay` has no effect (the constructor overwrites it with
  `wing_sigma_cutoff * wing_sigma`). `tau_0_traj` is hardcoded to `0.0` in `init_simulation_parameters`. The
  heatmap snapshot time is keyed off `wing_sigma_cutoff*wing_sigma` for that reason. Revisit all three together.
- `export_field_heatmap_z0`'s window is `±beam_cylinder_radius` and ignores `beam_center_x/y`. Only when the
  radius is `0` does it fall back to `field_heatmap_*` (plane wave) or `±2 w0` (LG).

## Config file format

Flat `key value [unit]` lines. `IoUtils::convert_unit_to_number` resolves units (`lambda`, `pi`, `mc`,
`cycles_adim`, `omega_laser`, `w0`, `a.u.`, ...) against the laser's own wavelength and frequency. Check that
function before adding a unit.
- **An unrecognized unit only warns** on stderr and assumes `1.0`, so a typo fails silently.
- **A key with a blank value is dropped**, so `config.at(...)` then throws `std::out_of_range`. Always give every
  key a real value.
- `w0` is valid only with `laser_type=laguerre_gauss` (it resolves through `laser_lg_w0`'s own unit).
- Unknown keys are ignored, so stale keys from old configs (e.g. `laser_nx`, `print_field_in_canonical_frame`)
  don't error.

## Debug mode

`debug=true` (use `config/debug.cfg`) requires `beam_particle_count=1`, a 1x1 rectangular detector, and
`radiation_formula=simplified`; `main.cpp` throws otherwise. It evaluates at the fundamental only and writes
`debug_integrand.dat` (per-tau `LR`/`SR`/`BR` bivector terms, trapezoid weights applied as in production) and
`debug_exponent.dat` (the shared phase factor). The normal pipeline still runs, so the tau-sum can be
cross-checked against `radiation_field.dat` once `general_factor` is divided out. Compare at the exact
fundamental (`dense_frequency_spectrum=false`, `N_harmonics=1`, row 0): the sums then agree to print precision
(~1e-7). A dense-scan row that is only near `ω/ω₁ = 1` differs by O(1) because of the long-trajectory phase.

**`debug/debug_radiation.cpp` deliberately duplicates `compute_radiation`'s per-tau math instead of sharing it.
Any change to the prefactors, weights, or boundary term must be mirrored there**, or the cross-check silently goes
stale.

`debug/plot_exponent.py` and `particle/plot_trajectory.py` plot against `tau/T`, reading the laser period from
that run's own `config.cfg` snapshot via `get_laser_period`. `debug/plot_integrand.py` uses raw `tau`.

## Outputs and plotting conventions

- `run_log.txt` summarizes every parameter and derived quantity (atomic units, `lambda`/`w0`, and SI): the
  fundamental and its shift ratio, `q`, and the frequency list. It re-reads config keys the way the factories do
  (`read_scaled`) rather than adding getters to the physics classes. `plot_observables.py` appends its
  screen-integrated tables to it, replacing only its own section by exact header match (`_find_run_log_sections`).
  Don't reintroduce a "truncate after first marker" approach.
- `electron.dat` holds up to 10 electrons tagged by `electron_id`, with trailing `a0..a3` acceleration columns.
  `plot_trajectory.py` draws one panel per four-vector component, since sharing an axis hid real variation, with
  position in `lambda` and momentum in `mc`.
- `radiation/plot_field.py` (and `plot_spherical_components.py`) draw 2x2 Re/Im/Abs/Phase panels. Re and Im share
  the `±|F|max` scale. Rectangular, circular, and spherical detectors all render as `pcolormesh` on the native grid.
  On circular detectors Re/Im/Abs use `gouraud` shading on cell centers, while **Phase stays flat-shaded** because
  interpolating across the ±π wrap is misleading. Colorbar alignment requires freezing the layout (`fig.canvas.draw()`,
  `set_layout_engine(None)`, `cbar.ax.set_axes_locator(None)`) before `set_position`. `constrained_layout` alone
  misaligns them.
- **Spherical detectors near the full 4π:** stereographic projection is singular at θ=π (`nan` in
  `detector_stereographic.dat` is expected). `get_spherical_plot_grid` falls back to a `(θ, φ)` map.
- Heatmaps of LG runs add secondary `x/w0`, `y/w0` axes (`utils/w0_axes_utils.py`, a leaf module that
  deliberately avoids importing `plot_field` to prevent circular imports).
- `plot_observables.py` computes the SAM/OAM/TAM densities `S_z`/`L_z`/`J_z`, the fluxes
  `Sigma_zz`/`Lambda_zz`/`Flux_tot`, and energy `u`/`P_z`, per `theory/numerical_calculation_of_angular_momentum.md`.
  It requires a rectangular or circular screen **facing ±Oz** (it raises otherwise). The `L_z` operator is
  `x∂y−y∂x` (rectangular, `np.gradient`) or a periodic `∂φ` (circular). The angular-momentum quantities divide by
  the true `omega` (not `omega/omega_1`). Sanity checks it prints: each flux/density ratio ≈ `±c` (holds to ~3e-5),
  and `Lambda_zz/P_z` ≈ `m/omega_N` (magnitude OK, sign open, see Open issues).
- `get_faraday_tensor`-style extraction: `E_i = c·F^{i0}`, `B = (F^{32}, F^{13}, F^{21})`.

## Physics notes worth knowing before comparing results

- **Rectangular vs. spherical screens differ for real** (exact `R`, no far-field approximation). They agree only
  when the Fresnel number `N_F = a²/(Dλ) << 1`. Growing `D` at a fixed *angular* window makes `N_F` larger, not
  smaller. Compare `Re/Im(F01)`, not `|F01|`.
- Per contribution, `B` is exactly ⊥ `n0` but `E` is not: `|B_r|/|B_θ|` ~ 1e-5 while `|E_r|/|E_θ|` ~ 0.1. That is
  a useful regression check for `plot_spherical_components.py`.
- The independent Python reference is `~/Dropbox/work/bin/python/Superradiant_Thomson`, and `config/config.cfg`
  is matched to its `main.py` `INPUTS`. Several past bugs were found by cross-checking against it (missing
  trapezoid weight, `e_0` vs `q_0` sign, envelope `σ` convention, `<a²>`); see `theory/dev_notes.md`.
- `MathUtils::mirror_antisymmetric_in_place` is dead code.

Full investigation history, measurements, and rationale for everything above: `theory/dev_notes.md`.
