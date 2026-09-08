# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

End-to-end today (`src/app/main.cpp`): build a config-driven laser pulse and detector screen, generate an electron
beam, integrate the first up to 10 electrons' trajectories for plotting, then run the coherent radiation spectrum
calculation (`Simulation::run_simulation`) over the whole beam and export `.dat` files for plotting. The
per-electron physics `run_simulation` delegates to, `Radiation::compute_radiation`, accumulates the full
antisymmetric Faraday bivector tensor (long-range and short-range amplitude terms, normalized) per
frequency/screen point — see "Known gaps" below for what's still missing.

## Commands

Build (out-of-source build enforced by CMake):
```
cmake -B build/ && cmake --build build/
```

Configure, build, and run in one step (binary lands at `bin/coherent_thomson_solver`):
```
python3 py_scripts/compile_and_run.py [config_file]
```
Running the binary directly requires a config file argument:
```
./bin/coherent_thomson_solver config/coherent_thomson.cfg
```

When testing a code change (not validating physics), run with a drastically reduced
`beam_particle_count` (e.g. `50` instead of the repo default `2000`) — simulation time scales with
the electron count and a full run can take tens of seconds to minutes, dominating iteration time for
no benefit while just checking that something builds/runs/exports the right files. Copy the config
rather than editing `config/coherent_thomson.cfg` in place, unless the task is specifically about
changing the default config.

Visualize `.dat` outputs (no arguments — each locates the most recent `<output_folder>/YYYYMMDD_HHMMSS` run via
`py_scripts/run_output_utils.py`, so plotting always targets the last solver run):
```
python3 py_scripts/plot_laser_field.py
python3 py_scripts/plot_detector_stereographic.py
python3 py_scripts/plot_electron_trajectory.py
python3 py_scripts/plot_radiation_field.py
python3 py_scripts/plot_detector_scatter.py       # only has output if plot_detector_scatter=true in the config
python3 py_scripts/plot_electron_beam_scatter.py  # only has output if plot_beam_scatter=true in the config
python3 py_scripts/plot_field_heatmap_z0.py       # only has output if plot_field_heatmap=true in the config
python3 py_scripts/plot_point_spectrum.py <long|short> <mu> <nu>  # meaningful output only if dense_frequency_spectrum=true
python3 py_scripts/plot_debug_integrand.py <long|short> <mu> <nu> # meaningful output only if debug=true
python3 py_scripts/plot_debug_exponent.py                        # meaningful output only if debug=true
python3 py_scripts/plot_angular_momentum_flux.py                 # rectangular/circular detectors only, see its module docstring
python3 py_scripts/plot_spherical_field_components.py <long|short> <E|B> <r|theta|phi>  # spherical detectors only
```

Build types (`-DCMAKE_BUILD_TYPE=...`, default `Release`): `Release` (`-O3 -march=native -mtune=native`), `Debug`
(`-g`), `RelWithDebInfo`, `MinSizeRel`. C++20. `-Wall -Wextra -Wpedantic` on all types except
`RelWithDebInfo`/`MinSizeRel`.

Run `clang-format -i` on touched files before committing (`.clang-format` at repo root: Google style, 120 cols,
2-space indent). `py_scripts/` is type-checked against a `.venv` interpreter via `pyrefly.toml`/`pyrightconfig.json`
(no enforced lint command yet, just editor/CLI type-checker config).

There is no test suite yet.

Sync the repo to a remote build/run server (excludes `build/`, `src/build/`, `bin/`, `.venv/`, `.cache/`,
`.claude/`, `.vscode/`, and the local `.code-workspace` file — `py_scripts/` itself *is* synced):
```
./sync_to_remote.sh user@remote:/path/to/CoherentThomson/
```

## Architecture

### Directory / module strategy

- `src/app/` — the executable. Owns `main.cpp` only; all logic lives in the core library it links against.
- `src/core/` — built as static library `coherent_thomson_core`. Public headers under `src/core/include/`, mirrored
  by a `.cpp` under `src/core/<module>/` for types needing out-of-line definitions.
  - `math_utils/`, `phys_utils/`, `io_utils/` are **header-only** (no `.cpp`).
  - `laser/`, `detector/` follow a **class + factory** split: `<thing>.hpp`/`.cpp` declares/defines the class,
    `<thing>_factory.hpp` only *declares* `create_<thing>` (defined in `<thing>_factory.cpp`, where all config-key
    parsing and unit handling happens). `particle/` mostly follows this, but its factory
    (`generate_cylinder_beam`) takes a precomputed rotation matrix + `tau_0`/`d_tau`/`N_tau` as scalars instead of
    a built object, and returns a `std::vector<Electron>`.
  - Each module with plottable output has a matching `<module>_plotter.hpp`/`.cpp` exporter, called from
    `main.cpp` after the object is built, writing into the same per-run directory from
    `IoUtils::make_run_output_directory`.
  - `simulation/` is a single free-function module (no class+factory pattern): `init_simulation_parameters(config,
    laser)` and `run_simulation(config, laser, detector, electron_beam, frequencies_list, num_threads)`.
  - `radiation/` holds `compute_radiation(electron, laser, frequencies_list, detector, field)`, called once per
    electron from inside `run_simulation`'s per-thread loop, plus a `radiation_plotter.hpp`/`.cpp` exporter
    (`plot_radiation_field`) following the same plottable-output pattern as the other modules.
- Adding a new physics component: follow the class + `.cpp` + `create_<thing>(const ConfigMap&)` factory pattern
  above, add the new `.cpp` to `add_library(...)` in `src/core/CMakeLists.txt`, and wire the factory call into
  `main.cpp`.

### Namespace structure

All core-library code lives under `Core`; subdirectories of `src/core/` map to sub-namespaces
(lower_snake_case folder → CamelCase namespace):

| Directory | Namespace | Contents |
|---|---|---|
| `math_utils/` | `Core::MathUtils` | `FourVector<T>`/`FourTensor<T>` (+ `Real*`/`Complex*` aliases), Minkowski contractions, 3D rotation helpers, constants |
| `phys_utils/` | `Core::PhysUtils::AtomicUnits` | Physical constants in atomic units |
| `io_utils/` | `Core::IoUtils` (`ConfigMap` alias lives in `Core`) | Config parsing, unit conversion, per-key laser/beam accessors, `CylinderBeamParams`, `make_run_output_directory` |
| `particle/` | `Core::Particle` | `Electron` (RK4 Lorentz-force integrator, always records its trajectory), `generate_cylinder_beam`/`generate_electron`, `plot_particle_trajectory` |
| `laser/` | `Core::Laser` | `LaserField` base (Gaussian-flat-top temporal envelope, direction/polarization — `zeta_1`/`zeta_2` are complex, e.g. `(1,0)`/`(0,1)` for circular — caches a 4x4 `rotation_matrix`, from which `epsilon_1`/`epsilon_2`/`unity_n` are derived as its columns; `get_faraday_tensor` is a single non-virtual implementation shared by every derived type, built on the pure-virtual `complex_amplitude` customization point) + `PlaneWaveLaser`/`LaguerreGaussLaser` derived types (each implementing `complex_amplitude`, already scaled by `E0_c`, returning `std::tuple<Complex, Complex, Complex>` — `{amplitude, d/dx_loc, d/dy_loc}`; `PlaneWaveLaser` has no transverse profile so its derivatives are always `{0, 0}` — see "Known gaps"), `create_laser` (returns `std::unique_ptr<LaserField>`, dispatches on `laser_type`), `export_field_vs_phase`, `export_field_heatmap_z0` (canonical-frame z=0 transverse snapshot — see "Known gaps") |
| `detector/` | `Core::Detector` | `Detector_2D` base + `RectangularDetector`/`SphericalDetector`/`CircularDetector` (each built orthogonal to its own canonical-frame direction, then rotated together with the laser via its shared 4x4 `rotation_matrix`), `create_detector` (takes the `LaserField`), `plot_detector` |
| `simulation/` | `Core::Simulation` | `init_simulation_parameters`, `Faraday`/`RadiationField` (full 4x4 tensor, post-reduction) + `PackedFaraday`/`PackedRadiationField` (6-element packed bivector, accumulation-time) + `run_simulation` (multithreaded, partitions beam across `num_threads`) |
| `radiation/` | `Core::Radiation` | `compute_radiation` — one electron's contribution to `Simulation::PackedRadiationField` (packed antisymmetric Faraday bivector, long/short-range amplitudes, summed over trajectory points/screen points/frequencies); `plot_radiation_field` exporter |
| `debug/` | `Core::Debug` | `export_radiation_integrand`/`export_radiation_phase` — diagnostic-only, deliberately duplicated (not shared) reimplementation of `compute_radiation`'s per-tau math; see "Debug mode" below |

### Config file format

`config/coherent_thomson.cfg` is a flat `key value [unit]` text format covering detector geometry, laser
frequency/envelope/direction/polarization, beam particle count/geometry, initial momentum distribution, radiation
spectrum range, and a trajectory-print-frequency knob. Numeric values may carry a unit suffix (`lambda`, `pi`,
`mc`, `cycles_adim`, `omega_laser`, `a.u.`, ...) resolved by `IoUtils::convert_unit_to_number` against the laser's
own wavelength/frequency — check that function before adding a new unit keyword; an unrecognized suffix falls
through to a `std::cerr` warning and an assumed value of `1.0` rather than throwing, so a typo'd unit fails silently
(only visible as a startup warning) instead of erroring out. A key with a blank value (nothing but a trailing
comment) is dropped entirely by the parser rather than stored empty, so `config.at(...)` throws
`std::out_of_range` for it — always give every key a real value, even a placeholder.

One exception to the single-number-plus-unit convention: the laser's polarization coefficients are complex
(`Core::MathUtils::Complex`), so each is split into two plain-number keys instead of one `key value [unit]` line —
`laser_zeta_1_re`/`laser_zeta_1_im` and `laser_zeta_2_re`/`laser_zeta_2_im` — read via
`IoUtils::get_complex_config_value`/`get_laser_zeta` (`io_utils.hpp`), not `convert_unit_to_number`.

### Debug mode (per-tau radiation integrand diagnostics)

`config/coherent_thomson_debug.cfg` is a ready-made config satisfying the constraints below (`debug=true`,
`beam_particle_count=1`, single-point rectangular detector) — use it (or a copy) instead of hand-editing the
default config when exercising this path.

`debug/` (`Core::Debug`, `debug_radiation.hpp`/`.cpp`) is a diagnostic-only module, deliberately kept isolated from
`radiation.cpp`'s production accumulation path: `export_radiation_integrand`/`export_radiation_phase`
**reimplement (duplicate, not share) `compute_radiation`'s per-tau `n0`/`R`/`amp_long`/`amp_short`/bivector-term
math**, so this path can never be affected by, or accidentally affect, the real simulation. Enabled via the
`debug` config key (default `false`); when `true`, `main.cpp` requires (and throws at startup otherwise)
`beam_particle_count=1`, `detector_type=rectangular`, and `rectangular_detector_Nx=rectangular_detector_Ny=1`, and
forces the diagnostic to evaluate at the fundamental only (`sim_par.fundamental_frequency`, already in `k =
omega/c` units, not raw `omega` — see the `radiation_field.dat` `omega` bullet below), independent of
`dense_frequency_spectrum`/`N_harmonics`/`N_harmonics_min`/`N_omega`. Purely additive: the normal pipeline (including
`run_simulation`/`radiation_field.dat`) still runs, so the integrand's own tau-sum can be cross-checked against
the coherently-summed field (confirmed to agree to the file's printed precision, both long-range and short-range,
once `general_factor` and any canonical-frame rotation are accounted for).

Two output files, both one row per trajectory point (`tau`), for the single electron/screen point:
- `debug_integrand.dat` (`export_radiation_integrand`): the 6 independent upper-triangle Faraday bivector
  components' long-range/short-range contribution at that tau — the raw terms `compute_radiation` sums over tau,
  before the final reduction. Visualize with `py_scripts/plot_debug_integrand.py <long|short> <mu> <nu>` (`mu >
  nu` resolved via `F^{nu mu} = -F^{mu nu}`; `mu == nu` rejected, since the diagonal is identically zero).
- `debug_exponent.dat` (`export_radiation_phase`): the phase factor `exp(i*(x[0]+R)*k)` common to every component
  in `debug_integrand.dat` (all 6 x long/short) — factored out into its own file rather than repeated 12x per row.
  Visualize with `py_scripts/plot_debug_exponent.py`.

Both of those plotting scripts, and `plot_electron_trajectory.py`, plot against `tau/T` rather than raw `tau` (`T
= 2*pi/omega`, the laser period) — `run_output_utils.get_laser_period(run_dir)` reads `laser_frequency` back out
of that specific run's own `config.cfg` (written by `IoUtils::copy_config_to_run_directory`), not the live repo
config, so the axis stays correct even if `config/coherent_thomson.cfg` has since changed. `plot_debug_integrand.py`
still plots against raw `tau`, not `tau/T`.

### Known gaps / TODOs worth knowing before touching related code

- **`LaguerreGaussLaser`'s transverse-mode physics still isn't physics-reviewed against the reference formula**
  (Allen, Beijersbergen, Spreeuw & Woerdman, Phys. Rev. A 45, 8185 (1992); Siegman, "Lasers", ch. 17).
  `complex_amplitude` (`laser_field.cpp`) returns `{amplitude, d/dx_loc, d/dy_loc}` via the product rule on
  `amplitude = prefactor * V * C_n * hypergeometric_val`; derivatives were verified against finite differences
  (~`1e-8` to `1e-14` relative error) across several `p`/`l` combinations. The radial profile uses
  `MathUtils::hypergeometric_1F1_neg_int_a(-p, n+1, u)` (`n = |l|`) via the identity `L_p^n(x) = binom(p+n,p) *
  1F1(-p; n+1; x)` — **not** `MathUtils::generalized_laguerre`, which remains implemented and tested via its own
  recurrence but is unused elsewhere. `Npn`'s normalization is a **deliberately custom (non-unit-power) constant**
  matching this project's own convention, not the standard `sqrt(2*p!/(pi*(p+|l|)!))` LG normalization. The
  azimuthal factor is built as the exact polynomial `(x_loc + i*sign(l)*y_loc)^|l|` (equal to `rho^|l| *
  exp(i*l*azimuth)`) so amplitude and derivatives stay smooth exactly on-axis instead of hitting a removable but
  awkward `1/rho` singularity. `x_loc`/`y_loc`/`z_loc` are the electron's position rotated into the laser's
  canonical frame, assuming the beam waist sits at that frame's `z=0` origin (tied to the `tau_0_traj=0.0`
  hardcode below) — check that the beam actually starts near the waist for whatever config you're using.
- **`get_faraday_tensor` builds genuine `Ez`/`Bz` from the `dx`/`dy` amplitude derivatives** (the `div(E)=0`/
  `div(B)=0` condition), consistent with `complex_amplitude`'s `exp(-i*phi)` carrier-sign convention (flipping one
  without the other breaks the `+i/k_wave` coefficient's sign). This is only a first-order paraxial construction
  (Lax et al. 1975; Davis 1979): the `div(E)` residual is `~1e-6` relative for the fundamental Gaussian but grows
  to `~1e-4`–`3e-2` for higher-order `p`/`|l|` modes (worst on/near axis) — an expected limitation, not a bug, but
  it does impart a small spurious force on off-axis electrons during RK4 trajectory integration for higher-order
  modes.
- **Polarization coefficients `zeta_1`/`zeta_2` are complex, not real.** `get_faraday_tensor` computes `E =
  Real(epsilon_1*zeta_1*amplitude + epsilon_2*zeta_2*amplitude)`. Linear polarization along `epsilon_1` is
  `zeta_1=(1,0)`, `zeta_2=(0,0)`; circular is `zeta_1=(1,0)`, `zeta_2=(0,1)`. **`create_laser`
  (`laser_factory.cpp`) normalizes `zeta_1`/`zeta_2` right after reading them from the config**, so `|zeta_1|^2 +
  |zeta_2|^2 = 1` always holds regardless of the raw values' sum (e.g. the circular convention above sums to `2`
  unnormalized) — without this, switching polarization would silently change the field's amplitude/intensity as a
  side effect of the polarization choice, instead of `a0` alone controlling field strength. Throws if `zeta_1`/
  `zeta_2` are both exactly zero.
- **`Radiation::compute_radiation` accumulates the full per-electron radiation physics into a packed
  representation.** For every `(tau, screen_point)` pair it builds the null vector `n0` from the electron-to-screen
  separation, then per frequency accumulates `amp_long`/`amp_short` times the six independent bivector terms
  `n^alpha u^beta - n^beta u^alpha` into a `Simulation::PackedRadiationField` (`ComplexBivector`, 6-element packed
  form) — real amplitude physics, just stored packed until reduction. Built in the **lab frame** (matching
  `laser_nx/ny/nz`), not the canonical frame `init_simulation_parameters` uses (see the `k1`/`p`/`n2` bullet
  below). `run_simulation` gives each thread its own accumulator, sums them once all threads join, then unpacks
  each summed `PackedFaraday` into a full 4x4 `Faraday` tensor via `MathUtils::unpack_bivector`. When the config
  key `print_field_in_canonical_frame` is `true` (default), the tensor is then rotated back into the laser's
  canonical frame via `rotate_tensor(inverse_rotation_tensor(laser.get_rotation_matrix()), ...)`, so plots don't
  depend on `laser_nx/ny/nz`; `false` leaves it in the lab frame instead. Either way the final tensors are scaled
  by `general_factor` (`= 1/(2*pi*c)`). `MathUtils::mirror_antisymmetric_in_place` is dead code (predates
  `unpack_bivector`, nothing calls it). `num_threads` (config key, default `0` = all hardware threads) controls
  `run_simulation`'s beam partitioning.
- **`compute_radiation`'s loop order and phase-factor evaluation are deliberately tuned for performance** — this
  is the dominant cost of a run. The loop nests screen point (`i_d`) outer, trajectory point (`i_tau`) inner,
  reversed from the naive order, so each screen point's tau-sum accumulates in cache/register-resident local
  tensors and the (potentially tens-of-MB) `field` array is written once per `(i_d, i_freq)` instead of once per
  `(i_tau, i_d, i_freq)` triple. Since `frequencies_list` is always built as an evenly-spaced (arithmetic
  progression) list — consecutive harmonics by default (see the `omega_min`/`omega_max` bullet below), or a
  linear scan when `dense_frequency_spectrum=true` — the phase factor `exp(i*phase*freq)` for `frequencies_list[0]`
  and for the constant spacing between entries are each computed via one `std::polar` call, and every other
  entry's phase factor follows from cheap complex multiplications (guarded by a per-electron
  `frequencies_are_evenly_spaced` check, falling back to per-frequency `std::polar` otherwise) — keep this in sync
  with the `omega_min`/`omega_max` gap below if you touch either.
  Measured ~2.6x wall-clock speedup from these two changes together on a fine-detector-grid config.
- **`tau_0_traj` is hardcoded to `0.0`** in `Simulation::init_simulation_parameters` rather than derived from the
  pulse's actual physical start — every electron starts its proper-time grid at `tau=0` regardless of the pulse's
  leading Gaussian wing / `laser_delay` shift. Flagged in-code with `// hardcoded, to be modified`.
- **`laser_delay` currently has no effect**: `LaserField`'s constructor unconditionally overwrites `delay` with
  `wing_sigma_cutoff * wing_sigma`, discarding the parsed `laser_delay`. Because of this,
  `export_field_heatmap_z0`'s snapshot time (`t = wing_sigma_cutoff * wing_sigma / omega`, computed in
  `main.cpp`) is deliberately keyed off `wing_sigma_cutoff`/`wing_sigma` rather than `laser_delay` — revisit both
  together if the delay bug is fixed. The heatmap's x/y window is also laser-type-dependent: `plane_wave` (no
  transverse profile) uses the configured `field_heatmap_x/y_min/max`; `laguerre_gauss` ignores those and uses
  `+-2 * laser_lg_w0` instead, since `w0` sets the mode's actual transverse scale.
- All electrons in the generated beam get `compute_trajectory` run on them (inside `run_simulation` via
  `compute_radiation`); `main.cpp` additionally runs it once more directly on the first `min(10,
  electron_beam.size())` electrons (indices `0..N-1`) purely to export their trajectories to `electron.dat`.
  `Particle::plot_particle_trajectory` now takes `(electron_beam, electron_indices, filepath)` and writes every
  selected electron's rows into one file tagged with a leading `electron_id` column (blank-line-separated blocks;
  `pandas.read_csv` skips blank lines automatically). `plot_electron_trajectory.py` groups by `electron_id` and
  encodes electron identity by color (fixed `tab10` order, never cycled past its 10 slots — matching the `main.cpp`
  cap) and the x/y/z (or p1/p2/p3) component by linestyle, since color is already spent on electron identity.
- **Each recorded `Electron::State` also carries the electron's exact 4-acceleration** (`du^mu/dtau =
  (q_0/m_0) F^{mu nu} u_nu`), exported as `electron.dat`'s trailing `a0 a1 a2 a3` columns (not currently
  plotted by `plot_electron_trajectory.py`). `Electron::update_state`'s RK4 stage `k1` is already the exact
  derivative at the state being stepped away from (`trajectory.back()` on entry), so it's reused to fill in
  that state's acceleration at no extra `Faraday`-tensor-evaluation cost; the very last state in a trajectory
  has no following `update_state()` call to supply its `k1`, so `Electron::compute_trajectory` fills it with
  one extra explicit `compute_derivative()` call after the loop — the only additional evaluation this scheme
  costs for the whole trajectory. Any `State` read back before `compute_trajectory()` has run holds a
  zero-initialized, not physically meaningful, acceleration.
- The whole beam is generated and held in memory upfront rather than per-thread/on-the-fly inside
  `run_simulation`; a deliberate temporary simplification until the radiation calculation is validated,
  with on-the-fly generation planned as a later memory optimization.
- **`PhysUtils::non_linear_Thomson_formula` (`phys_utils.hpp`) itself is flagged in-code as still needing review**
  (`// TO BE DOCUMENTED IN CLAUDE.md ; the agreement with the implementation should be checked`). It computes the
  lab-frame frequency of the `N`-th relativistic-Doppler harmonic radiated by a particle of four-momentum `p1`,
  seen along direction `n2`, for incident light of four-momentum `k1`: `omega2 = omega1 * N * contract(p1,
  n1)/contract(p1, n2)` where `omega1 = k1[0]*c` and `n1 = k1/k1[0]`; the return value is `omega2/c` (an `omega/c`
  convention, not raw angular frequency — see the `omega_min`/`omega_max` bullet below). It is exactly linear in
  `N`, which is what makes the `radiation_field.dat` `omega`-column normalization (see below) exact once both the
  harmonics and the fundamental share the same momentum argument.
- **`omega_min`/`omega_max` are only honored when `dense_frequency_spectrum` (config key, default `false`) is
  `true`.** By default `init_simulation_parameters` still builds `frequencies_list` as `N_harmonics` consecutive
  harmonics of the nonlinear Thomson formula (`PhysUtils::non_linear_Thomson_formula`) starting at harmonic index
  `N_harmonics_min` (config key; `1`, the fundamental, reproduces the old "first `N_harmonics` harmonics"
  behavior — `IoUtils::get_number_of_harmonics_min`), ignoring `omega_min`/`omega_max` entirely — appropriate for
  imaging over a whole detector screen, where you want the harmonic peaks' locations, not fine resolution between
  them. `N_harmonics_min` only shifts which harmonics are computed; it has no effect on `fundamental_frequency`
  (always the true, dressed-momentum fundamental at `N=1`, see the `k1`/`q`/`n2` bullet below), so
  `radiation_field.dat`'s exported `omega` column still reads e.g. `5.0` for the first row when
  `N_harmonics_min=5` — see that bullet's normalization note. `Radiation::compute_radiation`'s fast phase-factor
  recurrence (see the loop-order bullet above) was generalized to any evenly-spaced `frequencies_list` (constant
  step between consecutive entries, checked from `frequencies_list[0]`/`frequencies_list[1]`, not that the list
  start at harmonic `1`), so it stays in effect for `N_harmonics_min != 1` and, incidentally, for
  `dense_frequency_spectrum=true`'s linear scan too — verified to reproduce the direct per-frequency evaluation
  bit-for-bit in both cases. When `dense_frequency_spectrum=true`, `frequencies_list` is instead
  a plain linear scan of `N_omega` points from `omega_min` to `omega_max` (`IoUtils::get_omega_range`, real omega,
  divided by `c` to match `non_linear_Thomson_formula`'s own `omega/c` return convention) — `N_harmonics` and
  `N_omega` are deliberately separate config keys (`IoUtils::get_number_of_harmonics`/`get_number_of_frequencies`)
  so switching `dense_frequency_spectrum` doesn't silently reinterpret whichever count was already configured for
  the other mode. This mode is meant for probing the *shape* of a single Thomson line (its width, not just its
  peak position) at one screen point. It does indirectly depend on `fundamental_frequency` (see the `k1`/`q`/`n2`
  bullet below): `omega_min`/`omega_max` are each multiplied by `frequency_scaling_factor = fundamental_frequency *
  c / laser.get_omega()` before building the linear scan, so the configured range re-centers on the actual
  (possibly nonlinearly-shifted) fundamental instead of assuming the emitted frequency is an exact multiple of
  `laser_omega` — no new `Detector` type or plotting
  machinery: `Radiation::compute_radiation`/`Simulation::run_simulation`/`Radiation::plot_radiation_field` are
  already fully generic over both the frequency list's spacing and the detector's point count (every `Detector_2D`
  subclass already guards `N==1` and collapses to one exact point when its grid counts are set to `1`), so getting
  "one point, many frequencies" is purely a config choice — set the configured detector's grid to a single point
  (e.g. a `spherical` detector with `spherical_detector_N_theta=spherical_detector_N_phi=1` and matching
  `_min`/`_max` angles) and turn on `dense_frequency_spectrum`. If `dense_frequency_spectrum=true` but the detector
  has more than one point (a dense scan over a full imaging grid is a very expensive footgun, not something to fail
  silently into), `main.cpp` prints a non-fatal `std::cerr` error and then calls `detector->restrict_to_first_point()`
  (`Detector_2D::restrict_to_first_point`, `detector.hpp`) to collapse it to just index 0 before `run_simulation`
  runs — safe post-construction since `get_row_coordinate`/`get_col_coordinate` at index 0 only depend on
  per-detector-type spacing fixed at construction (`dx`/`dy`, `d_phi`, ...), not on the grid dimensions it
  overwrites. The detector geometry plots (`plot_detector`/`plot_detector_scatter`, called earlier in `main.cpp`)
  still reflect the full originally-configured grid, since the restriction happens only right before the actual
  (expensive) simulation. Visualize the result with
  `py_scripts/plot_point_spectrum.py` (`<long|short> <mu> <nu>`) — a line plot of `radiation_field.dat`'s
  `F^{mu nu}` vs. `omega`, the 1D counterpart of `plot_radiation_field.py`'s per-frequency 2D field-map PNGs (the
  wrong plot shape once frequency, not screen position, is the interesting axis).
- **`radiation_field.dat`'s exported `omega` column is in units of the fundamental, not raw atomic-unit omega.**
  `Radiation::plot_radiation_field` divides every `frequencies_list` entry by
  `simulation_parameters::fundamental_frequency` (`= PhysUtils::non_linear_Thomson_formula(k1, q, n2, 1)` — note
  the *dressed* momentum `q`, not bare `p`, see the `k1`/`q`/`n2` bullet below — computed once in
  `init_simulation_parameters` regardless of `dense_frequency_spectrum`) before writing it, so a
  value of `1.0` always means "the fundamental" and `3.0` always means "third harmonic" — including for the
  default harmonics mode, where the per-harmonic list is now also built from the dressed `q`, not bare `p`
  (`frequencies_list[i] = non_linear_Thomson_formula(k1, q, n2, N_harmonics_min + i)`). Since
  `non_linear_Thomson_formula` is exactly linear in its harmonic index `N` (`omega2 = omega1 * N *
  contract(p1,n1)/contract(p1,n2)`, `phys_utils.hpp`), and both calls now share the same `q`,
  `frequencies_list[i] / fundamental_frequency == N_harmonics_min + i` holds *exactly*, including for `a0 != 0` —
  this used to only hold approximately back when harmonics used bare `p` against a `q`-based
  `fundamental_frequency`; that mismatch was fixed by the "corrected the dressed momentum" /
  "changed the dressed momentum" commits. This normalizes against the
  *actual*, possibly Doppler-shifted fundamental for the configured observation direction/electron momentum, not
  the bare `laser_frequency` — the two coincide only when the beam's average momentum is zero and `a0` is
  negligible, so don't assume `omega=1.0` corresponds to `laser_frequency` for a moving beam, an off-axis detector
  direction, or a strong pulse. Both `plot_radiation_field.py` and `plot_point_spectrum.py` label this axis
  `$\omega/\omega_1$` accordingly.
- **`k1`, `p`/`q`, and `n2` in `init_simulation_parameters` are deliberately evaluated in the canonical frame**
  (laser along `Oz`), not the rotated lab frame — `non_linear_Thomson_formula` only combines its arguments through
  Minkowski contractions, invariant under a *common* rotation, so evaluating pre-rotation gives the same result
  without rotating anything. (Different from the radiated field itself, which *is* built in the lab frame by
  `compute_radiation` and rotated back explicitly — see above.) `k1` is fixed along canonical `Oz` scaled by
  `laser.get_omega() / c`; `n2` is the detector's own canonical-frame direction via
  `IoUtils::get_detector_direction_angles(config)`, not the electron's direction of motion — building `k1` from
  `laser.get_unity_n()` or `n2` from `average_px/py/pz` directly would mix frames and give wrong frequencies
  whenever the laser's configured direction isn't along `Oz`. **`q` is a ponderomotively-dressed momentum**, `q =
  p + (mc)^2*xi^2/(2*contract(p, k1)) * k1` with `xi = laser.get_a0()`, `mc = m_0*c` — used both to compute
  `fundamental_frequency` (accounting for the nonlinear frequency shift of the Thomson fundamental at high `a0`)
  and, as of the "changed/corrected the dressed momentum" commits, for the per-harmonic `frequencies_list`
  entries in the default (non-dense) mode too (previously those used bare `p`, which made the
  `frequencies_list[i]/fundamental_frequency` normalization only approximate for `a0 != 0` — see the
  `radiation_field.dat`'s exported `omega` column bullet above).
  The `(mc)^2` prefactor is required for unit consistency (`xi` is dimensionless, and `mc = 137.036` in these
  atomic units, so omitting it would make the correction term negligibly small regardless of `a0`). **RESOLVED**:
  the `2` denominator (`<a^2> = xi^2/2`) is confirmed correct — it fixed the backward-detector spectral-peak
  mismatch this formula was debugged against, and no longer needs re-deriving.
- **The detector has its own direction (`detector_direction_theta`/`detector_direction_phi`), independent of the
  laser's, but shares the laser's rotation.** `create_detector` passes both the laser's 4x4 `rotation_matrix` and
  the detector's own local direction into `Detector_2D`, which builds a 3x3 `local_rotation` orthogonal to that
  direction and composes the two per-point in `to_lab_frame` (local → shared canonical frame → lab frame) — so
  the detector rotates together with the laser instead of being locked to point exactly along it. `(0.0 pi, 0.0
  pi)` (the config default) points the detector straight along the laser.
- **`CircularDetector`'s radial grid is spaced in equal-*area* steps, not equal-distance**: `r_i = sqrt(R_min^2 +
  i*(R_max^2-R_min^2)/(N_R-1))`, so each ring encloses the same annular area despite fixed `N_phi` per ring —
  linear `r` spacing would make point density diverge as `1/r` near the center and collapse all `N_phi` points at
  `i=0` onto the origin when `R_min=0`. Anything reading `radiation_field.dat`'s circular-detector coordinates must
  reconstruct `r` from `i` via this formula, not assume linear spacing. `plot_radiation_field.py` mirrors it
  exactly and renders all three detector types as a `pcolormesh` over their native grid rather than a scatter —
  for `CircularDetector`/`SphericalDetector` this keeps the `phi=0`/`phi=2*pi` seam continuous for helical/vortex
  patterns; for `RectangularDetector` (grid already Cartesian-monotonic, so continuity isn't the issue) it instead
  avoids a coarse `Nx`/`Ny` grid rendering as sparse colored dots instead of a filled screen.
- **The beam cylinder's spatial axis is derived from the beam's own mean momentum direction, not fixed along
  canonical `Oz`.** `generate_cylinder_beam` computes a rotation from `average_px/py/pz` once per beam and applies
  it to each electron's local cylinder point before the shared laser rotation (falls back to identity if the mean
  momentum is numerically zero) — a beam whose mean momentum isn't along canonical `Oz` gets different
  per-electron realizations, not just different statistics, than a naive canonical-`Oz` cylinder would.
  `beam_center_x/y/z` (config keys, `lambda` units) offset each electron's raw cylindrical-Cartesian point
  (`generate_electron`, `electron_factory.cpp`) **before** both this mean-momentum rotation and the laser rotation
  — so the configured center is itself expressed in, and rotates along with, the same canonical frame as
  `beam_cylinder_radius`/`height`, not a fixed lab-frame shift applied after the beam is built. Translating
  before `beam_axis_rotation` rather than after is intentional, not an oversight: since that rotation is built to
  map local `+Oz` exactly onto the mean momentum direction, `beam_center_z` alone already lands as a pure
  translation along the beam's direction of motion (useful for staging the beam upstream so it reaches the laser
  focus region in sync with the pulse) — but this whole arrangement (including its interaction with the
  hardcoded `tau_0_traj=0.0` and the `laser_delay` bug below) is planned for revisiting, so don't assume the
  current scheme is final.
- **A rectangular detector and a spherical detector covering "the same" angular window will *not* generally show
  the same radiation pattern — this is real physics, not a bug.** `compute_radiation` uses the *exact*
  electron-to-screen distance `R` (no far-field linearization) in both amplitude falloff and phase. A rectangular
  screen's `R` varies across its area (`R = sqrt(D^2+x^2+y^2)`); a spherical screen's `R` is constant by
  construction — so the rectangular screen picks up an extra quadratic ("Fresnel") phase the spherical one never
  sees. The governing quantity is the **Fresnel number** `N_F = a^2/(D*lambda)` (`a` = rectangular half-width, `D`
  = distance): the two only converge once `N_F << 1` (far field); at `N_F` of order 1 or larger they visibly
  differ. The repo's example config (`x_min/x_max = -250/250 lambda`, `distance = 50000 lambda`) gives `N_F =
  1.25` — near field, not far field — so don't expect a rectangular-vs-spherical comparison to agree there without
  changing the config. **Common pitfall**: growing `rectangular_detector_distance` while holding the *angular*
  window fixed makes `N_F` grow, not shrink (`N_F = theta^2 * D/lambda`) — to actually approach the far field,
  hold the rectangular screen's *absolute* half-width fixed and increase `D` until `N_F << 1` (e.g. `D >=
  500000 lambda` for `N_F <= 0.125`), then match the spherical detector's `spherical_detector_theta_max` to the
  new, smaller angular window at that `D`. When comparing, use `Re(F01)`/`Im(F01)`, not `|F01|` — since the beam
  is tiny relative to the screen, the extra Fresnel phase is common to every electron's contribution and cancels
  out of the coherent sum's magnitude, only showing up in phase (`|F01|` correlates >0.9998 in both regimes and
  isn't a useful diagnostic here). The fundamental (`i_omega=0`) also isn't a useful test this close to the beam
  axis — its angular variation is ~6 orders of magnitude below its constant offset, below the double-precision
  noise floor.
- **`theory/angular_momentum_flux_density.md` derives the spectral angular-momentum flux density along `Oz`**,
  split into long-range/short-range field cross-terms (`(ll)`, `(ls)`, `(sl)`, plus the doc's implicitly-omitted
  `(ss)` needed for the four to sum back to the total) matching the long/short-range amplitude split
  `compute_radiation` already uses for the Faraday bivector. Implemented in Python only, in
  `py_scripts/plot_angular_momentum_flux.py`, entirely as post-processing of an existing run's
  `radiation_field.dat` — no C++ code computes or exports it. Restricted to `rectangular`/`circular` detectors
  (the formula assumes one flat transverse plane with a shared normal, which a `spherical` detector's points don't
  generally satisfy) and generalizes "`Oz`" to the detector's own local normal (`detector_direction_theta/phi`),
  reusing `Core::MathUtils::rotation_matrix_from_direction`'s exact three-case logic (reimplemented in Python) to
  rotate the exported Faraday tensor from the canonical/lab frame back into the detector's own local frame before
  applying the formula — see the script's module docstring for the full reasoning and its scope limits.
- **`py_scripts/plot_spherical_field_components.py` projects the exported Faraday tensor onto canonical-frame
  spherical components** (`E_r`/`E_theta`/`E_phi`/`B_r`/`B_theta`/`B_phi`), for `spherical`-detector runs only —
  again Python-only post-processing of `radiation_field.dat`, no C++ output involved. Unlike
  `plot_angular_momentum_flux.py`'s flat-screen restriction, a spherical detector is the natural fit here: every
  screen point already has its own observation direction, obtained by rotating its local (`theta`, `phi`) —
  `Core::Detector::SphericalDetector`'s own cone-point construction, generally relative to the detector's own axis
  via `detector_direction_theta/phi`, not necessarily canonical `Oz` — into the canonical frame via
  `get_detector_local_rotation`, then building the standard orthonormal `(r_hat, theta_hat, phi_hat)` basis at that
  direction and dotting it into the (already-canonical, or canonical-ized from lab frame via
  `get_field_to_canonical_rotation`) Cartesian `E`/`B`. Both new scripts share this rotation machinery, defined
  once in `plot_angular_momentum_flux.py` (`get_detector_local_rotation`, `get_laser_lab_rotation`,
  `get_field_to_canonical_rotation`, `extract_rotated_faraday_fields`) and imported by
  `plot_spherical_field_components.py` rather than duplicated.
  **Non-obvious physics sanity check surfaced while validating this script**: `compute_radiation`'s per-`(tau,
  screen point)` bivector term is built once from `n0`/`u` and shared, unscaled, between the long-range and
  short-range complex amplitude prefactors (`radiation.cpp`'s `add_bivector_term`) — so for a single contribution,
  `B` (built purely from the spatial `F^{jk}` block, i.e. an `n0 x u`-like construction) is *exactly* orthogonal to
  `n0` (triple product identity, `n0 . (n0 x u) = 0`), while `E` (built from the mixed `F^{i0}` block) is generally
  *not* exactly transverse to `n0` — unlike the textbook Lienard-Wiechert acceleration field. This was confirmed
  numerically on a 30-electron spherical-detector test run: `|B_r|/|B_theta|` came out ~`1e-5` (consistent with
  pure numerical/finite-beam-size residual) while `|E_r|/|E_theta|` was ~`0.1` (a real, nonzero effect of this
  formula, not a bug) — a useful regression check if this script (or `compute_radiation` itself) is ever modified:
  `B_r` should stay pinned near zero; `E_r` should not.
- **A spherical detector covering (close to) the full 4*pi sphere breaks a naive stereographic-projection
  heatmap plot** — a real bug, found and fixed while validating `plot_spherical_field_components.py` against a
  `spherical_detector_theta_max=1.0 pi` config. The stereographic formula both
  `Core::Detector::SphericalDetector::get_stereographic_projection` (`detector.cpp`, feeding
  `detector_stereographic.dat`) and its Python mirror `get_spherical_cell_edges`
  (`py_scripts/plot_radiation_field.py`) use — `rho = R*sin(theta)/(1+cos(theta))` — maps `theta=pi` (the pole
  antipodal to the projection's own reference pole) to `0/0`; this is an inherent property of stereographic
  projection (no single finite 2D chart can cover an entire sphere), not a rounding-error bug, so
  `detector_stereographic.dat` correctly (if silently) contains `nan` rows for any grid ring that reaches `theta=pi`
  — expected, not itself a defect. The actual bug was downstream: `plot_radiation_field.py`'s
  `plot_radiation_component` fed those `nan`-containing edges straight into `pcolormesh`, which raises
  (`x and y arguments to pcolormesh cannot have non-finite values`) instead of degrading gracefully. Fixed via two
  new functions in `plot_radiation_field.py`: `spherical_projection_is_well_defined` (checks whether
  `spherical_detector_theta_max` comes within a small tolerance of `pi`) and `get_spherical_plot_grid`, which uses
  the stereographic projection when that holds and otherwise falls back to a plain `(theta, phi)` rectangular map
  (own explicit cell corners, `aspect='auto'` since the axes are angles, not a shared length scale) — the standard
  way to visualize near-full-sphere angular data, since it has no polar singularity. Both
  `plot_radiation_component` (`plot_radiation_field.py`) and the heatmap branch of
  `plot_spherical_field_components.py` now call `get_spherical_plot_grid` instead of `get_spherical_cell_edges`
  directly, so both degrade gracefully together. `plot_detector_stereographic.py` was not touched — it renders
  with `plt.scatter`, which already drops `nan` points silently rather than crashing (so a full-sphere config just
  shows a scatter plot missing the exact-pole ring, not an error).
- **OPEN VALIDATION GAP: a single electron at rest at the origin, observed with a full-4*pi spherical detector,
  should reproduce the classical Thomson differential radiation distribution** (`dP/dOmega` proportional to
  `1+cos^2(theta)` for the repo's default circular polarization — `laser_zeta_1`/`zeta_2` giving `zeta_1=(1,0)`,
  `zeta_2=(0,1)` — independent of azimuthal `phi`, in the weak-field/small-`a0` dipole limit; `sin^2` about the
  polarization axis instead, for linear polarization). This is the standard benchmark any Thomson-scattering
  solver should reproduce for a single free electron (`beam_particle_count=1`, `beam_cylinder_radius`/
  `beam_cylinder_height`/`beam_center_x/y/z=0`, zero `average_p*`/`sigma_p*`), and as of this note it has been
  tried and does **not** come out matching that expected angular pattern — not yet root-caused (candidates worth
  checking first: whether `general_factor`/`amp_long`/`amp_short`'s normalization, the long+short recombination,
  or the observation-direction convention feeding `compute_radiation` actually reduce to the standard dipole
  formula in this weak-field/single-electron/large-`R` limit). Flagging this here since it's a fundamental
  correctness check that should hold before trusting the solver's absolute intensities/angular patterns for
  anything more complex (a coherent beam, higher `a0`, etc.) — revisit before relying on those.
- **`compute_radiation`'s per-`(tau, screen point, frequency)` integrand is factored into isolated
  `long_range_prefactor`/`short_range_prefactor`/`radiation_phase_argument` functions** (`radiation.cpp`, anonymous
  namespace), in preparation for the above gap: the current PREFACT formulas (multiplying the shared geometric
  bivector term `n0^alpha u^beta - n0^beta u^alpha`) were derived via integration by parts of the standard
  radiation integral and are suspected as the root cause, so they're kept swappable independently of the
  loop/phase machinery around them. `radiation_phase_argument(x, R) = x[0] + R` is the frequency-independent "rest
  of the exponent" in `PREFACT * exp(i*freq*phase_argument)` — frequency is applied by the caller, not this
  function, so it must stay frequency-independent for `compute_radiation`'s evenly-spaced-frequency phase-factor
  recurrence (see the loop-order bullet above) to remain valid regardless of how the PREFACT terms change.
  `debug/debug_radiation.cpp`'s `export_radiation_integrand`/`export_radiation_phase` deliberately duplicate (per
  the "Debug mode" section above) this same per-tau math inline rather than calling these functions — if the
  PREFACT formula is changed, that duplicated diagnostic path needs a matching update, or its cross-check against
  `compute_radiation`'s coherently-summed field (see "Debug mode" above) will silently go stale.
- **`theory/FT_Faraday_tensor-direct_and_simplified_forms.md` documents two independently-derived, analytically
  equivalent closed forms for the Fourier-transformed radiation field** (a "direct" form, straight FT of the
  Lienard-Wiechert field, vs. a "simplified" form via integration by parts on Jackson's compact form — the same
  integration-by-parts approach `compute_radiation`'s PREFACT formulas came from), each split into its own
  long-range/short-range pieces — **the two derivations' `F_l`/`F_s` splits do not agree term-by-term with each
  other or with `compute_radiation`'s; only each form's own `F_l+F_s` total is a valid cross-check target**.
  **Both forms are now implemented in `compute_radiation`**, selected via the `radiation_formula` config key
  (`"simplified"`, the default, or `"direct"`; `Simulation::run_simulation` throws at startup for any other
  value). The simplified form is what the anonymous-namespace `long_range_prefactor`/`short_range_prefactor`
  functions above compute (an explicit frequency factor on the long-range term, no acceleration needed); the
  direct form (`long_range_prefactor_direct`/`short_range_prefactor_direct`/`direct_long_range_tensor_term`,
  `radiation.cpp`) needs the electron's 4-acceleration (`Particle::Electron::State::acceleration`, already stored
  per trajectory point — see the acceleration bullet above) and has **no explicit frequency factor outside the
  shared phase** (`radiation_phase_argument` is identical between the two forms — confirmed against the theory
  doc's "Common notation" section — so the evenly-spaced-frequency phase-factor recurrence stays valid for both
  formulas unchanged). The direct form's short-range PREFACT also carries an extra `c^2` factor relative to the
  long-range term (baked into `short_range_prefactor_direct` itself, not into `run_simulation`'s uniform
  `general_factor`, since it's specific to this one formula — see the theory doc's Form 1 constants). `debug=true`
  only supports `radiation_formula=simplified`: `debug/debug_radiation.cpp`'s diagnostic still only reimplements
  the simplified form's per-tau math (unchanged by this addition), so `main.cpp` throws at startup if both are set
  together, rather than silently producing a debug export that can't cross-check a direct-formula run.
  **Fixing a pre-existing bug found while implementing this**: the short-range PREFACT's `n0·u` four-dot was
  previously computed as `contract(n0, u) * R` (an extra, erroneous factor of `R`, present in both
  `radiation.cpp` and its `debug_radiation.cpp` duplicate before this fix) rather than the bare `contract(n0, u)`
  the simplified-form formula calls for — this made the short-range term scale as `1/R^3` instead of the
  documented `1/R^2` near-field scaling. Fixed in both files; per the regime tested so far the short-range term
  doesn't contribute either way, so this fix alone does not resolve the OPEN VALIDATION GAP above — that
  remains open, and is now suspected to sit in the long-range term instead.
- **Cross-checking `radiation_formula=simplified` against `radiation_formula=direct` (same beam/laser/detector
  config, differing only in that key) found and fixed a genuine sign bug in the simplified form's `F_s` term, but
  also surfaced a second, still-open discrepancy in `F_l` that this fix does not touch.** Comparing the two forms'
  total field (`LR+SR`) component-by-component: 5 of the 6 independent Faraday-tensor components (`F01`, `F02`,
  `F12`, `F13`, `F23`) agree to ~4 significant figures between the two formulas, as expected; `F03` (equivalently
  `F30`) disagreed by roughly 1-2 orders of magnitude with an unrelated phase/sign, isolated to that one component
  pair. Because `SR` is ~6-7 orders of magnitude below `LR` in every config tested so far (confirmed again during
  this investigation), "only the total `F_l+F_s` needs to match, not the individual pieces" (the caveat in the
  theory doc above) is moot in practice: total is dominated entirely by `F_l`, so an `F_l`-vs-`F_l` mismatch shows
  up directly in the total, and an `F_s` bug of any size can't visibly move it.
  - **Sign bug (fixed).** Hand-rederiving the simplified form's integration-by-parts step (`theory/
    FT_Faraday_tensor-direct_and_simplified_forms.md`, Form 2) found `d/dτ(1/|R_0|) = +(n_{R_0}·u)/|R_0|^2`
    (three-vector dot), whereas the derivation this repo's `FT-simplified.tex` was built from had that term with
    the opposite sign — flipping the sign of the entire derived `F_s` term. Confirmed against the actual
    `FT-simplified.tex` source (since deleted from the repo — it was a working scratch file, not needed once the
    result was transcribed into the theory doc, and the bug is now documented there under "Sign bug (fixed)").
    Fixed in `radiation.cpp`'s `short_range_prefactor` (now returns the negated value) and in the theory doc's
    Form 2 `F_s` formula and its "shared tensor factor" note. **Verified this fix alone does not resolve the `F03`
    discrepancy above** — re-running both formulas after the fix still shows the same-sized `F03` mismatch, exactly
    because `SR`'s magnitude is too small relative to `LR` for its sign to matter to the total. So this is a real,
    independent bug (worth having fixed before relying on `F_s` on its own — see the angular-momentum-flux bullet
    below, where `(ls)`/`(sl)` cross-terms depend on `F_s`'s sign directly, unlike the total field), but it is not
    the explanation for the `F03` anomaly.
  - **`F03`/`F30` discrepancy — OPEN, not yet fixed.** Leading hypothesis: the simplified form's IBP derivation
    drops a boundary term, `[e^{ikφ(τ)} · n^{αβ}(τ) / (d_u(τ)·|R_0(τ)|)]`, evaluated at the τ-integration limits
    (`d_u = n_{R_0}·u`, the four-dot). The theory doc's Form 2 justifies dropping it as "contributes only as
    ω→0", but its actual suppression mechanism is the explicit `1/|R_0|` factor — i.e. it requires the electron to
    have travelled far enough that its distance to the (fixed, finite) detector point has grown large relative to
    its interaction-region value. For a realistic detector distance and pulse-timescale trajectory, the electron's
    total excursion is a tiny fraction of `|R_0|` regardless of how long the trajectory is integrated, so `1/|R_0|`
    barely changes — meaning this term may not actually be small at any practically affordable `wing_sigma_cutoff`,
    unlike the short-range term (whose own vanishing outside the pulse comes from a different, unrelated condition:
    the particle's velocity/acceleration settling to a constant, not from `1/|R_0|` growing). This was tested
    directly: increasing `wing_sigma_cutoff` did not shrink the `F03` gap, consistent with the boundary term being
    a fixed geometric omission rather than a finite-window artifact that integrates away. `F03` specifically (and
    not `F01`/`F02`/`F12`/`F13`/`F23`) is where this shows up because `F03` (`Ez`, the field component along the
    beam/propagation axis) is physically expected to be small (near-cancellation, by transversality) while the
    other components are not — so `F03` is the one component small enough for an uncounted boundary term to
    dominate its total, rather than being swamped by a genuinely large signal the way it is everywhere else.
    **Next step (not yet done): derive this boundary term's exact tensor form and evaluate it numerically at the
    electron's actual trajectory endpoints, to check directly whether it's actually the missing piece (zero
    residual once added back) or not** — before attempting any fix, since this is currently a hypothesis, not a
    confirmed root cause.
