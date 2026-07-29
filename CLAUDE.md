# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

End-to-end today (`src/app/main.cpp`): build a config-driven laser pulse and detector screen, generate an electron
beam, integrate one hardcoded electron's trajectory for plotting, then run the coherent radiation spectrum
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
```

Build types (`-DCMAKE_BUILD_TYPE=...`, default `Release`): `Release` (`-O3 -march=native -mtune=native`), `Debug`
(`-g`), `RelWithDebInfo`, `MinSizeRel`. C++20. `-Wall -Wextra -Wpedantic` on all types except
`RelWithDebInfo`/`MinSizeRel`.

Run `clang-format -i` on touched files before committing (`.clang-format` at repo root: Google style, 120 cols,
2-space indent).

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
| `particle/` | `Core::Particle` | `Electron` (RK4 Lorentz-force integrator, optional trajectory recording), `generate_cylinder_beam`/`generate_electron`, `plot_particle_trajectory` |
| `laser/` | `Core::Laser` | `LaserField` base (Gaussian-flat-top temporal envelope, direction/polarization — `zeta_1`/`zeta_2` are complex, e.g. `(1,0)`/`(0,1)` for circular — caches a 4x4 `rotation_matrix`, from which `epsilon_1`/`epsilon_2`/`unity_n` are derived as its columns; `get_faraday_tensor` is a single non-virtual implementation shared by every derived type, built on the pure-virtual `complex_amplitude` customization point) + `PlaneWaveLaser`/`LaguerreGaussLaser` derived types (each implementing `complex_amplitude`, already scaled by `E0_c`, returning `std::tuple<Complex, Complex, Complex>` — `{amplitude, d/dx_loc, d/dy_loc}`; `PlaneWaveLaser` has no transverse profile so its derivatives are always `{0, 0}` — see "Known gaps"), `create_laser` (returns `std::unique_ptr<LaserField>`, dispatches on `laser_type`), `export_field_vs_phase`, `export_field_heatmap_z0` (canonical-frame z=0 transverse snapshot — see "Known gaps") |
| `detector/` | `Core::Detector` | `Detector_2D` base + `RectangularDetector`/`SphericalDetector`/`CircularDetector` (each built orthogonal to its own canonical-frame direction, then rotated together with the laser via its shared 4x4 `rotation_matrix`), `create_detector` (takes the `LaserField`), `plot_detector` |
| `simulation/` | `Core::Simulation` | `init_simulation_parameters`, `Faraday`/`RadiationField` (full 4x4 tensor, post-reduction) + `PackedFaraday`/`PackedRadiationField` (6-element packed bivector, accumulation-time) + `run_simulation` (multithreaded, partitions beam across `num_threads`) |
| `radiation/` | `Core::Radiation` | `compute_radiation` — one electron's contribution to `Simulation::PackedRadiationField` (packed antisymmetric Faraday bivector, long/short-range amplitudes, summed over trajectory points/screen points/frequencies); `plot_radiation_field` exporter |

### Config file format

`config/coherent_thomson.cfg` is a flat `key value [unit]` text format covering detector geometry, laser
frequency/envelope/direction/polarization, beam particle count/geometry, initial momentum distribution, radiation
spectrum range, and a trajectory-print-frequency knob. Numeric values may carry a unit suffix (`lambda`, `pi`,
`mc`, `cycles_adim`, ...) resolved by `IoUtils::convert_unit_to_number` against the laser's own
wavelength/frequency — check that function before adding a new unit keyword. A key with a blank value (nothing but
a trailing comment) is dropped entirely by the parser rather than stored empty, so `config.at(...)` throws
`std::out_of_range` for it — always give every key a real value, even a placeholder.

One exception to the single-number-plus-unit convention: the laser's polarization coefficients are complex
(`Core::MathUtils::Complex`), so each is split into two plain-number keys instead of one `key value [unit]` line —
`laser_zeta_1_re`/`laser_zeta_1_im` and `laser_zeta_2_re`/`laser_zeta_2_im` — read via
`IoUtils::get_complex_config_value`/`get_laser_zeta` (`io_utils.hpp`), not `convert_unit_to_number`.

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
  `(i_tau, i_d, i_freq)` triple. Since `frequencies_list` is built as integer harmonics of one fundamental, the
  phase factor `exp(i*phase*freq)` is computed via one `std::polar` call plus cheap complex multiplications for
  higher harmonics (guarded by a per-electron `frequencies_are_harmonics` check, falling back to per-frequency
  `std::polar` otherwise) — keep this in sync with the `omega_min`/`omega_max` gap below if you touch either.
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
  `compute_radiation`); `main.cpp` additionally runs it once more directly on `electron_beam[2]` (falling back to
  `electron_beam[0]` if the beam has 2 or fewer electrons) purely to export its trajectory to `electron.dat`.
- The whole beam is generated and held in memory upfront rather than per-thread/on-the-fly inside
  `run_simulation`; a deliberate temporary simplification until the radiation calculation is validated,
  with on-the-fly generation planned as a later memory optimization.
- **`omega_min`/`omega_max` from the config are currently ignored.** `init_simulation_parameters` builds
  `frequencies_list` as the first `N_omega` harmonics of the nonlinear Thomson formula
  (`PhysUtils::non_linear_Thomson_formula`), regardless of the configured frequency range — flagged both in-code
  and in `config/coherent_thomson.cfg`'s comment on `omega_max`.
- **`k1`, `p`, and `n2` in `init_simulation_parameters` are deliberately evaluated in the canonical frame** (laser
  along `Oz`), not the rotated lab frame — `non_linear_Thomson_formula` only combines its arguments through
  Minkowski contractions, invariant under a *common* rotation, so evaluating pre-rotation gives the same result
  without rotating anything. (Different from the radiated field itself, which *is* built in the lab frame by
  `compute_radiation` and rotated back explicitly — see above.) `k1` is fixed along canonical `Oz` scaled by
  `laser.get_omega() / c`; `n2` is the detector's own canonical-frame direction via
  `IoUtils::get_detector_direction_angles(config)`, not the electron's direction of motion — building `k1` from
  `laser.get_unity_n()` or `n2` from `average_px/py/pz` directly would mix frames and give wrong frequencies
  whenever the laser's configured direction isn't along `Oz`.
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
  exactly and renders `CircularDetector`/`SphericalDetector` as a `pcolormesh` over the native grid (not a
  scatter), so the `phi=0`/`phi=2*pi` seam stays continuous for helical/vortex patterns; `RectangularDetector` is
  the only type still rendered as a scatter, since its grid is already Cartesian-monotonic.
- **The beam cylinder's spatial axis is derived from the beam's own mean momentum direction, not fixed along
  canonical `Oz`.** `generate_cylinder_beam` computes a rotation from `average_px/py/pz` once per beam and applies
  it to each electron's local cylinder point before the shared laser rotation (falls back to identity if the mean
  momentum is numerically zero) — a beam whose mean momentum isn't along canonical `Oz` gets different
  per-electron realizations, not just different statistics, than a naive canonical-`Oz` cylinder would.
  `beam_center_x/y/z` (config keys, `lambda` units) offset each electron's raw cylindrical-Cartesian point
  (`generate_electron`, `electron_factory.cpp`) **before** both this mean-momentum rotation and the laser rotation
  — so the configured center is itself expressed in, and rotates along with, the same canonical frame as
  `beam_cylinder_radius`/`height`, not a fixed lab-frame shift applied after the beam is built.
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
</content>
