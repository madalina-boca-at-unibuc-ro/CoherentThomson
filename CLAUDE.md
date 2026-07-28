# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

End-to-end today (`src/app/main.cpp`): build a config-driven laser pulse and detector screen, generate an electron
beam, integrate one hardcoded electron's trajectory for plotting, then run the coherent radiation spectrum
calculation (`Simulation::run_simulation`) over the whole beam and export `.dat` files for plotting. The
per-electron physics `run_simulation` delegates to, `Radiation::compute_radiation`, now accumulates the full
antisymmetric Faraday bivector tensor (long-range and short-range amplitude terms, normalized) per
frequency/screen point — see "Known gaps" below for what's still missing (the frequency-range and pulse-start
caveats).

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
    (`plot_radiation_field`) following the same plottable-output pattern as the other modules. `compute_radiation`
    computes the full antisymmetric Faraday bivector; `run_simulation` mirrors and normalizes it — see "Known
    gaps" below for the remaining frequency-range and pulse-start caveats.
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
| `laser/` | `Core::Laser` | `LaserField` base (Gaussian-flat-top temporal envelope, direction/polarization — `zeta_1`/`zeta_2` are complex, e.g. `(1,0)`/`(0,1)` for circular — caches a 4x4 `rotation_matrix`, from which `epsilon_1`/`epsilon_2`/`unity_n` are derived as its columns; `get_faraday_tensor` is pure virtual) + `PlaneWaveLaser`/`LaguerreGaussLaser` derived types (each with its own private `complex_amplitude`, already scaled by `E0_c`; `PlaneWaveLaser`'s returns a single `Complex`, `LaguerreGaussLaser`'s returns `{amplitude, d/dx_loc, d/dy_loc}` — see "Known gaps"), `create_laser` (returns `std::unique_ptr<LaserField>`, dispatches on `laser_type`), `export_field_vs_phase`,
`export_field_heatmap_z0` (canonical-frame z=0 transverse snapshot — see "Known gaps") |
| `detector/` | `Core::Detector` | `Detector_2D` base + `FlatDetector`/`SphericalDetector` (each built orthogonal to its own canonical-frame direction, then rotated together with the laser via its shared 4x4 `rotation_matrix`), `create_detector` (takes the `LaserField`), `plot_detector` |
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

- **`LaguerreGaussLaser`'s transverse-mode physics lives in `complex_amplitude` (`laser_field.cpp`), not
  `get_faraday_tensor`, and only its amplitude is actually consumed so far — the mode is still not physics-reviewed
  against the reference formula (Allen, Beijersbergen, Spreeuw & Woerdman, Phys. Rev. A 45, 8185 (1992); Siegman,
  "Lasers", ch. 17), and its longitudinal field components aren't wired in yet.** `complex_amplitude` returns
  `std::tuple<Complex, Complex, Complex>` — `{amplitude, d(amplitude)/dx_loc, d(amplitude)/dy_loc}` — computed
  together since they share almost all of their work. `x_loc`/`y_loc`/`z_loc` are the electron's position rotated
  into the laser's own canonical frame (`dot3(x_mu, epsilon_1)`/`dot3(x_mu, epsilon_2)`/`dot3(x_mu, unity_n)` —
  `epsilon_1`/`epsilon_2`/`unity_n` are themselves the columns of the laser's 4x4 `rotation_matrix`, i.e. the
  canonical frame's axes expressed in the lab frame, the same rotation `Detector`/`Particle` also share), assuming
  the beam waist sits at that canonical frame's origin (`z=0`); since `tau_0_traj` is hardcoded to `0.0` (see
  below) and the electron beam's own placement is independent of this assumption, check that the beam actually
  starts near the waist for whatever config you're using. The radial profile (associated Laguerre polynomial via
  `MathUtils::generalized_laguerre`, using the identity `dL_p^n/du = -L_{p-1}^{n+1}(u)` for its derivative), Gouy
  phase, and wavefront-curvature phase match the same closed form documented here previously; the azimuthal factor
  `rho^|l| * exp(i*l*azimuth)` is instead built as the exact complex polynomial `(x_loc + i*sign(l)*y_loc)^|l|`
  (identically equal to it, since `rho*exp(i*sign(l)*azimuth) = x_loc + i*sign(l)*y_loc`), so both the amplitude
  and its `x`/`y` derivatives stay smooth exactly on-axis (`rho=0`) instead of going through polar-coordinate
  expressions with a removable but awkward `1/rho` singularity there. The analytic derivatives were checked
  against central finite differences at an off-axis point (~1e-10 relative error) and confirmed smooth near the
  axis for `l=2` in a one-off scratch test (not committed to the repo — temporarily made `complex_amplitude`
  `public` to call it directly, then reverted). **`get_faraday_tensor` currently only consumes `std::get<0>` (the
  amplitude) and is otherwise structurally identical to `PlaneWaveLaser::get_faraday_tensor`**, projecting it
  straight onto `epsilon_1`/`epsilon_2` — so the LG mode still has no genuine field components along the
  propagation direction. Wiring `d(amplitude)/dx_loc`/`d(amplitude)/dy_loc` into `E_z`/`B_z` via the
  `div(E)=0`/`div(B)=0` condition (flagged in-code above `get_faraday_tensor`) is the next planned step, not yet
  done. **The mode is also missing its overall normalization constant** — `R_val` in `complex_amplitude` is built
  from `generalized_laguerre` and the Gaussian radial falloff alone, with no `sqrt(2*p!/(pi*(p+|l|)!))`-type
  prefactor applied, so the LG amplitude is not normalized consistently across different `p`/`l`/`w0` combinations
  (a plain `l=0, p=0` mode happens to reduce to an unnormalized Gaussian, which is why this hasn't surfaced yet).
  Separately, `MathUtils::hypergeometric_1F1_neg_int_a`/`hypergeometric_1F1_neg_int_a_derivative`
  (`math_utils.hpp`) already implement Kummer's confluent hypergeometric function for the exact
  non-positive-integer-`a` case the associated Laguerre polynomials need (`L_p^k(x) = binom(p+k,p) *
  1F1(-p; k+1; x)`), but neither is actually called anywhere yet — `generalized_laguerre`'s own three-term
  recurrence is what's live. Re-expressing the radial profile through the hypergeometric form (or otherwise wiring
  these helpers in) is a possible follow-up, not just the normalization fix, before relying on
  `laser_type = laguerre_gauss` for production physics.
- **Polarization coefficients `zeta_1`/`zeta_2` are complex (`Core::MathUtils::Complex`), not real, and both laser
  types build their field the same way from them.** `PlaneWaveLaser::get_faraday_tensor` and
  `LaguerreGaussLaser::get_faraday_tensor` both compute `E = Real(epsilon_1*zeta_1*amplitude +
  epsilon_2*zeta_2*amplitude)`, where `amplitude` is that laser type's own `complex_amplitude` (already scaled by
  `E0_c = a0*omega*m_e/|e|` — previously missing from `PlaneWaveLaser`, now folded directly into
  `complex_amplitude`'s `exp(...)` so no caller needs to apply it separately). Linear polarization along
  `epsilon_1` is `zeta_1=(1,0)`, `zeta_2=(0,0)`; circular is `zeta_1=(1,0)`, `zeta_2=(0,1)`. The repo's default
  config now sets linear polarization (`zeta_1=(1,0)`, `zeta_2=(0,0)`) — a behavior change from the old real-valued
  default (`zeta_1=zeta_2=1.0`), which under the previous cos/sin-quadrature formula actually traced out circular
  polarization, not linear. See the "Config file format" section above for how these are read from the config.
- **`Radiation::compute_radiation` computes the full per-electron radiation physics, into a packed representation.**
  For every `(tau, screen_point)` pair, `compute_radiation` builds the null vector `n0` from the electron-to-screen
  separation (`x - detector.get_point(i_d)`), then for every frequency accumulates `amp_long`/`amp_short` times the
  six independent upper-triangle bivector terms `n^alpha u^beta - n^beta u^alpha` into a `MathUtils::ComplexBivector`
  (a 6-element `std::array<Complex, 6>`, one entry per independent Faraday-tensor component) which is added into
  `field.field[i_freq][i_d].long_range` / `.short_range`, where `field` is a `Simulation::PackedRadiationField` —
  this is real amplitude physics, not a bare phase factor, just stored in packed (not full 4x4 tensor) form because
  the diagonal/lower-triangle are never touched until reduction is complete. Because `n0`/`u` are built from the
  electron's actual lab-frame trajectory and screen position, the resulting Faraday tensor comes out in the lab
  frame (rotated to match the laser's configured `laser_nx/ny/nz` direction), not the canonical frame the rest of
  `init_simulation_parameters` works in (see the `k1`/`p`/`n2` bullet below). `run_simulation` (`simulation.cpp`)
  gives each thread its own private `PackedRadiationField` accumulator, sums the per-thread packed fields once all
  threads join, then per `(i_omega, i_screen)` on the final reduced field (not per trajectory point) unpacks each
  summed `PackedFaraday` into a full `Simulation::Faraday` (4x4 `ComplexFourTensor`, zero diagonal, antisymmetric
  lower triangle) via `MathUtils::unpack_bivector`, then — **only when the config key
  `print_field_in_canonical_frame` is `true` (the default)** — rotates it from the lab frame back into the laser's
  own canonical frame (laser along `Oz`) via `MathUtils::rotate_tensor(MathUtils::inverse_rotation_tensor(laser.get_rotation_matrix()),
  ...)` (so plots are read in the same canonical orientation regardless of `laser_nx/ny/nz` —
  `inverse_rotation_tensor` inverts a rotation `FourTensor` via its `ud()` mixed form, valid because a spatial
  rotation is orthogonal; `rotate_tensor` then applies the full `T'^{alpha beta} = R^alpha_mu R^beta_nu T^{mu nu}`
  transformation law, verified by round-tripping forward+inverse on a random antisymmetric tensor and by checking
  it reproduces `epsilon_1 = contract(rotation_matrix, pol_dir1)`'s convention on a canonical-frame E-field, and
  by confirming empirically that two otherwise-identical runs differing only in `laser_nx/ny/nz` agree on
  `radiation_field.dat` to floating-point precision when this flag is `true` and clearly diverge when it's
  `false`). Setting `print_field_in_canonical_frame` to `false` instead leaves the field exactly as computed, in
  the lab frame (rotated to match `laser_nx/ny/nz`, not the canonical frame `k1`/`p`/`n2` use) — useful for
  inspecting the raw lab-frame field rather than the canonical-frame one. Either way, the final step multiplies
  every `long_range`/`short_range` tensor by `general_factor` (the `1/(2*pi*c)` physical prefactor) via
  `FourTensor::operator*=`. `MathUtils::mirror_antisymmetric_in_place` (`math_utils.hpp`) still
  exists but is **currently dead code** — it predates `unpack_bivector` and nothing calls it anymore; don't assume
  it's on the live path. `run_simulation` itself is fully wired into `main.cpp` and does partition the beam across
  threads correctly; `num_threads` is read from the config's `num_threads` key (`IoUtils::get_num_threads`) and
  defaults to `0` (= max hardware threads) in `config/coherent_thomson.cfg`.
- **`compute_radiation`'s loop order and phase-factor evaluation are deliberately tuned for performance, not just
  correctness — this is the dominant cost of a run and worth understanding before touching it.** The loop nests
  screen point (`i_d`) outermost and trajectory point (`i_tau`) innermost, the reverse of the naive "for each tau,
  splat onto every screen point" order: every trajectory point still contributes to every `(screen point,
  frequency)` pair either way (so the total arithmetic is unchanged), but with `i_tau` innermost each screen
  point's full tau-sum accumulates into small local `ComplexFourTensor`s (`local_long`/`local_short`) that stay
  cache/register-resident, and the `field` array (tens of MB for a fine detector — `N_omega * N_screen` Faraday
  structs) is only written once per `(i_d, i_freq)` pair after the tau loop finishes, instead of once per `(i_tau,
  i_d, i_freq)` triple. Separately, since `frequencies_list[i]` is built (see the `omega_min`/`omega_max` bullet
  below) as exact integer harmonics of a single fundamental, the phase factor `exp(i*phase*frequencies_list[i])`
  is obtained from one `std::polar` (cos/sin) evaluation for the fundamental plus cheap complex multiplications
  for the higher harmonics, rather than a separate transcendental call per frequency. This is guarded by a
  `frequencies_are_harmonics` check done once per electron (cost O(N_freq), negligible) rather than assumed
  outright, so a future fix to the `omega_min`/`omega_max` gap that makes frequencies non-harmonic safely falls
  back to the direct per-frequency `std::polar` evaluation instead of silently computing the wrong phase — if you
  touch either that gap or this recurrence, keep the two in sync. On a config with a fine detector grid and many
  trajectory points, these two changes together gave a measured ~2.6x wall-clock speedup on one test machine
  (verified bit-identical-to-rounding against the pre-optimization loop order by diffing `radiation_field.dat`);
  the relative win from the loop-order change specifically scales with how much the `field` array exceeds the
  CPU's cache size, so it matters more on larger detector grids or smaller-cache machines than it did on the
  machine it was measured on.
- **`tau_0_traj` is hardcoded to `0.0`** in `Simulation::init_simulation_parameters` rather than derived from the
  pulse's actual physical start — every electron currently starts its proper-time grid at `tau=0` regardless of the
  pulse's leading Gaussian wing / `laser_delay` shift. Flagged in-code with `// hardcoded, to be modified`.
- **`laser_delay` currently has no effect**: `LaserField`'s constructor (`laser_field.cpp`) unconditionally
  overwrites `delay` with `wing_sigma_cutoff * wing_sigma`, discarding the parsed `laser_delay` value passed in
  from `laser_factory.cpp`. Flagged in-code as a TODO to fold into the `tau_0_traj` fix above rather than keep two
  independent pulse-start-time mechanisms. **Because of this, `Laser::export_field_heatmap_z0`'s snapshot time
  (`t = laser_wing_sigma_cutoff * laser_wing_sigma / omega`, computed in `main.cpp`, gated by `plot_field_heatmap`
  in the config) is deliberately built from `laser_wing_sigma_cutoff`/`laser_wing_sigma` directly, not from
  `laser_delay`** — this matches `envelope(phi)`'s actual switch point from the leading Gaussian wing to the flat
  plateau (`phi = delay`, the *internal*, overwritten value) exactly, regardless of what `laser_delay` is set to
  (verified by setting `laser_delay` to a deliberately mismatched value and confirming the exported snapshot phase
  stayed locked to `wing_sigma_cutoff * wing_sigma`). If the `laser_delay` bug above is ever fixed so `delay`
  genuinely follows the parsed `laser_delay`, this snapshot-time formula in `main.cpp` should be revisited too.
  **The heat map's x/y window is also laser-type-dependent, decided in `main.cpp` before calling
  `export_field_heatmap_z0`**: for `laser_type = plane_wave` (which has no transverse profile at all, so any
  window just shows a uniform color) it stays fixed at the `field_heatmap_x/y_min/max` config values; for
  `laser_type = laguerre_gauss` those config values are ignored and the window is instead sized to the mode's own
  waist, `+-2 * laser_lg_w0` (via `IoUtils::get_laser_lg_params`), since `w0` — not any fixed length — is what
  actually sets the mode's transverse scale. Verified by exporting one heat map per `laser_type` with a
  non-default `laser_lg_w0` and confirming the LG run's exported `x` column spans exactly `+-2 * w0` while the
  plane-wave run's still spans the configured `field_heatmap_x_min/x_max`.
- All electrons in the generated beam now get `compute_trajectory` run on them (inside `run_simulation` via
  `compute_radiation`); `main.cpp` additionally runs it once more directly on `electron_beam[2]` (falling back to
  `electron_beam[0]` if the beam has 2 or fewer electrons) purely to export its trajectory to `electron.dat` for
  plotting.
- The whole beam is generated and held in memory upfront rather than per-thread/on-the-fly inside
  `run_simulation`; this is a deliberate temporary simplification until the radiation calculation is validated,
  with on-the-fly generation planned as a later memory optimization.
- **`omega_min`/`omega_max` from the config are currently ignored.** `init_simulation_parameters` builds
  `frequencies_list` as the first `N_omega` harmonics of the nonlinear Thomson formula
  (`PhysUtils::non_linear_Thomson_formula`), evaluated at the same set of values for the entire screen, regardless
  of the configured frequency range — flagged both in-code and in `config/coherent_thomson.cfg`'s comment on
  `omega_max`.
- **`k1`, `p`, and `n2` in `init_simulation_parameters` are deliberately evaluated in the canonical frame (laser
  along `Oz`), not the rotated lab frame.** `non_linear_Thomson_formula` only ever combines its arguments through
  Minkowski contractions, which are invariant under a *common* rotation of all arguments — so evaluating everything
  pre-rotation gives the same result without needing to rotate anything. (This is a different situation from the
  radiated field itself, which *is* built in the lab frame by `compute_radiation` and has to be rotated back
  explicitly at the end of `run_simulation` — see the `compute_radiation`/`PackedRadiationField` bullet above.)
  `k1` is fixed along canonical `Oz`
  (`create_unit_light_like_vector<double>(0.0, 0.0)`) scaled by `laser.get_omega() / c`; `n2` is the detector's own
  canonical-frame direction via `IoUtils::get_detector_direction_angles(config)` (shared with
  `detector_factory.cpp`), not the electron's direction of motion. Building `k1` from `laser.get_unity_n()` or `n2`
  from `average_px/py/pz` directly would mix frames and give wrong frequencies whenever the laser's configured
  direction isn't along `Oz`. Because `n2` is independent of `p1`, a massive particle always has `p0 > |p·n̂|`, so
  `non_linear_Thomson_formula`'s division is well-defined for every momentum/direction — no "electron at rest"
  special case needed.
- **The detector has its own direction (`detector_direction_theta`/`detector_direction_phi`), independent of the
  laser's, but shares the laser's rotation.** `create_detector` (`detector_factory.cpp`) takes the `LaserField`,
  reads the detector's local direction (same canonical frame as the beam), and passes both the laser's 4x4
  `get_rotation_matrix()` and the detector's own `dir_x/dir_y/dir_z` into `Detector_2D`'s constructor
  (`detector.cpp`), which builds its own 3x3 `local_rotation = rotation_matrix_from_direction(dir_x, dir_y, dir_z)`
  orthogonal to that local direction and stores the laser's 4x4 matrix as `lab_rotation`. `to_lab_frame` then
  composes the two per-point rather than precombining them into one matrix: `rotate3d(local_rotation, ...)` into
  the shared canonical frame, then `contract(lab_rotation, ...)` into the lab frame — so the detector rotates
  together with the laser instead of being locked to point exactly along it. `(0.0 pi, 0.0 pi)` (the default in
  `coherent_thomson.cfg`) points the detector straight along the laser.
- **The beam cylinder's spatial axis is derived from the beam's own mean momentum direction, not fixed along
  canonical `Oz`.** `generate_cylinder_beam` (`electron_factory.cpp`) computes `beam_axis_rotation =
  MathUtils::rotation_matrix_from_direction(average_px, average_py, average_pz)` once per beam and applies it to
  each electron's local cylinder point (height axis originally along canonical `Oz`) before the shared laser
  rotation; falls back to the identity rotation when the mean momentum is (numerically) zero. A beam config whose
  mean momentum isn't along canonical `Oz` gets different per-electron realizations (not just different
  statistics) than a naive canonical-`Oz` cylinder would.
- **A flat detector and a spherical detector covering "the same" angular window will *not* generally show the same
  radiation pattern, and this is real physics, not a bug — worth knowing before treating one as a correctness
  check for the other.** `compute_radiation` uses the *exact* electron-to-screen distance `R` (not a linearized
  far-field approximation) in both the amplitude falloff (`amp_long_0`/`amp_short_0`, `radiation.cpp`) and the
  phase (`phase_base = x[0] + R`). For a source point near the origin — the beam is only `beam_cylinder_radius =
  2.5 lambda` across, negligible next to the detector — a flat screen's `R` varies across its area (`R = sqrt(D^2 +
  x^2 + y^2)`), while a spherical screen's `R` is exactly constant (`= detector_radius`) at every point by
  construction. So the flat screen picks up an extra quadratic ("Fresnel") phase term the spherical one never
  sees, of size `a^2/(2D)` at the screen edge (`a` = flat screen half-width, `D = detector_distance`). The
  governing quantity is the **Fresnel number** `N_F = a^2/(D*lambda)`: the two detectors' patterns only converge
  once `N_F << 1` (the Fraunhofer/far-field regime); when `N_F` is order 1 or larger (the Fresnel/near-field
  regime) they will visibly differ, and neither is "wrong" — it's the same exact-`R` physics evaluated on two
  different screen shapes.
  - The repo's own example numbers (`detector_x_min/x_max = -250/250 lambda`, `detector_distance = 50000 lambda`)
    give `N_F = 250^2 / 50000 = 1.25` — order 1, i.e. squarely in the *near*-field regime, not the far field.
    Concretely, the corner-to-center path difference is `~1.25 lambda` — over a full wavelength — which is enough
    to imprint visible extra spiral fringe rings on the flat screen (a helical/vortex phase combined with a
    Fresnel chirp is a classic way to get spiral fringes) that a spherical detector, having zero path-length error
    at any distance by construction, will not show. **At these settings, do not expect a flat-vs-spherical
    comparison (e.g. of `E_x`) to agree** — a visible mismatch here is the expected outcome, not a bug to chase.
  - **Common pitfall: growing `detector_distance` while keeping the *angular* window fixed does not approach the
    far field — it moves further from it.** If the flat screen's half-width `a` is scaled up proportionally with
    `D` to preserve a fixed angular window `theta ~= a/D` (which is what happens if you pick the spherical
    detector's `detector_theta_max` to match the flat screen's *current* angular window/area and then grow `D`),
    then `N_F = a^2/(D*lambda) = theta^2 * D/lambda` — this *grows* with `D`, not shrinks. Matching the two
    detectors' *area* or *angular window* at a new distance does not by itself put you in the far field.
  - **To actually construct a converging flat/spherical comparison**, hold the flat screen's *absolute* half-width
    `a` fixed (e.g. keep `a = 250 lambda`) and increase `D` until `N_F = a^2/(D*lambda) << 1` — e.g. `D >= 500000
    lambda` gets `N_F <= 0.125`, `D ~= 5000000 lambda` gets `N_F ~= 0.0125`. This necessarily *shrinks* the angular
    window (`theta_max = a/D`) as `D` grows, so the spherical detector's `detector_theta_max` must be set to match
    that new, smaller window at the new `D` — not the original window from the near-field config.
  - **Confirmed experimentally** in a past session, comparing a near-field config pair (flat/spherical detector,
    `N_F~1.25`) against a far-field pair (`N_F~0.0125`) with a one-off comparison script (none of these
    exploratory `.cfg`/`.py` files are checked into the repo — recreate them from the base config by adjusting
    `detector_distance`/`detector_x_min`/`detector_x_max`/`detector_theta_max` per the bullets above if you need to
    redo this check): at `N_F~1.25` the flat screen's `Re(F01)` showed a visible two-armed spiral fringe pattern
    that the spherical screen (a plain dipole-like pattern, no spiral) did not share (correlation ~0.21–0.49 across
    harmonics); at `N_F~0.0125` the two converged to the same smooth pattern (correlation ~0.985–0.9999 for the
    1st/2nd harmonics). **Compare `Re(F01)`/`Im(F01)`, not `|F01|`**: since the beam (`2.5 lambda`) is tiny next to
    the screen, the extra Fresnel phase is (to leading order) common to every electron's contribution at a given
    screen point, so it cancels out of the magnitude of the coherent sum and only shows up in the phase — `|F01|`
    correlates >0.9998 in *both* regimes and is not a useful diagnostic here. The fundamental (`i_omega=0`) doesn't
    show this convergence in either regime because, this close to the beam axis, its angular variation is ~6
    orders of magnitude below its constant offset — below the double-precision noise floor, not a real discrepancy
    — so it isn't a meaningful test at this `theta_max`; the 1st/2nd harmonics are the ones that actually confirm
    the crossover.
