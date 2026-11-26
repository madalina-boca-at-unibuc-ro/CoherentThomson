# Plotting notes

Conventions and gotchas of the scripts in `py_scripts/`. Usage and commands are in `AGENTS.md`.

## Rendering

- `radiation/plot_field.py` and `radiation/plot_spherical_components.py` draw 2x2 Re/Im/Abs/Phase panels. Re and
  Im share the `±|F|max` scale.
- Rectangular, circular, and spherical detectors all render as `pcolormesh` on the native grid. On circular
  detectors Re/Im/Abs use `gouraud` shading on cell centers, while **Phase stays flat-shaded**, because
  interpolating across the ±π wrap is misleading.
- Colorbar alignment requires freezing the layout (`fig.canvas.draw()`, `set_layout_engine(None)`,
  `cbar.ax.set_axes_locator(None)`) before `set_position`. `constrained_layout` alone misaligns them.
- The laser heatmap of LG runs adds secondary `x/w0`, `y/w0` axes (`utils/w0_axes_utils.py`). That module is a
  leaf: it deliberately doesn't import `plot_field`, to avoid circular imports. Radiation screen plots don't use
  it: at the screen `x/w0` grows with the distance.

## Axes and units

- Radiation screen heatmaps (`plot_field`, `plot_observables`, `plot_spherical_components`) are **plotted in
  angle** (`plot_field.get_angular_axes`): `x/D` for flat screens, `2x_proj/R` (= `2tan(θ/2)cosφ`) for the
  spherical stereographic grid, in rad/mrad/μrad picked from the screen's half-width. A screen grown with its
  distance then gets identical axes and panel sizes. Limits are set explicitly to the cell-center extent, so
  flat- and gouraud-shaded panels share one window.
- Screen *coordinates* are always in `lambda` (`plot_field.DETECTOR_AXES_UNIT`, matching `main.cpp`'s
  `detector_axes_unit`) for every detector type, whatever units the config gives; observables are computed and
  integrated in these units, and only converted to angle for plotting. Each detector bound is converted
  with its own unit (`plot_field.read_length`), as the C++ factory does, so mixed units (e.g. `R_min` in `lambda`,
  `R_max` in `w0`) are fine. Never take one key's unit for another key's value.
- Spherical detectors near the full 4π: stereographic projection is singular at θ=π (`nan` in
  `detector_stereographic.dat` is expected), and `get_spherical_plot_grid` falls back to a `(θ, φ)` map.
- `particle/plot_trajectory.py` plots against `tau/T`, reading the laser period from that run's own `config.cfg`
  snapshot via `get_laser_period`. It draws one panel per four-vector component (sharing an axis hid real
  variation), with position in `lambda` and momentum in `mc`.

## Helicity plots (spec 13)

- `radiation/field_components.py` holds the only copies of the total-field extraction (`get_E_total`, `get_B_total`),
  `helicity_components`, `local_helicity_components` and `read_laser_spin`; `plot_observables.py` and
  `lg_decomposition.py` re-export the old private names. `plot_field.py` imports it inside functions only:
  `faraday_frame_utils` (which it imports) imports `plot_field` at module level.
- `plot_field.plot_helicity_component(sign, ...)` renders E_+ or E_- in the 2x2 format through the same
  `_screen_grid`/`_render_2x2` as the F plots (the F PNGs were checked byte-identical after the split). Only the total
  field: a helicity of the long/short/boundary terms means nothing.
- Title: power fraction with `lg_decomposition.py`'s screen-area weights (equal to its `helicity_fraction`), "= E_σ
  (dominant)" or "(weak)" from the laser spin about +z (about the laser direction for `--incident`), and the basis.
- **Lab against local basis.** Lab: spin about +z everywhere (the AM fluxes and the selection rule are about z).
  Local (`--local-basis`): about each point's direction n̂ (`rotation_matrix_from_direction`), the photon's own
  helicity. They differ at O(θ²): 1e-15 at D = 27, a factor ~2 in the weak fraction (1e-8 level) for runs at rest.
  Opt-in; the title notes θ_max > 0.05 rad in the lab basis. Refused with `--incident`. Files get `_local`.

## Scaled units (`utils/scaled_parameters.py`, spec 11)

- The single source of θ_*, Z_*, Z/Z_*, θ_max/θ_*, D, N_c, κ and the per-harmonic quantities. It parses
  `run_log.txt` (what the solver used), so it is a leaf: it imports only `utils/run_output_utils.py`. Plot scripts
  convert angles with `to_scaled_angles` or `theta_star`; never recompute k₁ or θ_* in a script.
- θ_max: the farthest screen point from the axis for a line (rectangular `Nx` or `Ny` = 1, circular `N_phi` = 1), the
  inscribed half-width for a 2D rectangle, R_max for a circle. Harmonics: the log's per-frequency N (spec 16), or for a
  dense spectrum the bands that overlap the window; the replica and the angular detuning use each harmonic's largest
  ω/ω₁.
- `peak_a2_fraction` (`plot_spectrum_map.py`) wraps its `lg_peak_intensity`; the run-log section writer
  (`write_run_log_section`) and `read_harmonic_indices` live in `utils/run_output_utils.py`, re-exported from
  `plot_observables.py` under the old names.

## `plot_observables.py`

- Computes the SAM/OAM/TAM densities `S_z`/`L_z`/`J_z`, the fluxes `Sigma_zz`/`Lambda_zz`/`Flux_tot`, and the
  energy `u`/`P_z`, per `theory/physics/angular_momentum_numerics.md`.
- The `L_z` operator is `x∂y−y∂x` (rectangular, `np.gradient`) or a spectral (FFT) periodic `∂φ` (circular; exact for `|l| < (N_phi−1)/2`,
  unlike the centered difference used before 2026-10-29, which was 4% low at `l = 5`), taken on the
  screen's **local** grid and multiplied by `lz_sign = det(R_⊥)` (−1 for a screen facing −Oz, whose local frame is
  lab `(x, −y, −z)`), so every observable is about lab +z. `E`/`B` are lab components; don't drop `lz_sign`.
  `--incident` needs no correction (`lz_sign = 1`): `incident_field.dat` is entirely in the laser's own frame, for
  any laser or detector direction. Its screen (distance `Z`, beam-sized grid) is not the detector's: scripts read
  it from the run's `incident_config.cfg` through `utils/run_output_utils.field_config_path`, never from
  `config.cfg` directly (older runs without that file: the waist plane on the detector's grid).
- The theoretical `Lambda_zz/P_z`: `(N(m+σ)−σ)/omega_N` for the emitted harmonics (`σ = 2 Im(ζ₁*ζ₂)`, the laser
  spin; `N m/omega_N` for linear), with `m`, `σ` times `sign(n_laser,z)` (about lab +z); `m/omega` with
  `m = laser_lg_l` for `--incident`.
- The angular-momentum quantities divide by the true `omega` (not `omega/omega_1`), recovered from
  `fundamental_frequency` in `run_log.txt`.
- Sanity checks it prints: each flux/density ratio ≈ `±c` (holds to ~3e-5), and `Lambda_zz/P_z` ≈ the theory column.
- It appends its screen-integrated tables to `run_log.txt`, replacing only its own section by exact header match
  (`_find_run_log_sections`). Don't reintroduce a "truncate after first marker" approach.
