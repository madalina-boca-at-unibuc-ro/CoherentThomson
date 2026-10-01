# Plotting notes

Conventions and gotchas of the scripts in `py_scripts/`. Usage and commands are in `CLAUDE.md`.

## Rendering

- `radiation/plot_field.py` and `radiation/plot_spherical_components.py` draw 2x2 Re/Im/Abs/Phase panels. Re and
  Im share the `±|F|max` scale.
- Rectangular, circular, and spherical detectors all render as `pcolormesh` on the native grid. On circular
  detectors Re/Im/Abs use `gouraud` shading on cell centers, while **Phase stays flat-shaded**, because
  interpolating across the ±π wrap is misleading.
- Colorbar alignment requires freezing the layout (`fig.canvas.draw()`, `set_layout_engine(None)`,
  `cbar.ax.set_axes_locator(None)`) before `set_position`. `constrained_layout` alone misaligns them.
- Heatmaps of LG runs add secondary `x/w0`, `y/w0` axes (`utils/w0_axes_utils.py`). That module is a leaf: it
  deliberately doesn't import `plot_field`, to avoid circular imports.

## Axes and units

- Screen axes are always in `lambda` (`plot_field.DETECTOR_AXES_UNIT`, matching `main.cpp`'s
  `detector_axes_unit`) for every detector type, whatever units the config gives. Each detector bound is converted
  with its own unit (`plot_field.read_length`), as the C++ factory does, so mixed units (e.g. `R_min` in `lambda`,
  `R_max` in `w0`) are fine. Never take one key's unit for another key's value.
- Spherical detectors near the full 4π: stereographic projection is singular at θ=π (`nan` in
  `detector_stereographic.dat` is expected), and `get_spherical_plot_grid` falls back to a `(θ, φ)` map.
- `particle/plot_trajectory.py` plots against `tau/T`, reading the laser period from that run's own `config.cfg`
  snapshot via `get_laser_period`. It draws one panel per four-vector component (sharing an axis hid real
  variation), with position in `lambda` and momentum in `mc`.

## `plot_observables.py`

- Computes the SAM/OAM/TAM densities `S_z`/`L_z`/`J_z`, the fluxes `Sigma_zz`/`Lambda_zz`/`Flux_tot`, and the
  energy `u`/`P_z`, per `theory/numerical_calculation_of_angular_momentum.md`.
- The `L_z` operator is `x∂y−y∂x` (rectangular, `np.gradient`) or a periodic `∂φ` (circular).
- The angular-momentum quantities divide by the true `omega` (not `omega/omega_1`), recovered from
  `fundamental_frequency` in `run_log.txt`.
- Sanity checks it prints: each flux/density ratio ≈ `±c` (holds to ~3e-5), and `Lambda_zz/P_z` ≈ `m/omega_N`.
- It appends its screen-integrated tables to `run_log.txt`, replacing only its own section by exact header match
  (`_find_run_log_sections`). Don't reintroduce a "truncate after first marker" approach.
