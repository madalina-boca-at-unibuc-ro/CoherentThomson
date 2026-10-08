# CoherentThomson

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

Given a config-driven laser pulse (plane wave or Laguerre-Gauss mode, always propagating along Oz) and a detector
screen (rectangular, circular, or spherical), the solver generates an electron beam, integrates each electron's
relativistic trajectory in the laser field, and coherently sums the radiated field (as the antisymmetric Faraday
bivector, per frequency/screen point) across the whole beam, in parallel over electrons. Output is exported as
`.dat` files plus a `run_log.txt` into a per-run folder `~/<output_folder>/YYYYMMDD_HHMMSS/`, alongside Python
scripts for plotting the laser field, electron trajectories, detector geometry, and the resulting radiation.

## Build

Out-of-source build via CMake (>= 3.22), C++20:

```
cmake -B build/ && cmake --build build/
cmake --preset release && cmake --build --preset release   # same, via preset
```

The build is always `Release` (`-O3 -march=native`): `CMakeLists.txt` overrides any other `CMAKE_BUILD_TYPE`.
To build another type, pass `-DCOHERENT_THOMSON_ALLOW_NON_RELEASE=ON` (or use the `debug` preset); the option stays
cached until you pass `=OFF`. Because of `-march=native`, build on the machine that runs the binary.

The binary lands at `bin/coherent_thomson_solver` and requires a config file argument:

```
./bin/coherent_thomson_solver config/config.cfg
```

Or configure, build, and run in one step:

```
python3 py_scripts/compile_and_run.py [config_file]
```

## Visualizing output

Every script plots the most recent run by default; an optional trailing argument selects a run folder instead.
PNGs are written to `<run>/png_folder/<module>/`.

```
python3 py_scripts/laser/plot_field.py
python3 py_scripts/laser/plot_heatmap_z0.py        # only if plot_field_heatmap=true
python3 py_scripts/detector/plot_stereographic.py
python3 py_scripts/detector/plot_scatter.py        # only if plot_detector_scatter=true
python3 py_scripts/particle/plot_trajectory.py
python3 py_scripts/particle/plot_beam_scatter.py   # only if plot_beam_scatter=true
python3 py_scripts/radiation/plot_field.py <long|short|boundary|total> <mu> <nu> [--incident] [run_folder]
python3 py_scripts/radiation/plot_all_components.py [--incident]
python3 py_scripts/radiation/plot_point_spectrum.py <long|short|boundary> <mu> <nu>   # dense_frequency_spectrum=true
python3 py_scripts/radiation/plot_spherical_components.py <long|short|boundary> <E|B> <r|theta|phi>  # spherical detector
python3 py_scripts/radiation/plot_observables.py [--incident] [run_folder]   # angular momentum / energy observables
```

`--incident` plots the incident laser's analytic field (`incident_field.dat`, rectangular/circular detectors only)
instead of the scattered field. Plotting conventions: [py_scripts/PLOTTING_NOTES.md](py_scripts/PLOTTING_NOTES.md).

## Configuration

Configs are flat `key value [unit]` text files covering detector geometry and direction, laser
frequency/envelope/polarization/mode, beam particle count/geometry, initial momentum distribution, radiation
formula and spectrum settings, and thread count. Provided configs:

- `config/config.cfg`: the default, for current experiments.
- `config/config_cross_check.cfg`: parameter-matched to the independent Python reference implementation.
- `config/config_initial_momentum.cfg`.

`output_folder` is taken relative to `$HOME`, so the same config works on every machine the repo is checked out on.

## Development

See [CLAUDE.md](CLAUDE.md) for the architecture map, module/namespace layout, and known physics gaps, and
[theory/](theory/) for the formulas and implementation notes.

Run `clang-format -i` on touched files before committing (Google style, 120 cols, 2-space indent, see
`.clang-format`). There is no test suite yet.

## License

[MIT](LICENSE)
