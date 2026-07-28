# CoherentThomson

Simulates coherent (nonlinear) Thomson scattering of an intense laser pulse off a relativistic electron beam.

Given a config-driven laser pulse (plane wave or Laguerre-Gauss mode) and a detector screen (flat or spherical),
the solver generates an electron beam, integrates each electron's relativistic trajectory in the laser field, and
coherently sums the radiated field (as the antisymmetric Faraday bivector, per frequency/screen point) across the
whole beam. Output is exported as `.dat` files alongside a set of Python scripts for plotting the laser field,
electron trajectories, detector geometry, and the resulting radiation pattern.

## Build

Out-of-source build via CMake, C++20:

```
cmake -B build/ && cmake --build build/
```

The binary lands at `bin/coherent_thomson_solver` and requires a config file argument:

```
./bin/coherent_thomson_solver config/coherent_thomson.cfg
```

Or configure, build, and run in one step:

```
python3 py_scripts/compile_and_run.py [config_file]
```

Build types (`-DCMAKE_BUILD_TYPE=...`, default `Release`): `Release`, `Debug`, `RelWithDebInfo`, `MinSizeRel`.

## Visualizing output

Each plotting script takes no arguments — it locates the most recent `<output_folder>/YYYYMMDD_HHMMSS` run
automatically, so plotting always targets the last solver run:

```
python3 py_scripts/plot_laser_field.py
python3 py_scripts/plot_detector_stereographic.py
python3 py_scripts/plot_electron_trajectory.py
python3 py_scripts/plot_radiation_field.py
python3 py_scripts/plot_detector_scatter.py       # only if plot_detector_scatter=true in the config
python3 py_scripts/plot_electron_beam_scatter.py  # only if plot_beam_scatter=true in the config
python3 py_scripts/plot_field_heatmap_z0.py       # only if plot_field_heatmap=true in the config
```

## Configuration

`config/coherent_thomson.cfg` is a flat `key value [unit]` text format covering detector geometry, laser
frequency/envelope/direction/polarization, beam particle count/geometry, initial momentum distribution, and
radiation spectrum settings.

## Development

See [CLAUDE.md](CLAUDE.md) for the full architecture map, module/namespace layout, and a detailed list of known
physics gaps and TODOs in the current implementation.

Run `clang-format -i` on touched files before committing (Google style, 120 cols, 2-space indent, see
`.clang-format`). There is no test suite yet.

## License

[MIT](LICENSE)
