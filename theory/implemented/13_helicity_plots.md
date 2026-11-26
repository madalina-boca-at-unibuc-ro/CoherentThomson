# 13: Helicity-component plots E₊, E₋ (implemented)

Spec: `theory/specs/13_helicity_plots.md` (approved 2026-11-26: diagnostics only, axes stay in angle until spec 11 is
used by these plots; §8 open point: local basis opt-in).

- **Status:** done 2026-11-26; checkpoint (§7) shown to the user.
- **Commit:** code in `235a502` (committed by the user), docs in the commit after `0559d60`.

## Code

- `py_scripts/radiation/field_components.py` (new leaf of the radiation scripts): `get_total_field`, `get_E_total`,
  `get_B_total`, `helicity_components` (moved), `local_helicity_components` (new), `read_laser_spin` (moved).
  `plot_observables.py` and `lg_decomposition.py` re-export the old names.
- `plot_field.py`: `plot_radiation_component` split into `_read_field_file`, `_screen_grid`, `_png_dir`,
  `_render_2x2`; new `plot_helicity_component(sign, filepath, field_transform, png_subdir, local_basis, title_note)`
  and CLI `plot_field.py helicity <plus|minus> [run] [--incident] [--local-basis]`. Files
  `radiation_field_total_E{plus,minus}[_local]_omega{i}.png`.
- `plot_field_carrier_removed.py`: `plus|minus`, both by default, and `--local-basis`.
- `plot_all_components.py`: after the F plots, E± raw and (emitted field) carrier-removed.

## Deviations and judgment calls

- Power fraction with `lg_decomposition.py`'s screen-area weights (`get_screen_area_weights_au`) instead of the
  spec's "uniform": that is what makes it equal `helicity_fraction` exactly.
- `--local-basis` with `--incident` is refused (laser frame); local-basis files get `_local`.
- `--local-basis` added to `plot_field_carrier_removed.py` too (user's request, not in the spec).
- Spec 12's `run_batch.sh` does not exist yet: its calls wait for spec 12.

## Checks (2026-11-26; run C = `progress_report_20261029/runs_A-D_lattice/20261031_094054`, copied to
`~/output_tests/spec13/`)

| § | Check | Result |
|---|---|---|
| 6.1 | refactor no-op | 22 PNGs, `run_log.txt` sections, `lg_decomposition.dat` byte-identical |
| 6.2 | fractions vs `lg_decomposition.py` | E₊ 5.499e-13 (N = 1), 1.254e-12 (N = 2), equal; E₋ dominant (σ = −1) |
| 6.3 | winding of E₋ | l = −2 (N = 1), −5 (N = 2), azimuthal purity 1.000000 |
| 6.4 | `--incident` | E₊ dominant, fraction 1, about the laser direction; \|l_L\| = 2 lobes |
| 6.5 | local vs lab | run C: 8e-16; lens run at rest: weak fraction 5.5e-8 (lab) vs 3.3e-8 (local), N = 1 |
| 6.6 | either sign | no run with ζ₂ = −i; with it written into a copy of the config, σ = +1 about +z (label logic only) |

Paper-parameter test run (LG₀², a0 = 0.1, N_c = 30, D = 27, Z/Z_* = 2.81;
`reports/2026_thomson_paper/prelim/fig_maps/10um_a0_0p1_LG02.cfg`): E₊ fraction 4.4e-12 / 4.0e-12, l = −2 / −5.
Scripts: `~/Dropbox/work/AI-outputs/Claude/2026-11-26/spec13_checks.py`.
