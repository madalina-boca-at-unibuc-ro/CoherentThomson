# 11: Scaled parameters of a run (implemented)

Spec: `theory/specs/11_scaled_parameters.md` (approved 2026-11-26 with the amendments of its §10).

- **Status:** done 2026-11-26; checkpoint (§8) reviewed by the user.
- **Commit:** see `git log`.

## Code

- `py_scripts/utils/scaled_parameters.py`: `scaled_parameters(run_dir)` → `ScaledParameters` (§4 + §10:
  `kappa`, per harmonic `band_low`, `band_high`, `N_fringes`, `omega_max`), `to_scaled_angles`, `dominant_component`,
  `check_screen_group`, `lg_peak_intensity`, report/table output, CLI. Parses `run_log.txt` by section and label
  (two-space top-level lines, so the incident screen's indented `distance` is not taken).
- `utils/run_output_utils.py`: `find_run_log_sections`, `write_run_log_section`, `read_harmonic_indices` moved from
  `plot_observables.py` (re-exported there under the old names).
- `plot_spectrum_map.peak_a2_fraction` delegates to `lg_peak_intensity` (one copy of the LG peak formula).

## Deviations and judgment calls

- Everything is read from `run_log.txt`; `config.cfg` turned out unnecessary (the log has the frequency list, N per
  frequency, and the dense window).
- Head-on check allows p = 0 (electrons at rest), as the lens run needs.
- `N_c` from `wing_sigma (phase)` / 2π; `None` (warning) with a flat top.
- Dense spectra: harmonic N is listed when its band overlaps the window; its `omega_max` is the window top capped at
  the boundary halfway to the next band. Without κ (no LG), round(ω/ω₁) ± 0.5.
- θ = x/Z in `to_scaled_angles`, as `plot_field`'s angular axes.
- `linewidth` changed 2026-11-26 (after the spec): $2\sqrt{2\ln2}/(2\pi\sqrt N N_c)\simeq0.375/(\sqrt N N_c)$, the
  relative FWHM of the spectral intensity (manuscript Eq. `pulse_linewidth`), replacing the spec's $1/(NN_c)$; it also
  changes `angular_detuning_over_linewidth` (×2.67 at N = 1, ×3.77 at N = 2).

## Checks (2026-11-26)

| Check | Result |
|---|---|
| Helper move, byte-identity | `run_log.txt` identical after `plot_observables` + `lg_decomposition` (run C) and `plot_spectrum_map` (a line run) |
| Run A / lattice (`20261029_132340`, `20261031_093730`) | Z/Z_* 2.808, θ_max/θ_* 5.342 |
| Run B / lattice (`…143739`, `…094539`) | 2.850, 8.773 |
| Lens run at rest (`20261030_111559`) | Z/Z_* 10.24 (notes: θ_max 3.9 θ_*, Fresnel term 0.23 rad) |
| p_z = 13.483 mc, a0 = 0.1 | D 27.003, k1/kL 725.536 (log 725.536) |
| Runs C, D (LG22) | σ −1, m −2, L_N −2 / −5, `minus` |
| Line run `10um_a0_1p0_N1` | θ_max/θ_* 16.0, Z/Z_* 2.85, k1/kL 486.11; aliasing margin 1.22 (warning: replica at 19.5 θ_*) |
| y-line (`spec15/s5_yline`), 2D rectangle (`~/output_tests/spec11/rect`) | θ_max from the line / inscribed half-width |
| `--write-log` twice | one section, rest unchanged |

Outputs: `~/Dropbox/work/AI-outputs/Claude/2026-11-26/spec11_*`.
