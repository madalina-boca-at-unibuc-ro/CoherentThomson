# Spectral densities were 2π too small (fixed 2026-11-21)

## Summary

Every spectral density computed by the post-processing was too small by exactly `2π`: energy, Poynting flux, spin
and orbital AM densities and their fluxes. The solver's field (`radiation_field.dat`, `incident_field.dat`) was
always correct. Anything that is a ratio of two such densities is unchanged: AM per photon, flux/density = c,
helicity and LG fractions, purities, fidelities. Found while checking the post-processing against the manuscript
(`~/Dropbox/pub/ms/2026-Thomson/TeX/sections/02-model.tex`, Eqs. energy, am_density, am_flux), which has the correct
prefactors.

## The convention and where the factor was lost

The solver uses (`theory/FT_Faraday_tensor-direct_and_simplified_forms-v2.md`)

$$\tilde F(\omega,\mathbf x)=\frac1{2\pi}\int dt\,e^{i\omega t}F(t,\mathbf x),\qquad
F(t)=\int d\omega\,e^{-i\omega t}\tilde F(\omega).$$

The `1/(2π)` is in `general_factor` (`simulation.cpp`). The Thomson-dipole benchmark
(`theory/future_dev/tests/thomson_dipole_check.py`, the same convention) confirms the absolute amplitude to ~1e-4.

For real `A(t)`, `B(t)`, Parseval's theorem in this convention is

$$\int dt\,A_i(t)B_j(t)=2\pi\int_{-\infty}^{\infty}d\omega\,\tilde A_i\tilde B_j^*
=4\pi\int_0^{\infty}d\omega\,\Re[\tilde A_i\tilde B_j^*].$$

`numerical_calculation_of_angular_momentum.md` wrote the time integral as `2πδ(ω − ω′)` and then dropped the `2π` on
the next line, giving `2 Re[…]` instead of `4π Re[…]` for the one-sided density (Eq. bilinear-spectral-density).
Every formula in that note, and both scripts, were built on that line, so the factor is the same `2π` everywhere.

## Corrected formulas (circular screen; rectangular: `∂_φ → x∂_y − y∂_x`)

| Quantity | Before | Now |
|---|---|---|
| `du/dω` | `ε₀(\|E\|² + c²\|B\|²)` | `2π ε₀(\|E\|² + c²\|B\|²)` |
| `dP_z/dω` | `2 ε₀c² Re[E_x B_y* − E_y B_x*]` | `4π ε₀c² Re[…]` |
| `dS_z/dω` | `(4 ε₀/ω) Im[E_x* E_y]` | `(8π ε₀/ω) Im[…]` |
| `dL_z/dω` | `(2 ε₀/ω) Σ Im[E_i* ∂_φ E_i]` | `(4π ε₀/ω) Σ Im[…]` |
| `dΣ_zz/dω` | `(2 ε₀c²/ω) Im[B_x* E_x + B_y* E_y − B_z* E_z]` | `(4π ε₀c²/ω) Im[…]` |
| `dΛ_zz/dω` | `(2 ε₀c²/ω) Im[B_y* ∂_φE_x − B_x* ∂_φE_y + B_z* E_z]` | `(4π ε₀c²/ω) Im[…]` |

The photon number `dN_γ/dω = ∫dS (dP_z/dω)/(ħω)` scales with them.

## What changed

- `py_scripts/radiation/plot_observables.py`: the six prefactors and the docstring formulas.
- `py_scripts/radiation/plot_point_spectrum_observables.py`: the four prefactors and the docstring.
- `theory/numerical_calculation_of_angular_momentum.md`: the Parseval step, Eq. (bilinear-spectral-density) and
  every formula derived from it, with a dated note at the Parseval step.

No other script uses an absolute density: `lg_decomposition.py`, `summarize_runs.py` and `compare_lattice_random.py`
use ratios only. The C++ solver is unchanged.

## Check

`plot_observables.py` before (`HEAD`, `a738ee8`) and after, on copies of run `~/output/20261031_130825` (circular
screen, polar lattice, N = 1, 2):
- every value in "Screen-integrated observables" changes by a factor of 6.283185 (2π to the 7 printed digits);
- "Flux/density ratios" and "OAM-to-energy ratios" are unchanged, to every printed digit.

`plot_point_spectrum_observables.py` got the identical change. It compiles but has not been run, because no
`dense_frequency_spectrum` run was at hand.

## Consequences for existing results

- **`run_log.txt` of past runs** keeps the old absolute values in the sections written by `plot_observables.py` and
  `plot_point_spectrum_observables.py`. Re-run the scripts (they replace their own sections) before quoting an
  absolute energy, AM or photon number from those runs.
- **Scientific report (2026-10-29), manuscript plan blocks:** all quoted numbers are per photon, ratios or fractions,
  so none changes. The colour scales of the report's Figs. 3–5 (energy, spin and orbital flux maps, pasted from
  `plot_observables.py`'s PNGs by `make_figures.py`) are 2π too small; their shapes are unaffected. Regenerate them
  only if the report is reissued.
