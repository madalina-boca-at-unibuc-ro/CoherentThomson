# Diagnostic plots of the helicity components E₊, E₋: spec

> **Status: APPROVED 2026-11-26** (diagnostics only; axes stay in angle until spec 11). Test runs: run C =
> `~/output/progress_report_20261029/runs_A-D_lattice/20261031_094054`, lens run =
> `~/output/progress_report_20261029/lens_gaussian_pulse/20261030_111559`. Open point §8: local basis opt-in.

## 1. Purpose

`plot_field.py` draws the 2×2 Re/Im/|·|/arg diagnostic panels of the Faraday components $F^{\mu\nu}$. The paper's
statements are about the helicity components: the selection rule, pure helicity, LG content. So the same diagnostic
is wanted for

$$E_\pm=\frac{E_x\mp iE_y}{\sqrt2},\qquad E_i=cF^{i0}\ \text{(total field, lab frame)},$$

which carry spin ±1 about +z. These are diagnostics. The paper figures (one complex-colour panel per case) are a
separate layer, specified later.

## 2. Definitions

- **Lab basis, the default.** The definition above, identical to `lg_decomposition.helicity_components`, which stays
  the only implementation (§3).
- **Local basis (`--local-basis`).** Project onto the transverse basis of the local direction
  $\hat n$ = (screen point)/|screen point|: $E_1=\mathbf E\cdot R\hat x$, $E_2=\mathbf E\cdot R\hat y$ with
  $R$ = `rotation_matrix_from_direction(n̂)`, and $E_\pm=(E_1\mp iE_2)/\sqrt2$. This is the projection
  `plot_point_spectrum_observables.py` uses.
  - It agrees with the lab basis to $O(\theta^2)$: about 10⁻⁸ at D = 27.
  - It matters for runs at rest, where $\theta_{\max}\approx0.18$ rad.
- **Total field only.** The long/short/boundary split has no meaning for a helicity projection.
- **Dominant component:** $E_\sigma$, with $\sigma=2\,\Im(\zeta_1^*\zeta_2)\,\mathrm{sign}(n_{L,z})$ the laser spin about
  +z (`plot_observables._read_laser_spin`). The plots label it, so a run of either sign reads correctly.
- **Power fraction** of each helicity at each frequency:
  $P_\pm=\sum_i w_i|E_{\pm,i}|^2/\sum_i w_i(|E_{+,i}|^2+|E_{-,i}|^2)$. The weights are uniform on the equal-area
  circular grid and on the rectangular grid. This matches `lg_decomposition.py`'s `helicity_fraction`.
- **Incident field (`--incident`).** `incident_field.dat` is in the laser's own frame, so the helicities are about
  the laser direction. Titles say so.

## 3. Shared field module (refactor first)

`plot_observables.py` imports `plot_field.py`, so `plot_field.py` cannot import the field helpers from it. A new
**leaf** module, `py_scripts/radiation/field_components.py`, takes:
- `get_total_field`, `get_E_total`, `get_B_total`, moved from `plot_observables.py` (`_get_total_field`,
  `_get_E_total`, `_get_B_total`);
- `helicity_components`, moved from `lg_decomposition.py`;
- `local_helicity_components(E, directions)`, new (§2);
- `read_laser_spin`, moved from `plot_observables.py`.

The old names are re-exported from their modules, so every existing import keeps working.

**Check:** `run_log.txt` sections from `plot_observables.py` and `lg_decomposition.py` are byte-identical before and
after the move, on one run.

## 4. Changes to the scripts

- **`plot_field.py`:** split `plot_radiation_component` into
  1. getting the complex field per frequency (the F column logic as now), and
  2. rendering the 2×2 figure (unchanged).

  Then add `plot_helicity_component(sign, filepath, field_transform=None, png_subdir=None, local_basis=False)`,
  which renders the same 2×2 figure.
  - **CLI:** `plot_field.py helicity <plus|minus> [run_folder] [--incident] [--local-basis]`.
  - **Panel titles:** `Re(E_+)` and so on.
  - **Figure title:** frequency, power fraction (e.g. "power fraction 3.1e−12"), "= E_σ (dominant)" or "(weak)", and
    "lab basis" or "local basis".
  - **Files:** `png_folder/radiation/emitted/radiation_field_total_Eplus_omega{i}.png` (and `Eminus`).
  - **No new colour or axis conventions:** Re and Im share ±|E|max, the phase is flat-shaded, and the axes are in
    angle. When `scaled_parameters.py` exists, the axes switch to $\theta/\theta_*$ (as for every screen plot).
- **`plot_field_carrier_removed.py`** (was `plot_processed.py radial_phase`): `plus|minus` applies the existing carrier
  removal $e^{-ik_N(R-Z)}$ to $E_\pm$; the default (all components) includes both. **This is the version to look at:** without it, Re and Im are dominated by
  fringes.
- **`plot_all_components.py`:** after the 24 F plots, also the two helicities (raw), plus
  `plot_field_carrier_removed.py plus|minus`.
- **`run_batch.sh` of the paper campaign** (`theory/specs/12_paper_campaign.md` §6) calls `plot_field_carrier_removed.py plus`
  and `… minus` for each run.

## 5. What not to do

- No helicity plots of the long, short or boundary terms.
- No second copy of the helicity definition or of the spin sign.
- No change to the existing F plots. Their PNGs must stay identical (§6).

## 6. Verification

1. **Refactor is a no-op:**
   - `plot_field.py total 0 1` and `plot_field_carrier_removed.py total 0 1` give byte-identical PNGs before and
     after, on one run (same matplotlib version);
   - the `run_log.txt` sections are identical (§3).
2. **Consistency with the decomposition:** on report run C (LG22, $a_0=0.1$; lattice `20261031_094054`), the power
   fractions in the titles equal `lg_decomposition.py`'s `helicity_fraction` to 10⁻¹². `E_-` is labelled dominant
   ($\sigma=-1$ there).
3. **Phase winding:** on the same run, the arg panel of `E_-` at N = 1 winds $l=-2$ around the axis, and $l=-5$ at
   N = 2. Read off the plot, and cross-checked by the azimuthal spectrum in `lg_decomposition.py`.
4. **Incident field:** `--incident` on an LG$_0^2$ run gives one helicity with fraction ≈ 1 and phase winding $l_L$.
5. **Local basis:** on a lens run at rest (`20261030_111559`), lab and local fractions differ at the expected
   $O(\theta^2)$. On a D = 27 run they agree to about 10⁻⁸.
6. **Either sign:** on a run with $\zeta_2=-i$ (e.g. `reference_config`, if run), `E_+` is labelled dominant.

## 7. Checkpoint

Implement §3, check §6.1, then §4. Stop and show the `E_-` and `E_+` figures (raw and carrier-removed) for run C, with
§6.2–6.3. Then update `PLOTTING_NOTES.md` (helicity plots, the shared module, lab against local basis) and
`CLAUDE.md`'s command list.

## 8. Open point

- **Whether the local basis should be the default** for runs at rest only (decided automatically from
  $\theta_{\max}$), or always opt-in. I suggest opt-in, with a note in the figure title when
  $\theta_{\max}>0.05$ rad and the lab basis is used.
