# Scaled parameters of a run: spec

> **Status: DRAFT, NOT APPROVED FOR IMPLEMENTATION.** Written 2026-11-21. Discuss first; implement only after the user
> approves, then stop at the checkpoint (§8) for review.

## 1. Purpose and scope

The configs stay as they are: physical inputs in a.u., SI, `lambda` or `w0`. The scaling the paper claims is tested by
giving the solver different physical numbers (10 µm against 800 nm), so the inputs must not assume it.

This spec adds **one Python module that computes, from a finished run folder, every scaled quantity the paper uses**.
The plots and the paper-figure layer (separate spec, later) take their units from it, and nowhere else.

- **In scope:** the module; a command-line entry that prints the quantities, writes a `run_log.txt` section and
  tabulates several runs (the Supplement S4 table of all runs); the consistency and validity checks in §5.
- **Out of scope:**
  - any change to the solver, the configs or the output files;
  - plot changes (the figure spec);
  - collective parameters ($\zeta$), which need the physical charge, not part of a run.

## 2. Conventions (they must match the manuscript)

- **Head-on geometry:** laser along −z, beam along +z, screen normal +z. The module checks this (§5) and computes the
  scaled quantities only for it.
- **Doppler factor:** $D=(p^0-\hat{\mathbf n}_L\cdot\mathbf p)/mc$, the light-front momentum of the mean initial
  momentum. Head-on, this is $D=\gamma(1+\beta)$.
- **Reference dressing:** $A_{\rm ref}=a_0^2/2$, the solver's (`PhysUtils::dressed_momentum`). If the solver later gets
  an $A_{\rm ref}$ key, read it from `run_log.txt` instead.
- **Reference frequencies:** $k_{1,\rm ref}=k_LD^2/(1+A_{\rm ref})$ and $k_{N,\rm ref}=N\,k_{1,\rm ref}$. These are the
  solver's `frequencies_list`, and `radiation_field.dat`'s ω column is $\omega/\omega_{1,\rm ref}$.
- **Screen scales, defined at N = 1 for every harmonic** (decided 2026-11-21; the manuscript's Sect. III must say so):
  $$Z_*=\tfrac12 k_{1,\rm ref}w_0^2,\qquad \theta_*=\frac{2}{k_{1,\rm ref}w_0},\qquad \theta_*Z_*=w_0 .$$
  Angles are always given as $\theta/\theta_*$ and waist positions as $z_0/Z_*$, with the same $\theta_*$ and $Z_*$ for
  every N.
- **Pulse length:** $N_c=\tau_L/2\pi$ = `laser_wing_sigma` in cycles, for the pure Gaussian envelope
  (`laser_flat_duration 0`).
- **Helicity and AM about lab +z** (as `plot_observables.py`):
  - $\sigma=2\,\Im(\zeta_1^*\zeta_2)\,\mathrm{sign}(n_{L,z})$ and $m=l_L\,\mathrm{sign}(n_{L,z})$;
  - $J_N=N(m+\sigma)$ and $L_N=J_N-\sigma$;
  - the dominant emitted component is $E_\sigma$: `minus` for $\sigma=-1$, `plus` for $\sigma=+1$.

## 3. Inputs

- **`run_log.txt` is the primary source.** It holds what the solver actually used, not only what was requested: the
  lattice spacing after adjustment, the LG-table radius, the particle count, `fundamental_frequency`, $\hat{\mathbf n}$
  and $\zeta_{1,2}$ after normalization.
- **`config.cfg`** only for what the log lacks (e.g. `N_harmonics`, `trajectory_NT`), read with the existing
  `read_config_value`.
- **Parsing:** read the log by label, through regular expressions on its `key  value unit` lines (a.u. column).
- **Missing fields:** a missing field gives `None` and a warning, not an exception. Older runs lack some lines, and
  random beams have no lattice spacing.
- **Detectors:** circular and rectangular screens are supported. For a rectangular one, $R_{\max}$ is the inscribed
  half-width (the convention `lg_decomposition.py` uses). Spherical screens: the module reports the non-screen
  quantities and sets the screen ones to `None`.

## 4. Quantities

All lengths are converted internally to units of $\lambda_L$. The physical values are kept for the "physical" block.

**Physics parameters** (the inputs the scaling says matter)

| Name | Definition |
|---|---|
| `a0`, `N_c`, `a0sq_Nc` | $a_0$, $N_c$, $a_0^2N_c$ |
| `p`, `l`, `sigma`, `m` | LG indices (laser frame) and lab $\sigma$, $m$ |
| `aligned` | `sign(m) == sign(sigma)` (`None` for $m=0$) |
| `w0_over_lambda` | $w_0/\lambda_L$ |
| `R_over_w0`, `d_over_w0`, `n_electrons` | sheet radius, lattice spacing (actual), sites |
| `D`, `gamma`, `beta` | from `average_p` |
| `A_ref` | $a_0^2/2$ |

**Screen and frequency scales**

| Name | Definition |
|---|---|
| `k1_over_kL` | $D^2/(1+A_{\rm ref})$; compared with the log's `fundamental / laser omega` (§5) |
| `lambda_1` | $\lambda_L/(k_1/k_L)$ |
| `Z_star`, `theta_star` | as in §2 (in $\lambda_L$ and rad) |
| `Z_over_Zstar` | screen distance $/Z_*$ |
| `theta_max`, `theta_max_over_thetastar` | $\arctan(R_{\max}/Z)$, and $/\theta_*$ |
| `harmonics` | the list of N evaluated, from `N_harmonics_min` and `N_harmonics` |

**Per harmonic N** (one row per N)

| Name | Definition |
|---|---|
| `J_N`, `L_N`, `component` | expected AM per photon; `plus` or `minus` |
| `linewidth` | $1/(NN_c)$ |
| `replica_over_thetastar` | lattice grating replica at $\pi w_0/(Nd)$ in units of $\theta_*$ |
| `aliasing_margin` | `replica_over_thetastar / theta_max_over_thetastar` (must be > 1; warn below 1.25) |
| `angular_detuning` | $(\theta_{\max}/\theta_*)^2\,/\,[(k_1/k_L)(k_Lw_0)^2]$, the relative shift of the line at the screen edge |
| `angular_detuning_over_linewidth` | the previous divided by `linewidth` |

**Validity diagnostics** (the Supplement S4 checks)

| Name | Definition |
|---|---|
| `epsilon` | $2/(Dk_Lw_0)$: longitudinal fields and single-electron pattern |
| `xi` | $[\beta/(1+\beta)]\,N_c\lambda_L/z_R$, $z_R=\tfrac12k_Lw_0^2$: path through the focus, in Rayleigh lengths |
| `fresnel_phase_error` | $k_{N_{\max},\rm ref}\,(\rho_{\max}+R)^4/(8Z^3)$: next term of the Fresnel phase over the screen |
| `paraxial_angle` | `theta_max` in rad, flagged above 0.1 |

**Physical block** (for the physical columns of the configuration table)

| Name | Definition |
|---|---|
| `lambda_L_nm`, `lambda_1_nm`, `photon_energy_eV` | $\lambda_L$, $\lambda_1$, $\hbar\omega_{1,\rm ref}$ |
| `w0_um`, `R_um` | waist and sheet radius |
| `sigma_field_fs`, `fwhm_intensity_fs` | $N_cT_L$ and $\sqrt{2\ln2}\,N_cT_L$ |
| `t_int_fs` | $0.59\,N_cT_L$: intensity FWHM seen by an electron, in lab time (head-on) |
| `kinetic_energy_MeV` | $(\gamma-1)mc^2$ |

## 5. Checks (printed, and listed in the run-log section)

| Check | Condition | Action |
|---|---|---|
| Geometry | $\hat{\mathbf n}_L\cdot\hat z<-1+10^{-9}$, $\mathbf p\parallel\hat z$, screen normal $+\hat z$ | otherwise: warning, and the scaled block is not written |
| Envelope | `flat_duration` = 0 | otherwise: warning, `N_c` = `None` |
| Frequency | `k1_over_kL` against the log's `fundamental / laser omega`, to $10^{-6}$ | warning (catches an $A_{\rm ref}$ or geometry mismatch) |
| Aliasing | `aliasing_margin` > 1.25 for every N | warning naming N |
| Containment | `theta_max_over_thetastar` ≥ 5 | note (the screen may clip, as `lg_decomposition.py`'s edge check measures) |
| Paraxial and Fresnel | `paraxial_angle` < 0.1 rad; `fresnel_phase_error` < 0.1 rad | note |
| Screen group | optional `--expect L` or `--expect H`: $Z/Z_*$ within 1% of 2.81 or 2.85, $\theta_{\max}/\theta_*$ of 5.35 or 8.75 | warning: catches a config whose screen was not set for the paper |

## 6. Interface

**Module:** `py_scripts/utils/scaled_parameters.py`. It is a leaf module, like `w0_axes_utils.py`: it must not import
`plot_field` or `plot_observables`, so any script can import it.

```python
sp = scaled_parameters(run_dir)          # dataclass: the fields of §4, per-harmonic rows in sp.harmonic_rows
sp.warnings                              # list of strings from §5
to_scaled_angles(x, y, sp)               # screen coordinates (lambda) -> (theta_x/theta_*, theta_y/theta_*)
dominant_component(sp)                   # 'plus' or 'minus'
```

**Command line:**

```
python3 py_scripts/utils/scaled_parameters.py [run_folder ...] [--write-log] [--expect L|H]
                                              [--table md|csv] [--output FILE]
```

- **One run:** prints all blocks.
- **`--write-log`:** writes the section `Scaled parameters (scaled_parameters.py)` into `run_log.txt`, replacing
  it on a re-run.
- **Several runs:** one row per run and harmonic, in Markdown or CSV. This is the source for the S4 table; LaTeX can
  come with the figure layer.

**Shared helper:** move `_write_run_log_section` and `_find_run_log_sections` from `plot_observables.py` to
`utils/run_output_utils.py`, and import them in `plot_observables.py` and `plot_point_spectrum_observables.py` (the
latter imports them from `plot_observables.py` today). Behaviour stays byte-identical: check by diffing `run_log.txt`
before and after on one run.

## 7. Verification (acceptance)

The module must reproduce numbers already documented elsewhere:

| Run | Quantity | Expected | Source |
|---|---|---|---|
| Report run A (`20261029_132340`; lattice `20261031_093730`) | `Z_over_Zstar`, `theta_max_over_thetastar` | 2.808, 5.342 | recomputed for this spec; the manuscript rounds to 2.81, 5.35 |
| Report run B (`…143739` / `…094539`) | same | 2.850, 8.773 | same; the manuscript quotes 2.85, 8.75 |
| Lens run, $a_0=0.1$ at rest (`20261030_111559`) | `Z_over_Zstar` | 10.24 | manuscript lens table (10.2) |
| Any run with $p_z=13.483\,mc$, $a_0=0.1$ | `D`, `k1_over_kL` | 27.00, 725.5 (log: 725.536) | `scaling_laws.md`, run logs |
| `~/output/20261031_130825` (laptop) | `Z_over_Zstar`, `theta_max_over_thetastar` | 1.404, 3.561 | computed for this spec |
| Report runs C, D (LG22, $l=2$, $\zeta_2=+i$) | `sigma`, `m`, `L_N` (N = 1, 2), `component` | −1, −2, −2 and −5, `minus` | report Sec. 3.4 |
| `reference_config` (if run) | `aligned`, `component` | `False`, `plus` | config: $l=2$, $\zeta_2=-i$ |

Also check one rectangular-screen run and one older run with missing log fields (warnings, no crash). The server runs
are not on the laptops: check those there, or against the numbers above.

## 8. Checkpoint

Implement §4–§6 and run §7. Then stop and show:
- the printed blocks for runs A and B;
- the multi-run table for the report runs;
- the `run_log.txt` diff of the helper move.

Then update `CLAUDE.md` (Plotting section: the module and its role as the single source of scaled units) and
`PLOTTING_NOTES.md`.

## 9. Open points

- **The $A_{\rm ref}$ decision** (after the spectrum runs). If it becomes per-mode, the solver needs a config key, and
  this module reads the value from the log.
- **Whether `to_scaled_angles` should also remove the carrier phase** $e^{-ik_N(R-Z)}$, or leave that to the figure
  layer. I suggest the figure layer, because it is a display choice.
- **`xi` uses $N_c$, while `scaling_laws.md` used $n_{\rm eff}=N_c\sqrt{\pi/2}$.** Pick one before the S4 table and
  say which in the supplement.
