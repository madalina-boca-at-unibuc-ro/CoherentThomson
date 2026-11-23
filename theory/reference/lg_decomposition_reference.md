# `lg_decomposition.py`: reference

What `py_scripts/radiation/lg_decomposition.py` computes, its command-line parameters and their defaults, and what
it writes. Design history, verification runs and open points are in `theory/lg_decomposition_implementation.md`;
the physics plan is `theory/LG_mode_decomposition.md` (this script is its Stage 2).

The script post-processes one finished run. It changes nothing in the C++ output except appending (or replacing) its
own section of `run_log.txt`.

```
python3 py_scripts/radiation/lg_decomposition.py [run_folder] [--incident] [--harmonic N]
    [--component {plus,minus,Ex,Ey} ...] [--l L] [--p-min P] [--p-max P]
    [--w0 V1,V2,... [--w0-unit U]] [--n-w0 K]
    [--scan-z0] [--z0 V1,V2,... [--z0-unit U]] [--n-z0 K]
python3 py_scripts/radiation/lg_decomposition.py --self-test
```

## 1. Parameters

| Parameter | Default if omitted | Meaning |
|---|---|---|
| `run_folder` | the latest run in the output folder | Run directory holding `radiation_field.dat` (or `incident_field.dat`), `config.cfg` and `run_log.txt`. A folder, not a `.dat` path. |
| `--incident` | off: the emitted field `radiation_field.dat` | Analyse the incident laser field `incident_field.dat` instead (§6). |
| `--harmonic N` | every frequency in the file | Only the frequency with `round(ω/ω₁) = N`. |
| `--component C [C ...]` | `plus minus` | Field components to decompose, any of `plus`, `minus`, `Ex`, `Ey` (§3.1). |
| `--l L` | from angular-momentum conservation, `l = J − s` (§3.2) | Force one azimuthal index for every component. |
| `--p-min P` | `0` | Smallest radial index in the basis. |
| `--p-max P` | `3` | Largest radial index. Needs `0 ≤ p-min ≤ p-max`. |
| `--w0 V1,V2,...` | 41 log-spaced values over `[w_est/2, 2 w_est]` (§3.4) | Explicit trial basis waists, comma-separated, positive. |
| `--w0-unit U` | `lambda` | Unit of `--w0`: `lambda` (laser wavelength λ_L), `w0` (the laser's `laser_lg_w0`, LG lasers only), `a.u.`, or any other unit the config parser knows. |
| `--n-w0 K` | `41` | Number of points in the default `w0` scan. Ignored with `--w0`. |
| `--scan-z0` | off: basis waist fixed at `z0 = 0` | Also scan the basis waist position (§3.5). |
| `--z0 V1,V2,...` | none | Explicit trial waist positions, any sign. Turns the `z0` scan on by itself. |
| `--z0-unit U` | `lambda` | Unit of `--z0` (same choices as `--w0-unit`). |
| `--n-z0 K` | `41` | Number of points in the default `z0` scan. Rounded up to an odd number so that `z0 = 0` is included. Ignored with `--z0`. |
| `--self-test` | off | Run the synthetic checks (§7) and exit; no run folder needed. Exit code 1 if any check fails. |

### Fixed settings (module constants, not on the command line)

| Constant | Value | Role |
|---|---|---|
| `_W0_SCAN_FACTOR` | 2 | Default `w0` scan spans `[w_est/2, 2 w_est]`. |
| `_Z0_SCAN_ZR` | 3 | Default `z0` scan spans `±3 z_R` of the scan's central waist. |
| `_COND_WARN` | 1e8 | Condition number above which a fit's status reads `ill-conditioned`. |
| `_EDGE_BAND`, `_MAX_EDGE_FRACTION` | 0.9, 1e-2 | Containment check: power beyond 0.9 × the screen's inscribed radius, flagged above 1% (§4.3). |
| `_MAX_STEP_OVER_WIDTH` | circular 0.25, rectangular 0.5 | Resolution check thresholds (§4.4). |
| `_MIN_HELICITY_FRACTION_PLOT` | 1e-4 | Helicity components carrying less of the power get no coefficient figure. |
| `_MIN_T_FOR_CHECK` | 0.05 | Below this `|t|` the series-waist consistency check is skipped (§4.2). |
| `_SPIN_TOL` | 1e-9 | Laser spin must be within this of 0 or ±1, else `J` is not an integer and `--l` is required. |
| `_PREDICTION_CONTAINED` | 0.999 | Radius holding 99.9% of the predicted power, used for the "outside" share (§4.6). |
| `_SOURCE_RMAX_W0` | 8 | Source radius cap (laser waists) of the prediction when the run log gives no beam radius. |

## 2. Inputs and requirements

- **Field:** `radiation_field.dat` (default) or `incident_field.dat`. The total field is used, `E = c·F^{i0}`
  summed over the long-range, short-range and boundary terms.
- **Config:** the run's `config.cfg` (for `--incident`, `incident_config.cfg` when present). Read: detector type,
  distance and grid, laser polarization (`ζ₁, ζ₂`), `laser_lg_l`, laser direction, and for the prediction
  `laser_type`, `laser_lg_p`, `laser_lg_w0`, `average_px/py/pz`.
- **Run log:** `fundamental_frequency` (sets `k`) and the beam `cylinder radius` (prediction only).
- **Screen:** flat `rectangular` or `circular` only; spherical screens are rejected. For the emitted field the
  screen must face ±Oz (else an error). A laser not along ±Oz only prints a warning, because `J = N(m + σ)` about lab
  `z` then doesn't hold.
- **Frequencies:** `N = round(ω/ω₁)`; a warning is printed if `ω/ω₁` is not an integer. `--incident` uses `N = 1`.

## 3. What is computed

Everything is done separately for each selected frequency and component.

### 3.1 Field components

The default targets are the **helicity components** in the lab frame,

$$E_\pm = \frac{E_x \mp iE_y}{\sqrt2},\qquad \mathbf E = E_+\hat e_+ + E_-\hat e_-,\quad \hat e_\pm = \frac{\hat x \pm i\hat y}{\sqrt2},$$

which carry spin `s = ±1` about lab `+z` (`plus` ↔ `+1`, `minus` ↔ `−1`). Each has a single orbital index, so it
can be fitted with one `l`. `Ex` and `Ey` mix both helicities and are fitted with both `l = J − 1` and `l = J + 1`.

Coordinates are lab-frame `(x, y)` of each screen point (a screen facing −Oz has local `y` = −lab `y`), with
`r = √(x²+y²)`, `φ = atan2(y, x)`, and `s` = the detector distance (positive). For `--incident` everything is in the
laser's own frame instead.

### 3.2 Azimuthal index

Angular-momentum conservation per photon about lab `+z`:

$$J = N\,(m + \sigma),\qquad l = J - s,$$

with `m = laser_lg_l` and `σ = 2 Im(ζ₁* ζ₂)` (the laser's helicity), **both multiplied by `sign(n_laser,z)`**
because they are defined about the laser's own direction (−Oz in the baseline setup). `σ` must be 0 or ±1;
otherwise give `--l`. For `--incident`: `J = m + σ` without the sign flip. The summary prints `J` and both
`m, σ`.

### 3.3 LG basis

Unit-norm modes with the waist at `s = 0`, `exp(−iωt)` convention:

$$u_{p,l} = \sqrt{\frac{2\,p!}{\pi (p+|l|)!}}\;\frac{1}{w}\Big(\frac{\sqrt2\, r}{w}\Big)^{|l|} L_p^{|l|}\!\Big(\frac{2r^2}{w^2}\Big)\, e^{-r^2/w^2}\, e^{il\phi}\, e^{ikr^2/2R_c}\, e^{-i(2p+|l|+1)\psi},$$

`z_R = k w0²/2`, `w = w0 √(1+(s/z_R)²)`, `1/R_c = s/(s²+z_R²)`, `ψ = atan2(s, z_R)`, `k = (ω/ω₁)·ω₁/c`. They satisfy
`∫|u|² dA = 1` over the infinite plane at any `s`. The carrier `e^{iks}` is left out, so every coefficient carries a
common phase `ks`: relative phases between modes are physical, absolute ones are not. Because each mode carries its
own Gouy phase, the coefficients describe the beam **at its waist plane**, wherever the screen is.

The basis is `{u_{p,l}: p_min ≤ p ≤ p_max}` for the component's `l` (two `l` values for `Ex`/`Ey`), all with
one common `w0` (and `z0`).

### 3.4 Trial waists `w0`

Default: a geometric scan of `n_w0` values over `[w_est/2, 2 w_est]`. `w_est` is the waist of a `p = 0` mode with
the same second moment as `|E|²` on the screen:

$$w(s)^2 = \frac{2\langle r^2\rangle}{|l|+1},\qquad w_0^2 = \frac{2B^2}{w(s)^2 + \sqrt{w(s)^4 - 4B^2}},\quad B = \frac{2s}{k}$$

(the far-field branch). If the field is narrower than any Gaussian can be at `s`, `w0 = √B`. At `s = 0`,
`w0 = w(s)`. With `--w0`, the values are used as given.

### 3.5 Trial waist positions `z0` (`--scan-z0` / `--z0`)

The basis waist moves to `z0` along the propagation: the modes are evaluated at `s − z0`. Default grid:
`n_z0` (odd) values, linear over `±3 z_R`, with `z_R = k w_c²/2` and `w_c` the geometric mean of the first and last
trial `w0`. Every `(w0, z0)` pair is fitted, so the cost is `n_w0 × n_z0` fits per component and frequency (~50 s for
two components at two frequencies on a 256 × 64 screen with the defaults). Without the scan, `z0 = 0`.

### 3.6 Per trial `(w0, z0)`: projections and least squares

With `A_i` the screen's trapezoidal area weights and `‖E‖² = Σ A_i |E_i|²`:

**Projections** (the primary result):

$$a_{p,l} = \sum_i A_i\, u_{p,l}^*(r_i,\phi_i)\, E_i .$$

Since the modes are orthonormal on the infinite plane, if `E` vanishes outside the screen, `a_{p,l}` is **exactly**
the field's expansion coefficient, however much of the mode the finite screen clips. So `|a_{p,l}|²/‖E‖²` is that
mode's share of the field's power. Derived:

- **mode overlap** `max_{p,l} |a|²/‖E‖²`: the share of the whole field in the best single mode (the mode-matching
  efficiency; 1 only for a pure `LG_p^l` with this `w0`, `z0`);
- `Σ|a|²/‖E‖² ≤ 1`: the share of the field the fitted modes can reach.

**Weighted least squares** (reconstruction diagnostics):

$$\min_c \sum_i A_i\Big|E_i - \sum_{p,l} c_{p,l}\,u_{p,l}(r_i,\phi_i)\Big|^2 ,$$

solved by SVD on `√A`-scaled rows. Reported:

- `eps = √(Σ A|E − E_rec|² / ‖E‖²)`: relative residual of the best reconstruction **on the screen**;
- `captured = Σ|c|²/‖E‖²`: equals `1 − eps²` only when the modes are orthonormal on the screen;
- `cond`: condition number of the weighted design matrix (1 = orthonormal on the screen); `status` says
  `rank deficient` or `ill-conditioned` (`cond > 1e8`) when relevant.

`c = G⁻¹a`, with `G` the modes' overlap matrix on the screen. At small trial `w0` the modes are wide on the screen
(`w(s) ∝ 1/w0` in the far field), get clipped, and `G ≠ 1`. Then `c` no longer gives power fractions; `a` still does.

### 3.7 Reported `(w0, z0)`

The trial with the **largest mode overlap** `max |a|²/‖E‖²`. For a single `LG_p^l` beam that is its own waist and
waist position, whatever `p` and `l`. For a field that isn't one mode, it is the best single-mode approximation, and
the overlap value says how far from a single mode the field is. It doesn't depend on `p_max` (each `a` is a separate
integral). `AT SCAN EDGE` is printed if the chosen `w0` (or `z0`) is the first or last value of its scan: widen the
scan. The trial with the smallest `eps` is printed too, but only as a basis label: with several `p` free, `eps(w0)` is
flat and its minimum is set by noise.

## 4. Diagnostics (per frequency and component, at the reported point)

### 4.1 Helicity power fraction

`‖E_±‖²/(‖E_+‖² + ‖E_-‖²)`. Shows which helicity carries the harmonic; a component below 1e-4 gets no coefficient
figure (its fit is noise).

### 4.2 Series waist

Assumes the component is a single `LG_0^l` beam with its own waist `w` and waist position `z0` (any complex beam
parameter `q`). Its expansion in the basis of waist `w_b` is geometric,

$$\frac{a_{p+1}}{a_p} = t\,\sqrt{\frac{p+1+|l|}{p+1}},\qquad t = \frac{\alpha w_b^2 - 1}{\alpha w_b^2 + 1},$$

so `t` from `a_1/a_0` gives `q(s=0) = −i z_Rb (1−t)/(1+t)` (`z_Rb = k w_b²/2`), hence `z0 = −Re q` (along the
propagation, from the origin) and `w0 = √(2(−Im q)/k)`. Checks: `t` from `a_2/a_1` must equal `t` (skipped when
`|t| < 0.05`, where both ratios are noise), and the `w0_series_lambdaL`, `z0_series_lambdaL` columns must not
depend on the trial `w0`. Only for single-`l` components with `p = 0, 1` in the basis. `NaN` when `|t| ≥ 1` (not a
single-beam series), e.g. for a field dominated by `p ≥ 1`.

### 4.3 Containment

Share of `‖E‖²` beyond 0.9 × the screen's inscribed radius (`R_max` for a circular screen, the shorter half-width for
a rectangular one). Above 1% the summary prints **NOT CONTAINED**: the field probably continues past the screen, so
the projections (and the reported point) miss that part. About 0.2% comes from the uniform ~1% incoherent speckle
alone; more than that means coherent field at the edge.

### 4.4 Resolution

Grid step against the beam width at the reported point, `w(s − z0)`: the `r²` step over `w²` for a circular grid
(flag above 0.25), `max(dx, dy)/w` for a rectangular one (flag above 0.5). Flagged as **UNDER-RESOLVED**: nothing
in that block is meaningful. The trapezoid in `r²` is `O(h²)` for `l = 0, ±1`, so at a typical `r²` step of
`0.06 w²` the projections carry ~1e-3 relative error (least squares cancels most of it).

### 4.5 Spherical vs paraxial phase

`k r⁴/(8 s³)` at the radius holding 99% of the component's power: the difference between the exported field's exact
spherical phase and the paraxial `k r²/2R` of the modes. ~0.1 rad or more means the paraxial basis is no longer
adequate.

### 4.6 Azimuthal spectrum (circular screens, `plus`/`minus`)

FFT of each ring over its `N_phi − 1` unique angles, rings weighted by their area: the power fraction per `l`,
independent of the radial fit. Prints the fraction at the expected `l` and the three largest. It confirms `l` and
shows how much power no fixed-`l` fit can reach.

### 4.7 Coherent-source prediction (LG laser, electrons at rest)

A parameter-free model of the emitted component: the laser mode at its waist to the `N`-th power over the beam
cylinder (`r < R_beam` from `run_log.txt`, or `8 w0` if absent), Fresnel-propagated to the screen, times each
electron's near-axis angular factor `(r/s)^{|Δl|} e^{iΔlφ}` (`Δl = l − Nm`). Reported:

- **fidelity** `|⟨E_pred|E⟩|²/(‖E_pred‖²‖E‖²)`: the share of the field explained with only an overall complex scale
  free;
- the measured power outside the radius that holds 99.9% of the prediction's power.

Skipped (with the reason) unless `laser_type = laguerre_gauss` and `average_px = average_py = average_pz = 0`.
Not modelled: ponderomotive dressing, the incoherent speckle, a moving beam. For harmonics `N ≥ 2` of a structured
laser this, not the LG fit, is the relevant test: those harmonics are generally not single LG modes.

## 5. Outputs

All in the run folder (`_incident` suffix / `incident` subfolder for `--incident`).

- **`lg_decomposition.dat`**: one row per (frequency, component, `w0`, `z0`), space-separated. Columns:
  `i_omega, omega, N, component, l, w0_au, w0_lambdaL, z0_au, z0_lambdaL, eps, captured, cond, rank, norm_sq,
  proj_captured, overlap`, then per mode `c_p<p>_l<l>_re/_im` (least squares) and `a_p<p>_l<l>_re/_im`
  (projections), and, for single-`l` components, `w0_series_lambdaL, z0_series_lambdaL`. Coefficients are in the
  field's units × length (`norm_sq = ‖E‖²`).
- **`run_log.txt`**, section `LG decomposition (lg_decomposition.py)` (`-- incident beam` for `--incident`),
  replaced on every re-run: the conventions, then per frequency and component the reported point with its mode
  overlap, `Σ|a|²`, `eps`, `captured`, `cond` and status, every §4 diagnostic, and every mode's `a` with
  `|a|²/‖E‖²` (and the least-squares `|c|²/Σ|c|²` for comparison).
- **Figures** in `png_folder/radiation/lg_decomposition/{emitted,incident}/`:
  - `lg_coefficients_omega<i>_<component>.png`: `|a|²/‖E‖²` vs `w0` per mode with `Σ|a|²/‖E‖²` and the
    least-squares `eps²`, and the least-squares `captured` on its own right-hand log axis (it grows far above 1 at
    small trial `w0`, where the clipped modes aren't orthonormal on the screen, and would otherwise flatten the
    other curves); the Argand diagram of `a/‖E‖` at the reported point, rotated so the
    dominant mode is real (unit circle = the whole field in one mode); `arg(a/a_dominant)/π` vs `w0`, masked where
    `|a|²/‖E‖² < 1e-6`;
  - `lg_w0_z0_map_omega<i>_<component>.png` (`z0` scan only): mode overlap and least-squares `eps` over the
    `(w0, z0)` grid, reported point marked.

## 6. `--incident`

Decomposes the incident laser field (`incident_field.dat`) instead of the emitted one: a check of the laser mode
itself and of the basis conventions. The grid and components are in the laser's own frame, `N = 1`, `J = m + σ`, and
`s` is the detector distance (the plane the C++ writes the incident field on, `incident_config.cfg`), or 0 for runs
before 2026-10-27, which wrote it at the waist. A pure `LG_p^m` laser must come out as the single mode
`(p, m)` in `E_σ`, with `eps ~ 1e-7` (the `.dat` precision) and mode overlap ≈ 1 at the laser's `w0`, `z0 = 0`.

## 7. `--self-test`

Synthetic checks, no run needed:

1. normalization and orthogonality of the modes on a fine grid, at the waist and at `3 z_R`;
2. recovery of known coefficients on a +Oz and a −Oz screen, through the same lab-frame and helicity code the run
   analysis uses;
3. series waist of a single `LG_0^l` beam with a different waist and waist position;
4. the `(w0, z0)` scan picks a displaced `LG_1^2`'s true pair with mode overlap 1;
4b. a field contained on a screen that clips the basis: projections equal those on a wide screen (least squares
   doesn't);
5. Fresnel propagation behind the coherent-source prediction against the analytic modes.

## 8. Examples

```
python3 py_scripts/radiation/lg_decomposition.py                       # latest run, both helicities, every N
python3 py_scripts/radiation/lg_decomposition.py ~/output/20261027_213655 --scan-z0
python3 py_scripts/radiation/lg_decomposition.py --harmonic 2 --component minus --p-max 9
python3 py_scripts/radiation/lg_decomposition.py --w0 20,25,30 --w0-unit lambda --z0 -500,0,500
python3 py_scripts/radiation/lg_decomposition.py --incident             # check the laser mode
```
