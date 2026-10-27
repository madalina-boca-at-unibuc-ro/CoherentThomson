# LG decomposition of the screen field: implementation notes

Stage 2 of `theory/LG_mode_decomposition.md`. Task spec: `theory/lg_screen_decomposition_simple.md`. Code:
`py_scripts/radiation/lg_decomposition.py` (Python post-processing only; the C++ export is unchanged).

```
python3 py_scripts/radiation/lg_decomposition.py [--incident] [--harmonic N] [--component plus minus Ex Ey]
    [--l L] [--p-min P] [--p-max P] [--w0 20,25,30 [--w0-unit lambda|w0|a.u.]] [--n-w0 K] [run_folder]
python3 py_scripts/radiation/lg_decomposition.py --self-test
```

Defaults: components `plus minus`, `p = 0..3`, every frequency, and a 41-point log scan of `w0` over
`[1/2, 2] x` a second-moment estimate (below). Takes ~3 s on a 256 x 64 circular screen with two harmonics.

## What it computes

For each frequency, component and trial `w0`, weighted complex least squares on the screen grid,

$$\min_c \sum_i A_i\Big|E_i - \sum_{p,l} c_{p,l}\,u_{p,l}(r_i,\phi_i,s;w_0)\Big|^2 ,$$

with `A_i` the trapezoidal area weights of `plot_observables.get_screen_area_weights_au`. The rows are scaled
by `sqrt(A_i)` and solved with `numpy.linalg.lstsq` (SVD); the Gram matrix is never formed or inverted, and
orthogonality on the finite screen is not assumed. Reported:

- `c_{p,l}` (complex),
- `eps = sqrt(sum A|E - E_rec|^2 / sum A|E|^2)`, the spec's error,
- `captured = sum|c|^2 / sum A|E|^2`. With modes orthonormal on the screen this equals `1 - eps^2`. A large
  departure means they aren't: the screen clips the modes (see the incident test below: `captured = 2` on a
  screen of radius `2 w0` holding half of an `LG_2^2`) or the grid under-resolves them,
- `cond`, the weighted design matrix's condition number; rank deficiency or `cond > 1e8` is flagged in `status`.

**Reported `w0`** (`_choose_w0`, since 2026-10-27): among the interpretable trials (`|captured − (1 − eps²)| <
0.05`), the one where a single mode holds the largest share of `Σ|c|²`. For a single `LG_p^l` beam that is its own
waist, whatever `p` and `l`: at any other basis waist its power spreads over several `p`. Flagged `AT SCAN EDGE`
if it is the first or last value. The smallest-`eps` `w0` is printed on the next line as a basis label only.
It was the reported value before, but `eps(w0)` is flat (its residual is the other-`l` speckle, which no `w0`
reaches) and its minimum is set by noise: on run `20261027_194300` (at rest) it gave 17.6 λ for N = 1, where the
dominant-mode criterion gives 24.9 λ with 99.98% of `Σ|c|²` in `p = 0`, and 12.75 λ for N = 2, against 18.0 λ
(99.87% in `p = 0`; expected `w0/√2` = 17.7 λ).

## Basis and normalization

Unit-norm LG modes, `exp(-i omega t)` convention, waist at `s = 0`:

$$u_{p,l} = \sqrt{\frac{2\,p!}{\pi (p+|l|)!}}\frac{1}{w(s)}\Big(\frac{\sqrt2 r}{w}\Big)^{|l|}L_p^{|l|}\Big(\frac{2r^2}{w^2}\Big)e^{-r^2/w^2}\,e^{il\phi}\,e^{ikr^2/2R_c}\,e^{-i(2p+|l|+1)\psi},$$

`z_R = k w0^2/2`, `w = w0 sqrt(1+(s/z_R)^2)`, `1/R_c = s/(s^2+z_R^2)` (0 at the waist), `psi = atan2(s, z_R)`.
`L_p^a` by the three-term recurrence; the prefactor via `lgamma`. The laser's own C++ LG (custom `Npn`
normalization) is not reused.

- **The carrier `e^{iks}` is left out.** `s` is the same at every screen point, so it only adds `k s` to the phase
  of every coefficient. Relative phases between coefficients are meaningful; absolute ones are not.
- **Normalization does not affect `eps` or the reconstruction.** Least squares absorbs any per-column scale into
  the coefficient. It does matter for comparing `|c|^2` between different `p` and between `w0` trials, and for
  `captured`.
- `s` = the screen distance (`<type>_detector_distance`), positive: the emitted beam travels from the origin
  toward the screen in either direction. `s = 0` for `--incident`.
- `k = (omega/omega_1) * fundamental_frequency / c`, with `fundamental_frequency` read from `run_log.txt`.

## Components and azimuthal index (the two convention traps)

**Frame.** Everything is in the lab frame: lab `(x, y)` of each point (`screen_lab_xy`, the local grid rotated by
`get_detector_local_rotation`; on a −Oz screen local `y` = −lab `y`), `phi = atan2(y_lab, x_lab)`, and lab field
components (`E_i = c F^{i0}`, LR+SR+BR). So `l` and the spin are about **lab +z** whichever way the screen faces,
and `e^{il phi}` is an `L_z = l` eigenfunction regardless of the beam's propagation direction. `--incident` is
the exception: `incident_field.dat`'s grid and components are in the laser's own frame, so "lab" means the
laser frame there, and `s = 0`.

**Components.** A Cartesian component mixes both helicities, which have different `l`, so the default targets are
the helicity components

$$E_\pm = (E_x \mp iE_y)/\sqrt2 ,\qquad \mathbf E = E_+\hat e_+ + E_-\hat e_-,\ \hat e_\pm = (\hat x\pm i\hat y)/\sqrt2,$$

spin `±1` about lab +z (`S_z ∝ Im(E_x^* E_y) = |E_+|^2 − |E_-|^2`). Each is fitted separately with
`l = J − s`. `Ex`/`Ey` are available, and are fitted with both `l = J − 1` and `l = J + 1`.

**J.** `J = N(m + σ)`, with `m = laser_lg_l` and `σ = 2 Im(ζ₁* ζ₂)` (normalized) **times `sign(n_laser,z)`**, as in
`plot_observables.py`. The laser's `l` and spin are defined about its own direction, which is −Oz in the baseline
setup. Without the flip, the fitted `l` has the wrong sign and `eps ≈ 1`. For `--incident`, `J = m + σ` with no
flip. An elliptical `σ` (not 0 or ±1) makes `J` non-integer: the script stops and asks for `--l`. A laser not
along ±Oz only prints a warning.

**Azimuthal pre-check** (circular screens). Each ring is FFT'd in `φ` over its `N_phi − 1` unique columns
(`φ = 2π` duplicates `φ = 0`), with rings weighted by their area. This gives the power fraction per `l` for
`E_±` independently of the radial fit: it confirms `l` and shows how much power the fixed-`l` fit cannot reach.
On a −Oz screen the FFT index is multiplied by `det(R_⊥) = −1` to give lab `l`.

## Default w0 scan

Built around the waist of a `p = 0` mode with the same second moment of `|E|^2`, `w(s)^2 = 2<r^2>/(|l|+1)`
(`<r^2>` of `LG_p^l` is `(2p+|l|+1) w^2/2`). That is inverted on the far branch exactly as in Stage 1
(`w0^2 = 2B^2/(A + sqrt(A^2 − 4B^2))`, `B = 2s/k`); if the screen width is below the minimum possible
(`A < 2B`), it uses `w0 = sqrt(B)`. At `s = 0`, `w0 = w(s)`. Explicit values: `--w0 20,25,30 --w0-unit lambda`.

## Other diagnostics

- **Spherical vs paraxial phase.** The exported field has the exact `k(sqrt(s^2+r^2) − s)`; the modes have
  `k r^2/2R_c`. Their difference at the radius holding 99% of the component's power, `k r^4/(8 s^3)`, is printed.
  It was ~1e-7 rad on the runs below. If it ever reaches ~0.1 rad, the paraxial basis is no longer adequate.
- **Resolution.** The summary prints the grid step against the beam width `w(s)` at the reported `w0` (circular:
  `r²` step / `w²`; rectangular: `max(dx, dy)/w`) and flags it **UNDER-RESOLVED** above 0.25 / 0.5. This matters
  most for `--incident`: `incident_field.dat` uses the detector's grid at the laser waist, so a screen sized for
  the far field cannot sample the incident beam. Run `20261027_184830` (`R_max = 40 w0`, 256 equal-area rings) has
  its first ring at `2.5 w0`; the whole incident Gaussian sat on the centre point, and the fit came out rank 1 with
  equal coefficients and a 0.1 λ "waist". Fixed in the C++ on 2026-10-27: the incident screen is now at the
  detector distance with a beam-sized grid (`incident_config.cfg`, Verification 5).
- **Quadrature.** The trapezoid in `u = r^2` is only `O(h^2)` for `l = 0, ±1` (`|u|^2` has a nonzero `u`-slope at
  `r = 0`; relative error `~(2h/w^2)^2/12` with `h` the `r^2` step), and spectrally accurate for `|l| ≥ 2`. A
  coarse radial grid biases `captured` slightly for low `|l|`; target and basis are integrated with the same
  weights, so `eps` and `c` are much less affected.

## Outputs

- `<run>/lg_decomposition.dat` (`lg_decomposition_incident.dat` for `--incident`): one row per (frequency,
  component, w0) with `eps, captured, cond, rank`, `c_p<p>_l<l>_re/_im` columns, and, for single-`l` components,
  `w0_series_lambdaL, z0_series_lambdaL` (Series waist above).
- `run_log.txt` section `LG decomposition (lg_decomposition.py)` (`-- incident beam` suffix for `--incident`),
  replaced on re-run. It holds the reported `w0` per (frequency, component) with all coefficients, the series waist
  (with both `t`), the helicity power fraction, the azimuthal pre-check and the spherical-phase diagnostic.
- `png_folder/radiation/lg_decomposition/{emitted,incident}/lg_decomposition_omega<i>.png`: `eps(w0)` per
  component, and `|c|^2/sum|c|^2` at the reported `w0`.
- `.../lg_coefficients_omega<i>_<component>.png`: the whole coefficient set of one component, normalized by
  `||E||^2 = sum A|E|^2`. Components whose helicity power fraction is below `1e-4` are skipped (their fit is
  noise); `Ex`/`Ey` are always plotted. Three panels:
  1. `|c_{p,l}|^2/||E||^2` vs `w0` for every mode, with `eps^2` and `captured` overlaid (gray line: reported `w0`).
  2. Argand diagram of `c/||E||` at the reported `w0`, rotated so the dominant mode is real (gray unit circle: all of
     the field in one mode).
  3. `arg(c/c_dominant)/π` vs `w0`, masked where `|c|^2/||E||^2 < 1e-6`.

### Reading the coefficient plot

- **The coefficients don't depend on `s`.** Each mode carries its own Gouy phase `e^{-i(2p+|l|+1)ψ}`, so `c_{p,l}`
  are the beam's expansion coefficients at the waist plane, wherever the screen is. Their relative phases are
  physical; only the common `e^{iks}` was removed.
- **Small `w0` = clipping.** `w(s) ∝ 1/w0` in the far field, so below some `w0` the modes are wider than the
  screen and stop being orthonormal on it. `captured ≫ 1`, and `|c|^2` blow up with large cancelling values (left
  part of the scan). Only the region where `captured ≈ 1 − eps^2` is interpretable.
- **A displaced or resized Gaussian gives a geometric series.** A Gaussian with a different complex beam
  parameter, expanded in this basis, gives `c_p ∝ ξ^p` for some complex `ξ`. `ξ` is real for a pure waist-size
  mismatch and close to imaginary (`≈ −iδ/2z_R`) for a pure waist shift `δ`. So `|c_2/c_1| = |c_1/c_0|` and
  `arg c_2 − arg c_1 = arg c_1 − arg c_0`.
  - N = 1, run `20261025_002302`, at the best `w0`: `|c_1/c_0| ≈ 0.68` with `arg(c_1/c_0) ≈ 0.7π`. A geometric
    series would then predict `|c_2|^2 ≈ 0.14`, but the fit gives 0.015 (`|c_2/c_1| = 0.22`). The phase steps are
    similar (`0.68π`, then `0.80π`), but the magnitudes are not. This argues against "just a displaced
    and/or resized Gaussian" as the explanation of the `p = 1` content. Tentative: the series test assumes `eps`
    is small, and here it is 5%. Resolved by the `a0` test below: at `a0 = 0.1` it isn't a single Gaussian (the
    ponderomotive radial phase); at `a0 = 0.01` it is, with the laser's waist.

## Series waist

For a target that is one `LG_0^l` beam with its own waist `w` and waist position `z0` (any complex beam parameter
`q`), the coefficients in the waist-at-origin basis of waist `w_b` follow from the Laguerre generating function
`Σ_p L_p^a(x) t^p = (1−t)^{−a−1} e^{−xt/(1−t)}`. On the plane `s = 0` the target is `r^|l| e^{−α r²} e^{ilφ}`
(`α = −ik/2q`), so on the unit-norm modes

$$\frac{c_{p+1}}{c_p} = t\,\sqrt{\frac{p+1+|l|}{p+1}},\qquad t = \frac{\alpha w_b^2 - 1}{\alpha w_b^2 + 1},$$

independent of `s` (the basis carries the Gouy phases). `series_waist` takes `t` from `c_1/c_0` and inverts it:
`q(s=0) = −i z_Rb (1−t)/(1+t)` with `z_Rb = k w_b²/2`, giving `z0 = −Re q` (the waist position along the
propagation direction, from the origin toward the screen) and `w0 = sqrt(2(−Im q)/k)`. Real `t` means the waist
is at the origin: `w = w_b sqrt((1−t)/(1+t))`.

Two checks that the single-beam assumption holds, both reported:
- `t` from `c_2/c_1` (the summary's second `t`) must equal `t` from `c_1/c_0`. At the reported `w0` the basis
  is usually the beam's own waist, `c_1, c_2 → 0` and the ratio is noise, so for `|t| < 0.05` the summary skips
  it and the next check is the one to use;
- `w0_series`, `z0_series` (columns of `lg_decomposition.dat`, one per trial `w0`) must not depend on the basis
  waist, wherever `captured ≈ 1 − eps²`.

This replaces the smallest-`eps` `w0` as the waist estimate: that one is only a basis label, and with several `p`
free and a `w0`-independent residual (other-`l` content), `eps(w0)` is flat. Computed only for components fitted
with a single `l` and with `p = 0, 1` in the basis.

## Verification

1. `--self-test` (synthetic, no run needed), all passing:
   - `∫|u|^2 dA = 1` to 1e-4 for `(p,l) = (0,0), (1,0), (2,3), (0,−2), (4,1)` at `s = 0` and `s = 3 z_R`.
     Off-diagonal overlaps `<0,2|1,2>`, `<0,2|0,−2>`, `<1,1|3,1>` are below 1e-4. This needs a fine radial grid
     (6000 rings) because of the `l = 0, 1` quadrature note above; with 600 rings, `(1,0)` is off by 1e-2.
   - Known coefficients (`E_+`: `l = 2`, `p = 0, 1, 2`; `E_-`: `l = 4`) are recovered to 1e-15, with `eps ~ 1e-15`,
     on a +Oz **and a −Oz** screen. The field is synthesized from an independent lab-coordinate mapping, so the
     test fails if `screen_lab_xy` or the helicity convention is wrong. Fitting `Ex` with both `l` also gives
     `eps ~ 1e-15`.
   - Series waist: a single `LG_0^l` beam outside the basis (`l = 0`: `w = 1.3 w_b`, `z0 = 0` and `w = 0.8 w_b`,
     `z0 = 0.3 z_R`; `l = 2`: `w = 1.3 w_b`, `z0 = −0.4 z_R`), fitted with `p = 0..8` on a wide fine screen:
     `w0`, `z0` recovered to ~1e-8 and the two `t` agree, which fixes the signs of `z0` and the `|l|` factor.
2. **Incident field, run `20260929_002928`** (LG `p = 2, l = 2`, `w0 = 25 λ`, `ζ = (1, i)`, screen radius
   `50 λ`): all power is in `E_+` with `l = +2` (ring FFT: 1.0). At `w0 = 25 λ`, `c_{2,2}` is the only nonzero
   coefficient (the others are ~1e-9 of it), `eps = 1.1e-7` (the 6-digit `.dat` precision), and `eps` rises
   to 6e-5 at `w0 = 24, 26 λ`. So **the C++ `LaguerreGaussLaser` transverse mode at its waist equals the
   standard LG mode** up to normalization (part of CLAUDE.md's open "not physics-reviewed" issue; the Gouy and
   curvature terms are not tested here, since `s = 0`). `captured = 2.0`: the `50 λ = 2 w0` screen holds only
   half of this mode's norm.
3. **Emitted field, run `20261025_002302`** (Gaussian laser along −Oz, circular `ζ = (1, i)`, so `σ_lab = −1`,
   `γ ≈ 13.5`, 8192 electrons, forward screen at `1.6e7 λ_L`, `p = 0..3`):

   | N | J | E₋ power | ring FFT (E₋) | best w0 [λ_L] | eps | captured | |c|² share p=0,1,2,3 |
   |---|---|----------|---------------|---------------|-----|----------|----------------------|
   | 1 | −1 | 1 − 6e-13 | l=0: 0.997 | 15.8 | 0.052 | 0.992 | 0.67, 0.31, 0.015, 0.008 |
   | 2 | −2 | 1 − 7e-6 | l=−1: 0.882, l=0: 0.090 | 9.8 | 0.35 | 0.89 | 0.31, 0.45, 0.21, 0.03 |

   The predicted helicity and the sign of `l` are confirmed: lab `l = J + 1` in `E_-`.
4. **`a0` test, runs `20261027_181823` (`config_expansion.cfg`, `a0 = 0.1`) and `20261027_182248`
   (`config_expansion_a0_001.cfg`, `a0 = 0.01`)**, commit `e4c3589`, otherwise identical (seed 532, 8192
   electrons, flat disk of radius `3 w0`, forward screen at `4e6 λ_L` = `Z/z_R ≈ 2.8` for `w0 = 25 λ`,
   `R_max = 10 w0`). N = 1, `E_-`, `l = 0`:

   | | `a0 = 0.1` | `a0 = 0.01` |
   |---|---|---|
   | Stage 1 `w0_far` | 17.7 λ | 23.5 λ |
   | smallest-`eps` basis `w0` | 17.1 λ | 19.1 λ |
   | `|c|²` share `p = 0, 1, 2, 3` | 0.70, 0.29, 0.004, 0.008 | 0.925, 0.07, 0.004, 0.001 |
   | `t` from `c_1/c_0`; from `c_2/c_1` | `0.64 e^{0.66iπ}`; `0.12 e^{0.59iπ}` | `0.275 e^{−iπ}`; `0.238 e^{−0.99iπ}` |
   | series `w0` over the interpretable scan | 12.4 → 22.9 λ (drifts) | 25.2–25.4 λ |
   | series `z0` | ~9e5 λ_L | −2e3 to −9e3 λ_L (`< 0.006 z_R`) |
   | `eps²`; power outside `l = 0` | 0.0094; 0.0093 | 0.0090; 0.0090 |

   - **At `a0 = 0.01`, N = 1 is one Gaussian beam with the laser's waist (25 λ) at the origin**: `t` is real, the
     two `t` agree, and the series waist is stable across basis waists.
   - **At `a0 = 0.1` it isn't a single beam** (the two `t` disagree, the series waist drifts). The extra `p = 1`
     and the narrower apparent waist come from ponderomotive dressing: the frequencies are the on-axis dressed
     harmonics, so an electron at `r` is off resonance by `∝ a0² (1 − e^{−2r²/w0²})` and its emission carries a
     radial phase of order 1 rad between `r = 0` and `w0` (~65 cycles crossed). That acts as a Gaussian-shaped
     lens, not a pure curvature, hence no single-beam fit. It falls 100x at `a0 = 0.01`. An intensity profile
     entering squared (waist `w0/√2` = 17.7 λ, matching `a0 = 0.1`'s Stage 1 by coincidence) would not depend
     on `a0`, so it is ruled out.
   - **The residual `eps` is the other-`l` content in both runs** (`eps²` equals the power outside `l = 0`), so it
     doesn't depend on `a0` or `w0`. Its size (~1%) matches the speckle of random electron positions: the
     incoherent part carries `≈ Q/((2π)² n)` of the screen power (`Q` the screen's area in transverse `q`, `n` the
     electron surface density), ≈ 1.4% here. It also explains Stage 1's 23.5 λ at `a0 = 0.01`: a clean 25 λ
     Gaussian (74.2 λ spot) plus ~1% uniform background gives `w_eff` ≈ 77.6–79.4 λ, measured 77.9 λ.
   - **N = 2 at `a0 = 0.01` isn't usable yet**: it is ~`a0²` weaker than N = 1 and `trajectory_NT` is 100 (weak
     harmonics need finer, CLAUDE.md); its `l = −1` share drops to 0.71, `edge_frac` rises to 12% and the `plus`
     fit runs off the scan. Rerun with a finer `trajectory_NT` before interpreting it.
5. **Incident field at the detector distance** (2026-10-27). `incident_field.dat` is the plane `s = Z` downstream
   of the waist (`incident_config.cfg`; `--incident` uses `s = Z`). On `config_expansion_rest.cfg` (Gaussian,
   `w0 = 25 λ`, `Z = 2e4 λ = 10.2 z_R`, 128 x 33 screen of radius 952 λ): `eps = 1.0e-7`, all power in
   `E_+` with `l = 0`, series waist 25.000 λ at `z0 = 0.0009 λ`. This tests Gouy and curvature, which the waist
   check couldn't, and caught a curvature-sign bug in `LaguerreGaussLaser` (with the old sign `eps = 0.94`;
   `theory/dev_notes.md`).

## Open points

- **N = 1's `p = 1` content (31% at `a0 = 0.1`) is ponderomotive dressing**, not a displaced waist or the
  incoherent background (Verification 4). Not yet checked quantitatively: the predicted radial phase profile
  against the fitted field, and the scaling with `a0²` and with the pulse length.
- **Speckle:** ~1% of the screen power is incoherent, at every `l`. Halving it needs 2x the electrons (it scales as
  `1/n`); not tested.
- **N = 2 has 9% of `E_-` power at `l = 0`** instead of `−1`. That alone puts `eps` near 0.3, and it would lower
  `|L_z|` per photon by roughly that fraction. It may be related to the ~5% `Λ_zz/P_z` deficit at N = 2
  (CLAUDE.md open issue); not checked.
- `|c|^2` are not additive modal powers on a finite screen unless `captured ≈ 1 − eps^2`.
- `w0` is a basis label, scanned not optimized; with several `p` values free, `eps(w0)` is shallow, so the
  smallest-`eps` `w0` trades off against the `p` distribution. The reported `w0` (dominant mode) and the series
  waist are the waist estimates.
