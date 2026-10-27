# LG mode decomposition of the emitted field

Long-term goal: write the emitted field on the detector as a superposition of Laguerre-Gauss modes and identify
their parameters (waist, waist position, `p`, `l`, weights). This file tracks the stages as they are implemented.
Earlier stage: radial-phase removal, `theory/post_processing.md` (`plot_processed.py radial_phase`).

All of it is Python post-processing of one run (one detector distance `Z`); the C++ export is unchanged.

## Stage 1: Gaussian waist from the Poynting-flux second moment

Spec: `theory/gaussian_waist_postprocessing.md`. Implemented in `py_scripts/radiation/plot_observables.py`
(`estimate_gaussian_waist`, `append_waist_to_run_log`). It runs on every non-`--incident` call and writes the
section "Gaussian waist estimate (plot_observables.py)" to `run_log.txt` (also printed to stdout). Not run on
`--incident`: that field is sampled at the laser's own waist.

### Width

Per frequency,

$$w_{\rm eff}^2(Z) = 2\,\frac{\sum_p \Delta A_p\, r_p^2\, P_{n,p}}{\sum_p \Delta A_p\, P_{n,p}},$$

- `P_n` is the code's `dPz` (time-averaged Poynting flux along lab z, **not** `dSz`, which is the spin density),
  times the sign of the screen normal's z component, so a backward screen gives a positive total.
- `ΔA_p` are `get_screen_area_weights_au`'s trapezoidal weights. The circular grid includes both endpoints in
  `u = r²` and in `φ` (`0` and `2π`), so the spec's equal-weight "cell centre" sum does not apply.
- Works for circular and rectangular screens facing ±Oz (`r² = x² + y²` in the screen's local frame).
- For a Gaussian `|E|² ∝ exp(−2r²/w²)`, `⟨r²⟩ = w²/2`, hence the factor 2.

### Inversion

With `A = w_eff²`, `B = λ_s |Z| / π`, `λ_s = 2πc/ω` (`ω = (ω/ω_1)·fundamental_frequency` from `run_log.txt`):

$$w_{0,\rm far}^2 = \frac{2B^2}{A+\sqrt{A^2-4B^2}}\quad (Z > z_R),\qquad
w_{0,\rm near}^2 = \frac{A+\sqrt{A^2-4B^2}}{2}\quad (Z < z_R),\qquad z_R = \pi w_0^2/\lambda_s .$$

Both roots are reported; the data on one screen cannot choose between them. The wavefront curvature
`R(Z) = Z(1 + (z_R/Z)²)` could (the residual phase after `radial_phase` processing), not done yet.

### Reported columns and status

Lengths in units of the laser wavelength `λ_L` (caption gives `λ_L` in nm and `|Z|` in mm):
`lambda_s, w_eff, w_eff/|Z|, w0_far, z_R_far, Z/z_R_far, w0_near, z_R_near, Z/z_R_near, edge_frac, status`.

- `edge_frac`: share of the flux at `r ≥ 0.9 ×` the screen's inscribed radius. Large means the beam is clipped by
  the screen (w_eff too small) or a broad background dominates the second moment (w_eff too large).
- `status`:
  - `non-positive or non-finite flux/second moment`;
  - `Gaussian model incompatible`: `A² − 4B² < −1e-10 A²` (not clamped; smaller negatives are rounding → 0);
  - `roots close (Z/z_R,far < 2): ill-conditioned`: near `Z ≈ z_R` the discriminant → 0 and the waist is very
    sensitive to noise in `w_eff`.

### Validity

The model assumes a fundamental Gaussian (`M² = 1`) with its waist at the origin. It holds only for:

- a Gaussian incident laser (`laser_type laguerre_gauss`, `p = l = 0`); the section warns otherwise;
- the first harmonic. Harmonic `N` carries `L_N = N(m+σ) − σ` (CLAUDE.md), so with circular polarization
  `N ≥ 2` is an LG mode with `|l| = N − 1`. The section always carries this warning.

### First result

Run `20261025_002302` (Gaussian laser, `w0 = 25 λ_L`, circular, `γ ≈ 13.5`, 8192 electrons,
`|Z| = 1.6×10⁷ λ_L`, `R_max = 40 w0`), lengths in `λ_L`:

| N | w_eff | w0_far | Z/z_R,far | w0_near | Z/z_R,near | edge_frac |
|---|-------|--------|-----------|---------|------------|-----------|
| 1 | 416   | 16.9   | 24.6      | 415     | 0.04       | 0.4%      |
| 2 | 647   | 5.4    | 119       | 647     | 0.008      | 1.6%      |

N = 1: the far root is the physical one (deep far field). The inferred `w0 ≈ 16.9 λ_L` is smaller than the
laser's 25 `λ_L`, **not yet explained**. Candidate: the incoherent background (∝ `N_e`, broad) widens `w_eff` and
lowers the inferred waist. Test: two runs with different `beam_particle_count` (and seeds); the coherent part
scales as `N_e²`, the background as `N_e`. N = 2 is meaningless here (circular, `|l| = 1`).

## Stage 2: least-squares fit on LG modes at fixed l

Spec: `theory/lg_screen_decomposition_simple.md`. Implemented in `py_scripts/radiation/lg_decomposition.py`.
Implementation, conventions, verification and first results: `theory/lg_decomposition_implementation.md`.
It fits the helicity components `E_±` with unit-norm `LG_p^l` (waist at the origin, `l = J ∓ 1`,
`J = N(m+σ)` about lab +z) over a scan of trial `w0`. The fit keeps the complex phase, so it uses the wavefront
curvature that Stage 1 left unused. First results: N = 1 has `eps = 5%` at `w0 ≈ 15.8 λ_L`, with 31% of the
power in `p = 1`; N = 2 has 9% of its power at `l = 0`.

## Next stages (to do)

- **LG second moment:** for a pure `LG_p^l`, `w_eff² = (2p+|l|+1) w(Z)²` (`M² = 2p+|l|+1`); use the expected
  `l = L_N` per harmonic as a first generalization.
- **Waist position:** add `z_0` to the Stage 2 basis (`s = Z − z_0`) and fit `w_0, z_0` by minimizing `eps`.
- Separate coherent and incoherent contributions before fitting.
