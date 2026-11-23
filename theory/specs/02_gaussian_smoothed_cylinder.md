# Implementation specification: Gaussian-smoothed electron cylinder

## Objective

Add an **optional longitudinal electron-position distribution** to the existing C++ Thomson-scattering simulation: a uniform cylinder smoothed along its longitudinal (`z`) axis by Gaussian displacements. Preserve the existing transverse electron sampling, laser model, particle dynamics, Fourier transforms, detector calculations, and radiation summation.

This is an initial-position distribution only; do **not** multiply the computed radiation by an additional form factor. The radiation code already sums individual electron contributions coherently. The analytical form factor below is for validation.

## 1. Mathematical definition

Let `H > 0` be the length of the underlying uniform longitudinal interval, `sigma_z >= 0` the Gaussian smoothing width (standard deviation), and `z_center` its center.

Generate each electron's initial longitudinal position as

$$
z_j=z_{\mathrm{center}}+u_j+g_j,
\qquad
u_j\sim\mathcal U(-H/2,H/2),
\qquad
g_j\sim\mathcal N(0,\sigma_z^2),
$$

with independent `u_j` and `g_j` for each electron. The two random variables must also be independent of the transverse coordinates unless the existing model explicitly imposes correlations.

The normalized longitudinal probability density, for `sigma_z > 0`, is

$$
\rho_z(z)=\frac{1}{2H}\left[
\operatorname{erf}\left(\frac{z-z_{\mathrm{center}}+H/2}{\sqrt2\sigma_z}\right)
-\operatorname{erf}\left(\frac{z-z_{\mathrm{center}}-H/2}{\sqrt2\sigma_z}\right)
\right].
$$

Equivalently, this is the convolution of a uniform density on `[z_center-H/2,z_center+H/2]` with a zero-mean Gaussian of standard deviation `sigma_z`.

Limits:
- `sigma_z = 0`: original sharp-edged uniform cylinder, exactly.
- `0 < sigma_z << H`: approximately flat central region with smoothly rounded edges and Gaussian tails.
- `sigma_z` comparable to or greater than `H`: broad, increasingly Gaussian-like profile; the meaning of `H` is the **underlying** top-hat width, not the total spatial extent of the sampled bunch.

The exact mean and variance are

$$
\mathbb E[z]=z_{\mathrm{center}},\qquad
\operatorname{Var}(z)=\frac{H^2}{12}+\sigma_z^2.
$$

The distribution has unbounded tails. Do not clip positions to `[-H/2,H/2]`, since clipping would destroy the intended density and form factor. If the existing code imposes a spatial cutoff, flag it and document any truncation.

## 2. C++ implementation

Locate the existing initialization routine that samples initial electron positions. Add a distribution choice, for example:

- `uniform_cylinder` (existing behavior)
- `gaussian_smoothed_cylinder` (new behavior)

Expose `H`, `sigma_z`, and (if already supported) `z_center` as configuration parameters in the **same length units used by the simulation**. Preserve the existing RNG and seed handling if possible; do not reset or reseed inside the electron loop.

Minimal illustrative code (adapt names to the existing codebase):

```cpp
#include <random>
#include <stdexcept>

// rng is the existing seeded random engine, passed by reference.
double sample_smoothed_z(std::mt19937_64& rng,
                         double H,
                         double sigma_z,
                         double z_center = 0.0)
{
    if (!(H > 0.0) || sigma_z < 0.0)
        throw std::invalid_argument("Require H > 0 and sigma_z >= 0");

    std::uniform_real_distribution<double> uniform(-0.5 * H, 0.5 * H);
    double z = z_center + uniform(rng);
    if (sigma_z > 0.0) {
        std::normal_distribution<double> gaussian(0.0, sigma_z);
        z += gaussian(rng);
    }
    return z;
}
```

Use the existing transverse sampler **unchanged**: a uniform disk of radius `R_e = 75 lambda_L` in the current setup. No additional transverse Gaussian weighting, acceptance/rejection, or radial form factor is requested.

Preserve the existing `uniform_cylinder` code path for backward compatibility. When `sigma_z = 0`, the new sampler must reproduce the sharp-cylinder **distribution**; exact electron-by-electron reproducibility with a prior RNG stream depends on the existing RNG call order.

## 3. Analytical longitudinal form factor (validation only)

With the convention

$$
F_z(q_z)=\int_{-\infty}^{\infty}\rho_z(z)e^{iq_z z}\,dz,
$$

the exact result is

$$
\boxed{
F_z(q_z)=e^{iq_z z_{\mathrm{center}}}
\operatorname{sinc}\left(\frac{q_zH}{2}\right)
\exp\left(-\frac{q_z^2\sigma_z^2}{2}\right),
}
$$

where `sinc(x) = sin(x)/x`, including `sinc(0) = 1`. Thus

$$
\boxed{
|F_z(q_z)|^2=
\operatorname{sinc}^2\left(\frac{q_zH}{2}\right)
\exp(-q_z^2\sigma_z^2).
}
$$

For electrons initially at rest, linear Thomson scattering, and exact backscattering, take

$$
q_z=k_L+k_s\simeq 2k_L=\frac{4\pi}{\lambda_L}.
$$

At a selected scattered frequency use `k_s = omega_s/c`, with the appropriate sign convention for the propagation direction. The exact simulation may deviate from the simple factorized form because of finite-pulse, focused-field, nonlinear, and off-axis effects. The analytical expression is **not** a substitute for summing the electron fields.

For fixed `H`, the ratio of coherent longitudinal intensity factors to the sharp-cylinder model, wherever the sharp-cylinder factor is nonzero, is

$$
\frac{|F_z(q_z;H,\sigma_z)|^2}{|F_z(q_z;H,0)|^2}
=\exp(-q_z^2\sigma_z^2).
$$

At `H = n lambda_L/2` in ideal backward linear scattering, the sinc factor is zero for **both** the sharp and smoothed distributions; the ratio is then undefined. Gaussian smoothing does not move the sinc zeros.

## 4. Suggested initial scan

Retain the current parameters unless the existing project configuration specifies otherwise:

- `lambda_L` approximately `0.8 micrometers`;
- `R_e = 75 lambda_L`;
- `a0 = 0.2`;
- electrons initially at rest;
- observation in the backward direction;
- initial electron count around `N = 200000` for comparison with the recent run.

Use `H = 0.87 lambda_L` first, because the existing sharp-cylinder run at that width produced a clear coherent pattern. Test

$$
\sigma_z/\lambda_L\in\{0,0.02,0.05,0.10,0.15,0.20,0.25,0.30\}.
$$

For the linear backward-scattering estimate, the suppression factor relative to the sharp case is

| `sigma_z/lambda_L` | `exp[-16 pi^2 (sigma_z/lambda_L)^2]` |
|---:|---:|
| 0.00 | 1 |
| 0.02 | 0.939 |
| 0.05 | 0.674 |
| 0.10 | 0.206 |
| 0.15 | 0.0286 |
| 0.20 | 0.00181 |
| 0.25 | 0.0000517 |
| 0.30 | 0.000000673 |

Keep `N` fixed while scanning `sigma_z`, and use multiple independent random seeds when estimating ensemble averages. The predicted suppression applies to the **coherent term**, not automatically to the total measured intensity, which includes a shot-noise contribution.

## 5. Validation tests and outputs

1. **Distribution test:** Generate many `z` samples. Verify sample mean approaches `z_center`, and sample variance approaches `H^2/12 + sigma_z^2`. Compare a normalized histogram against the analytic density. Confirm the presence of tails beyond the underlying top-hat boundaries.
2. **Zero-smoothing test:** With `sigma_z = 0`, verify that the longitudinal distribution is uniform over the original interval, and that existing physical outputs remain consistent with the original uniform-cylinder implementation (allowing for RNG ordering differences).
3. **Empirical form factor:** For each run, compute
   $$
   \widehat F_z(q_z)=\frac1N\sum_{j=1}^N e^{iq_z z_j}.
   $$
   Compare its complex value and squared magnitude with the analytical `F_z` and `|F_z|^2`. Across independent realizations,
   $$
   \mathbb E[|\widehat F_z|^2]
   =\frac1N+\left(1-\frac1N\right)|F_z|^2.
   $$
   This finite-sample shot-noise floor is essential when the analytical form factor is tiny.
4. **Radiation diagnostic:** Save the total complex Fourier field on the detector, the integrated spectral radiation intensity/energy, and a fixed-pixel value for each `(H, sigma_z, N, seed)`. Keep the existing Fourier and radiation computation unchanged.
5. **Independent-seed test:** Repeat selected cases with independently seeded electron populations. Do not infer deterministic coherence from the visual similarity of nested samples that share the same seed.
6. **Optional control:** Repeat at `H = lambda_L` to test whether smoothing reduces residual coherence from nonlinear frequency shifts; do not expect smoothing to remove an exact sinc zero or to create coherence at that zero.

## 6. Acceptance criteria

- The new longitudinal distribution can be selected from configuration without changing the existing transverse distribution or trajectory/radiation algorithms.
- `sigma_z = 0` recovers the uniform-cylinder distribution.
- Sample mean, variance, and empirical form factor agree with analytical expectations within statistical uncertainty.
- Results are reproducible for a fixed seed and vary appropriately across independent seeds.
- Units and the interpretation of `H` and `sigma_z` are documented in the project's configuration or README.
- No form-factor multiplier is inserted into the computed radiation: electron contributions are still summed individually.

## 7. Out of scope

Do not implement microbunching, transverse Gaussian electron density, self-consistent electron-electron interactions, plasma response, or modifications to the laser and detector geometry in this change. Those can be investigated separately after this benchmark is validated.
