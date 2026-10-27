# Task: Estimate Gaussian waist from detector Poynting flux

Implement this **only in the Python post-processing** of a single simulation run. Do not change the C++ Faraday-tensor export or the field calculations. Each run has one fixed detector distance `Z`; do not implement comparisons between runs or convergence tests in `Z`.

## 1. Effective beam radius on the detector

Reuse the existing computed, time-averaged **normal Poynting flux** `S_z` (not `|E_x|²` or the magnitude of the Poynting vector). The detector is a flat plane normal to the `z` axis, with a polar grid uniform in `u = r²` and azimuth `phi`.

Define

\[
w_{\mathrm{eff}}^2(Z)=2\frac{\int_0^{2\pi}d\phi\int du\;u\,S_z(\sqrt{u},\phi,Z)}{\int_0^{2\pi}d\phi\int du\;S_z(\sqrt{u},\phi,Z)}.
\]

For equal-area grid **cell centers**, use

\[
w_{\mathrm{eff}}^2(Z)\approx 2\frac{\sum_{i,j}r_i^2 S_{z,ij}}{\sum_{i,j}S_{z,ij}}.
\]

If samples include integration endpoints or have unequal quadrature weights, use the appropriate weights instead. If detector radii are stored as angles `theta = r/Z`, first calculate `theta_eff² = 2 sum(theta_i² S_z)/sum(S_z)` with the same quadrature weights, then `w_eff = |Z| theta_eff`. Keep units consistent. Validate that the total forward flux and the resulting second moment are positive and finite.

## 2. Infer the waist at the origin

Assume a fundamental Gaussian beam with waist at `z = 0` and emitted wavelength `lambda_s` obtained from the run's emitted angular frequency (use the project's existing units/conversions):

\[
w_{\mathrm{eff}}^2(Z)=w_0^2+\left(\frac{\lambda_s Z}{\pi w_0}\right)^2.
\]

Set `A = w_eff²` and `B = lambda_s * |Z| / pi`. Solve for the **smaller** root (the branch with `Z > z_R`):

\[
w_0^2=\frac{A-\sqrt{A^2-4B^2}}{2}
=\frac{2B^2}{A+\sqrt{A^2-4B^2}},
\]

using the second, numerically stable expression. If `A² - 4B² < 0` beyond floating-point tolerance, report that the assumed Gaussian model is incompatible with the measured width; do not silently clamp a significant negative discriminant. Calculate `z_R = pi*w_0²/lambda_s` as an additional diagnostic.

## 3. Output and scope

For this run, report `Z`, `lambda_s`, `w_eff(Z)`, `w_eff(Z)/|Z|`, inferred `w_0`, and `z_R`, with clear units, in the existing Python post-processing output. Keep the implementation small and consistent with the existing code structure. Do not implement LG projections, waist optimization, or multi-distance analysis yet.
