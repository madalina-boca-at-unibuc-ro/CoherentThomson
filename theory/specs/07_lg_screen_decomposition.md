# Coding-agent task: first LG decomposition of one emitted-field component

## Scope — keep this first implementation simple

Implement a **post-processing** decomposition of **one selected complex Cartesian field component** of **one emitted harmonic**, on the **existing detector screen at its existing fixed distance**. Do not change the Thomson calculation, detector geometry, propagation, or radiation output.

For each trial waist `w0` listed in the existing input-file format, use the same `w0` for all modes in that particular expansion. Fit modes with radial indices `p_min <= p <= p_max` and one fixed azimuthal index. Save the **complex coefficients** and **one reconstruction error** for each trial waist. **Do not optimize** the waist, fit several components, propagate to another detector, or implement additional modal-purity diagnostics yet.

## Indices and input

Use the notation `p` for the radial index and `l` for the azimuthal index throughout.

For harmonic order `N`, incident azimuthal index `l_incident`, and circular polarization helicity `sigma = +/-1`, fix

\[
l_N=N(l_{\mathrm{incident}}+\sigma)-\sigma.
\]

Use this `l_N` for **every** mode in the expansion; do not vary it during fitting. This relation presumes the emitted polarization/helicity component follows the same-helicity convention. If the code's field components or helicity conventions differ, document the discrepancy rather than silently changing the formula.

Read from the existing input system (reuse existing parameters where available):

- `harmonic_N`, `l_incident`, `sigma` (or their existing equivalents);
- `field_component`: **one** component, for example `Ex` (not total intensity or a sum of components);
- `p_min`, `p_max`: inclusive integers, `0 <= p_min <= p_max`;
- `w0_values`: one or more positive trial waists, with units handled by the existing units parser.

Use the existing screen coordinates, integration weights, emitted harmonic wavelength, and output directory. For this first test, assume the LG waist plane is at the coordinate origin unless the project's geometry defines a different reference plane; state the convention explicitly.

## Explicit LG basis: normalized transverse modes

**Do not assume the LG fields already implemented for the incident laser have the normalization required here.** They may use a peak-field or other convention. Either implement a separate normalized basis function, or verify and correct the prefactor before using the existing function.

Let `s` denote distance measured **along the emitted beam's propagation direction** from its waist plane. In particular, for backscattering along negative laboratory `Z`, do **not** simply use a positive-laboratory-`Z` propagation formula unchanged. Define

\[
 k=\frac{2\pi}{\lambda_N},\qquad z_R=\frac{k w_0^2}{2}=\frac{\pi w_0^2}{\lambda_N},
\]
\[
 w(s)=w_0\sqrt{1+(s/z_R)^2},\quad
 R_c(s)=s\left[1+(z_R/s)^2\right],\quad
 \psi(s)=\arctan(s/z_R).
\]

At `s=0`, interpret `1/R_c(0)=0` by continuity.

With the complex time convention `exp(-i omega t)` and transverse azimuth `phi` defined in the **right-handed coordinate system aligned with the emitted propagation direction**, use

\[
\boxed{
\begin{aligned}
 u_{p,l}(r,\phi,s;w_0)={}&
 \sqrt{\frac{2p!}{\pi(p+|l|)!}}\,\frac{1}{w(s)}
 \left(\frac{\sqrt{2}r}{w(s)}\right)^{|l|}
 L_p^{|l|}\!\left(\frac{2r^2}{w^2(s)}\right)
 e^{-r^2/w^2(s)}\\
 &\times e^{i l\phi}
 e^{i k s}
 e^{i k r^2/[2R_c(s)]}
 e^{-i(2p+|l|+1)\psi(s)}.
\end{aligned}
}
\]

Here `L_p^{|l|}` is the generalized Laguerre polynomial. This convention satisfies the **unit transverse norm** on an infinite plane:

\[
 \int_0^{2\pi}\!d\phi\int_0^\infty\!r\,dr\;u_{p,l}^*(r,\phi,s;w_0)
 u_{q,m}(r,\phi,s;w_0)=\delta_{pq}\delta_{lm}.
\]

In particular, `u_00(r,0;w0) = sqrt(2/pi)/w0 * exp(-r^2/w0^2)` at the waist (up to the global carrier phase).

**Conventions to check against the code:** (1) sign of the azimuthal phase and definition of `phi`, (2) emitted propagation direction, (3) complex-field time convention, (4) the physical position of the waist plane, (5) whether the stored complex field has a global or transverse propagation phase removed. Use the **unmodified complex screen field** for this test. A spatially uniform carrier phase may be handled consistently in both target and basis if required, but do not flatten the transverse curvature phase. If the stored field is a laboratory Cartesian component, keep that same component throughout; do not silently reinterpret it as a helicity amplitude.

## Fit one complex component

For each trial `w0`, evaluate the modes `u_{p,l_N}` at the **actual detector plane** for all `p` from `p_min` through `p_max`.

Let `E_i` be the chosen complex emitted-field component at screen sample `i`, `A_i` the **existing physical area-integration weight**, and `B_{ip}=u_{p,l_N}(r_i,phi_i,s_screen;w0)`. Solve

\[
 \min_{\{c_p\}}\sum_i A_i\left|E_i-\sum_{p=p_{\min}}^{p_{\max}}c_p B_{ip}\right|^2.
\]

Use weighted **complex** least squares (`sqrt(A_i)` times rows of target and design matrix) with QR or SVD; **do not** assume exact LG orthogonality on the finite screen and do not explicitly invert the Gram matrix. Preserve both real and imaginary parts of each coefficient.

Reconstruct

\[
 E_i^{\mathrm{rec}}=\sum_{p=p_{\min}}^{p_{\max}}c_p B_{ip}
\]

and compute just the primary, area-weighted relative complex-field error

\[
\boxed{
 \epsilon(w_0)=\sqrt{\frac{\sum_i A_i|E_i-E_i^{\mathrm{rec}}|^2}
 {\sum_i A_i|E_i|^2}}.
}
\]

Use the existing screen-area weights, including the polar-grid Jacobian and endpoint conventions. If the chosen field component has zero norm, report that the relative error is undefined rather than dividing by zero.

## Minimal output

For **each** trial `w0`, write a small table containing:

- `w0`, `N`, `l_N`, selected field component, `p_min`, `p_max`;
- for each `p`: `Re(c_p)`, `Im(c_p)`, `abs(c_p)` (phase is optional);
- one value of `epsilon(w0)`.

Save the reconstructed complex component on the existing screen grid **only if this fits naturally into the existing output format**; it is useful for visual inspection but not a reason to build a new output subsystem. If the least-squares system is rank deficient or severely ill-conditioned, report a warning rather than interpreting the coefficients physically.

## Minimal verification before applying to Thomson data

1. Check the analytic prefactor numerically: on a sufficiently wide, adequately sampled plane, the area integral of `|u_{p,l}|^2` should be approximately 1 for a few low-order modes; check one off-diagonal overlap as well.
2. Construct a synthetic complex field from **one known normalized LG mode** on the same screen grid, fit it with its correct `w0` and `p` included, and verify that the coefficient and reconstruction error are recovered to numerical accuracy.
3. Then run the decomposition on the chosen Thomson harmonic and inspect the coefficient tables and errors. **Stop there.**

## Interpretation limits

The `w0` values label alternative **basis choices**; they are not assumed to be conserved from the incident laser or uniquely defined by the radiation. The coefficients have the normalization implied by the explicit unit-norm modes above. On a finite screen, `|c_p|^2` must not automatically be called an independently additive modal power. A nonzero residual can reflect the radial truncation, other azimuthal content, finite aperture, or deviations from the paraxial approximation. No multi-distance propagation test is requested.
