# Coding task: single-wavevector electron form-factor diagnostic

Add a lightweight diagnostic to the existing C++ Thomson-scattering code. **Do not modify the radiation calculation or electron generator.** Run this diagnostic on the generated initial electron ensemble, before time integration.

## Define one momentum transfer per run

1. Calculate the **mean initial electron momentum vector** over all generated electrons, and convert it to the corresponding mean velocity \(\boldsymbol\beta=\mathbf v/c\) (using the relativistic momentum–velocity relation). Use this same velocity for the entire diagnostic.
2. Use the **central direction of the detector/screen**, \(\mathbf n_s\), including any detector tilt, as the scattered-photon direction. Do not recalculate \(\mathbf q\) for each screen pixel.
3. Use the incident Thomson laser propagation direction \(\mathbf n_L\) and its central angular frequency \(\omega_L\). **Compute the scattered frequency using Lorentz four-vector scalar products, as already implemented elsewhere in the code.** With metric \((+,-,-,-)\), define the mean initial electron four-momentum \(P^\mu=(E/c,\mathbf p)\) and the photon four-wavevectors

   \[
   K_L^\mu=\frac{\omega_L}{c}(1,\mathbf n_L),\qquad
   K_s^\mu=\frac{\omega_s}{c}(1,\mathbf n_s).
   \]

   For the classical Thomson fundamental (neglecting photon recoil), impose the Lorentz-invariant relation

   \[
   K_L\cdot P=K_s\cdot P.
   \]

   To obtain \(\omega_s\) explicitly, introduce the null direction four-vectors \(N_L^\mu=(1,\mathbf n_L)\) and \(N_s^\mu=(1,\mathbf n_s)\):

   \[
   \boxed{\omega_s=\omega_L\frac{N_L\cdot P}{N_s\cdot P}.}
   \]

   **Reuse the existing four-vector and Minkowski scalar-product operations in the C++ code** rather than implementing a separate three-vector Doppler formula. This is equivalent to the usual Doppler relation but matches the code's existing conventions. Use the mean initial electron momentum to construct \(P\), and use the screen-center direction for \(N_s\).

4. Define the single run-level momentum-transfer **wavevector** (units of inverse length):

   \[
   \boxed{\mathbf q=\mathbf k_s-\mathbf k_L
   =\frac{\omega_s}{c}\mathbf n_s-\frac{\omega_L}{c}\mathbf n_L.}
   \]

   Use the code's consistent internal units. The sign convention does not affect the magnitude of the form factor, but retain and report this convention explicitly.

## Compute the form factor

Using the **initial electron positions** \(\mathbf r_j(0)\), evaluate

\[
A(\mathbf q)=\sum_{j=1}^{N}e^{i\mathbf q\cdot\mathbf r_j(0)},
\qquad
b(\mathbf q)=\frac{A(\mathbf q)}{N}.
\]

Report the complex value \(b\) (real and imaginary parts), its magnitude \(|b|\), and \(|b|^2\). Also report \(N\), the mean initial momentum/velocity, \(\omega_s\), \(\mathbf k_L\), \(\mathbf k_s\), and \(\mathbf q\), with units. The result should be written to the existing run log (and optionally a small diagnostic CSV). No per-electron output is needed.

## Validation

- All electrons at the same initial position: \(|b|=1\).
- Flat disk at \(z=0\), with head-on illumination and on-axis detection: \(|b|\approx1\) because \(q_x=q_y=0\).
- Gaussian longitudinal positions of rms width \(\sigma_z\), with \(q_x=q_y=0\): compare with \(|b|\approx\exp[-q_z^2\sigma_z^2/2]\), allowing for finite-\(N\) fluctuations.
- Test the existing cases \(N=8000,16000\), \(\sigma_z=0\) and \(8.82\times10^{-4}\lambda_T\), using the same run-level \(\mathbf q\) definition.

**Scope:** This is a diagnostic of the initial spatial form factor at the nominal Thomson fundamental, not a replacement for the full radiation calculation or a screen-dependent phase-matching calculation.
