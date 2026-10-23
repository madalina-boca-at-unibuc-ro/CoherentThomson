# Coding task: stabilize the external propagation phase for an axis-aligned detector only

## Goal and strict scope

Make a **minimal change** to the existing `long_distance_simplified` and `long_distance_direct` radiation implementations (including their `_approx` variants, where applicable). Change **only the constant propagation phasor outside the proper-time integral**. Do not refactor the detector geometry or introduce new coordinate rotations.

Apply this change **only when the detector plane is orthogonal to the laboratory `Oz` axis** (the screen is a plane of constant laboratory `z = Z`, including `Z < 0` for backscattering). For any tilted detector, **retain the current implementation exactly**, using `exp(i k R_j)`.

## Existing expression

For each electron $j$ and detector point, let $\mathbf r_{j0}$ be the reference position (the first stored trajectory sample) and $\mathbf X=(X,Y,Z)$ the screen point in laboratory coordinates. The existing external factor is

$$
\frac{e^{ikR_j}}{R_j},\qquad R_j=|\mathbf X-\mathbf r_{j0}|.
$$

The integration over proper time $\tau$, the transverse projections, the integration-by-parts boundary term, and the internal trajectory phase **must remain unchanged**.

## Change for a detector orthogonal to `Oz`

Remove the screen-independent phase $e^{ik|Z|}$ by replacing only

$$
e^{ikR_j}\longrightarrow e^{ik(R_j-|Z|)}.
$$

**Use $|Z|$, not $Z$.** In the backscattering configuration the screen may be at negative $Z$, but the reference propagation distance is positive.

To calculate the phase argument stably, avoid directly subtracting $R_j-|Z|$ and avoid subtracting two large squares. With $\mathbf r_{j0}=(x_j,y_j,z_j)$ and the screen at $(X,Y,Z)$, evaluate

$$
\boxed{\displaystyle
\Delta R_j=R_j-|Z|=
\frac{(X-x_j)^2+(Y-y_j)^2+z_j^2-2Zz_j}{R_j+|Z|}.}
$$

This identity is exact for both positive and negative $Z$. The existing $R_j$ and $1/R_j$ remain unchanged.

If the current code already has an unrotated screen-point coordinate and electron reference position in the same laboratory frame, use those existing values. Do **not** add a coordinate transformation just for this change.

Keep the initial-time phase, if factored separately by the current code, exactly as before. In particular, do not merge large and small phase arguments unnecessarily. The external propagation phasor becomes `exp(i*k*DeltaR_j)` rather than `exp(i*k*R_j)`.

## Detector-orientation guard

Use the code's **existing** detector orientation settings or established axis-alignment test. Apply the new phasor only when the detector is unambiguously perpendicular to laboratory `Oz` (including either normal direction). Do not infer alignment merely from a small numerical angle unless an appropriate tolerance is already defined. If uncertain, use the original path.

For tilted detectors, do not subtract $|Z|$, do not project it with $\cos\theta_d$, and do not alter the existing phase calculation.

## Do not change

- The proper-time integration, trapezoidal weights, or sample endpoints.
- The internal Fourier phase, stepped-phase recurrence, or re-anchoring.
- The boundary contribution from integration by parts.
- The two transverse projections or reconstruction of $F^{\alpha\beta}$.
- The electron-dependent direction $\mathbf n_j$, amplitude $1/R_j$, or the electron reference position.
- Detector coordinate rotations, detector orientation conventions, or other radiation formulas.

## Validation

1. **Axis-aligned forward and backward detectors:** compare old and new complex fields using

   $$
   F_{\mathrm{old}}(\mathbf X)=e^{ik|Z|}F_{\mathrm{new}}(\mathbf X)
   $$

   for each frequency (up to floating-point roundoff). In particular, test **negative $Z$**.
2. Verify unchanged intensities and energy fluxes. Check any angular-momentum diagnostics that use the field; the removed phase is constant across an axis-aligned screen.
3. **Tilted detector:** verify that the original code path and its results remain unchanged.
4. Test a source reference position with nonzero $z_j$ to check the sign of the $-2Zz_j$ term.
5. Keep the patch small and report exactly which function/lines changed. Do not make unrelated optimizations.

## Rationale

The propagation phase $kR_j$ can be extremely large at distant screens. Subtracting the constant reference distance $|Z|$ before evaluating the complex exponential reduces its argument without changing relative phases between electrons or points on an axis-aligned detector. This is **not** a new far-field approximation: the expression for $\Delta R_j$ is an exact algebraic identity.
