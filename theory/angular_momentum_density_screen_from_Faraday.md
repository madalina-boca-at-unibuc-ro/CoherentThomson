# Screen-level angular-momentum densities from the Faraday tensor

## Purpose

This note is an implementation specification for computing spatial maps of
the one-sided spectral density of the third component of electromagnetic
angular momentum on a detector screen parallel to the $x_1x_2$ plane.

The only field input is the complex Fourier-transformed Faraday tensor
$F^{\alpha\beta}(\omega,x_1,x_2,x_3)$, already evaluated at every screen
point. The required outputs are

$$
\frac{d\mathcal L^{\mathrm{int}}_{3,+}}{d\omega},\qquad
\frac{d\mathcal L^{\mathrm{spin}}_{3,+}}{d\omega},\qquad
\frac{d\mathcal L^{\mathrm{orb,res}}_{3,+}}{d\omega},\qquad
\frac{d\mathcal L^{\mathrm{orb,can}}_{3,+}}{d\omega}.
$$

All formulas are in SI units. They are exact and make no assumption about a
preferred propagation direction. The energy-flux centroid is assumed to be
at the origin:

$$
X_c(\omega)=Y_c(\omega)=0.
$$

Thus all moments and azimuthal derivatives below are taken about the
coordinate origin.

## Fourier and tensor conventions

The forward Fourier transform is

$$
F^{\alpha\beta}(\omega,{\bf x})
=\frac{1}{2\pi}\int_{-\infty}^{\infty}dt\,
e^{+i\omega t}F^{\alpha\beta}(t,{\bf x}),
$$

and the inverse transform is

$$
F^{\alpha\beta}(t,{\bf x})
=\int_{-\infty}^{\infty}d\omega\,
e^{-i\omega t}F^{\alpha\beta}(\omega,{\bf x}).
$$

All formulas below are one-sided and apply to $\omega>0$. The common
$4\pi$ factors in quadratic kinetic quantities contain the Parseval factor
$2\pi$ and the additional factor of two required to combine the positive and
negative frequencies of a real time-domain field. The canonical spin and
orbital expressions have the corresponding factor $2\pi/\omega$.

Use

$$
E_i=cF^{i0},\qquad
B_1=-F^{23},\qquad B_2=F^{13},\qquad B_3=-F^{12}.
$$

The Faraday tensor is antisymmetric,
$F^{\beta\alpha}=-F^{\alpha\beta}$, so the six independent complex inputs are

$$
F^{10},\ F^{20},\ F^{30},\ F^{12},\ F^{13},\ F^{23}.
$$

In the code-like expressions, `conj(z)`, `real(z)`, and `imag(z)` denote the
complex conjugate, real part, and imaginary part, respectively.

## 1. Intrinsic kinetic angular-momentum density

With the centroid at the origin, the third component is the moment of the
one-sided spectral linear-momentum density:

$$
\boxed{
\frac{d\mathcal L^{\mathrm{int}}_{3,+}}{d\omega}
=x_1\frac{d\mathcal G^2_+}{d\omega}
-x_2\frac{d\mathcal G^1_+}{d\omega}.
}
$$

The two required momentum-density components are

$$
\frac{d\mathcal G^1_+}{d\omega}
=-4\pi\epsilon_0c\operatorname{Re}\!\left[
F^{20}(F^{12})^*+F^{30}(F^{13})^*
\right],
$$

$$
\frac{d\mathcal G^2_+}{d\omega}
=4\pi\epsilon_0c\operatorname{Re}\!\left[
F^{10}(F^{12})^*-F^{30}(F^{23})^*
\right].
$$

Therefore the completely expanded expression is

$$
\boxed{
\begin{aligned}
\frac{d\mathcal L^{\mathrm{int}}_{3,+}}{d\omega}
=4\pi\epsilon_0c\Big\{&
x_1\operatorname{Re}\!\left[
F^{10}(F^{12})^*-F^{30}(F^{23})^*
\right]\\
&+x_2\operatorname{Re}\!\left[
F^{20}(F^{12})^*+F^{30}(F^{13})^*
\right]\Big\}.
\end{aligned}
}
$$

Code-like form:

```text
G1 = -4*pi*epsilon0*c * real(F20*conj(F12) + F30*conj(F13))
G2 =  4*pi*epsilon0*c * real(F10*conj(F12) - F30*conj(F23))

L3_int = x1*G2 - x2*G1
```

This is the intrinsic **kinetic** density because it is constructed from the
Poynting momentum density ${\bf g}={\bf S}/c^2$.

## 2. Spin angular-momentum density

The dual-symmetric spin density is

$$
\boxed{
\begin{aligned}
\frac{d\mathcal L^{\mathrm{spin}}_{3,+}}{d\omega}
=\frac{2\pi\epsilon_0c^2}{\omega}\operatorname{Im}\Big[&
(F^{10})^*F^{20}-(F^{20})^*F^{10}\\
&+(F^{13})^*F^{23}-(F^{23})^*F^{13}
\Big],\qquad \omega>0.
\end{aligned}
}
$$

Code-like form:

```text
spin_complex = (conj(F10)*F20 - conj(F20)*F10
                + conj(F13)*F23 - conj(F23)*F13)

L3_spin = (2*pi*epsilon0*c^2/omega) * imag(spin_complex)
```

Do not evaluate this formula at $\omega=0$. The zero-frequency contribution,
if needed, requires a separate limiting analysis.

## 3. Orbital density, form A: residual

One local orbital representative is obtained by subtracting the spin density
from the intrinsic kinetic density:

$$
\boxed{
\frac{d\mathcal L^{\mathrm{orb,res}}_{3,+}}{d\omega}
=\frac{d\mathcal L^{\mathrm{int}}_{3,+}}{d\omega}
-\frac{d\mathcal L^{\mathrm{spin}}_{3,+}}{d\omega}.
}
$$

Code-like form:

```text
L3_orb_res = L3_int - L3_spin
```

This definition enforces the pointwise identity

$$
\frac{d\mathcal L^{\mathrm{int}}_{3,+}}{d\omega}
=\frac{d\mathcal L^{\mathrm{spin}}_{3,+}}{d\omega}
+\frac{d\mathcal L^{\mathrm{orb,res}}_{3,+}}{d\omega}.
$$

## 4. Orbital density, form B: canonical azimuthal derivative

Introduce polar coordinates on the screen,

$$
x_1=\rho\cos\phi,\qquad
x_2=\rho\sin\phi,
$$

so that the generator of rotations about the third axis is

$$
\mathscr L_3=x_1\partial_2-x_2\partial_1
=\frac{\partial}{\partial\phi}.
$$

The dual-symmetric canonical orbital density is

$$
\boxed{
\frac{d\mathcal L^{\mathrm{orb,can}}_{3,+}}{d\omega}
=\frac{2\pi\epsilon_0c^2}{\omega}\operatorname{Im}\!\left[
\sum_{a=1}^{3}(F^{a0})^*\frac{\partial F^{a0}}{\partial\phi}
+\sum_{1\leq a<b\leq3}(F^{ab})^*
\frac{\partial F^{ab}}{\partial\phi}
\right],\qquad \omega>0.
}
$$

Expanded without index sums:

$$
\boxed{
\begin{aligned}
\frac{d\mathcal L^{\mathrm{orb,can}}_{3,+}}{d\omega}
=\frac{2\pi\epsilon_0c^2}{\omega}\operatorname{Im}\Big[&
(F^{10})^*\partial_\phi F^{10}
+(F^{20})^*\partial_\phi F^{20}
+(F^{30})^*\partial_\phi F^{30}\\
&+(F^{12})^*\partial_\phi F^{12}
+(F^{13})^*\partial_\phi F^{13}
+(F^{23})^*\partial_\phi F^{23}\Big].
\end{aligned}
}
$$

Code-like form:

```text
orb_can_complex = (conj(F10)*dphi_F10
                   + conj(F20)*dphi_F20
                   + conj(F30)*dphi_F30
                   + conj(F12)*dphi_F12
                   + conj(F13)*dphi_F13
                   + conj(F23)*dphi_F23)

L3_orb_can = (2*pi*epsilon0*c^2/omega) * imag(orb_can_complex)
```

The derivative is taken at fixed $\rho$, fixed $x_3$, and fixed $\omega$.

### Computing the azimuthal derivative on a Cartesian screen grid

It is not necessary to interpolate the field onto a polar grid. For every
component separately, use

$$
\boxed{
\partial_\phi F^{\alpha\beta}
=-x_2\,\partial_{x_1}F^{\alpha\beta}
+x_1\,\partial_{x_2}F^{\alpha\beta}.
}
$$

Code-like form:

```text
dphi_F = -x2*d_dx1(F) + x1*d_dx2(F)
```

Apply the spatial derivatives to the complex field itself, not separately to
its magnitude or phase. Use the same differentiation scheme for all six
components. At $\rho=0$, the Cartesian form above is preferable because the
polar angle is undefined; its exact value is zero when the Cartesian
derivatives remain finite.

## 5. Relation between the two orbital forms

The residual and canonical expressions are two different local
representatives:

```text
L3_orb_res  = L3_int - L3_spin
L3_orb_can  = azimuthal-derivative expression
```

They need not agree point by point. The kinetic and canonical momentum
densities differ by a spatial-divergence term. Their integrals over a complete
screen agree when the corresponding boundary term vanishes. Therefore:

- use `L3_orb_res` when an exact pointwise decomposition of `L3_int` is
  required;
- use `L3_orb_can` when the canonical orbital density associated with the
  azimuthal structure of the Fourier field is required;
- do not use pointwise equality between these two maps as a numerical test;
- if their screen integrals are compared, also test convergence with respect
  to screen size and inspect the field at the boundary.

## 6. Screen integration

For pixels with areas $\Delta A_p$, a screen integral of any density map is

$$
\int_A dA\,\frac{d\mathcal L^q_{3,+}}{d\omega}
\simeq
\sum_p\Delta A_p
\left.\frac{d\mathcal L^q_{3,+}}{d\omega}\right|_p,
$$

where $q$ can be `int`, `spin`, `orb,res`, or `orb,can`. These are integrals of
the angular-momentum **density on the screen**, not angular-momentum fluxes
through the screen. Fluxes $d\mathcal M_{33,+}/d\omega$ are specified in
`angular_momentum_flux_screen_from_Faraday.md`.

Use the same mask and quadrature weights for every map. An asymmetric mask can
invalidate the assumed condition $X_c=Y_c=0$.

## 7. Minimal implementation order

For every positive frequency:

1. Load the six complex Faraday components on the complete screen grid.
2. Compute `G1` and `G2`.
3. Compute `L3_int = x1*G2 - x2*G1`.
4. Compute `L3_spin`.
5. Compute `L3_orb_res = L3_int - L3_spin`.
6. Differentiate all six complex Faraday-component maps with respect to
   $x_1$ and $x_2$.
7. Form each `dphi_F` and compute `L3_orb_can`.
8. Store all four real-valued maps.
9. If required, integrate the maps using identical screen weights.

## 8. Numerical checks

- `L3_int`, `L3_spin`, `L3_orb_res`, and `L3_orb_can` must be real-valued.
- Verify `L3_int == L3_spin + L3_orb_res` pointwise to floating-point
  precision.
- Verify that the compact and expanded formulas for `L3_int` agree.
- Verify that the sum and expanded forms of `L3_orb_can` agree.
- Check convergence of $\partial_\phi F$ under grid refinement and with at
  least two spatial-differentiation schemes when possible.
- Do not differentiate across masked pixels or screen boundaries without an
  explicit one-sided or boundary stencil.
- Use physical $x_1,x_2$ coordinates with length units, not integer pixel
  indices.
- Apply the same Fourier normalization, units, component ordering, and screen
  mask at every frequency.

