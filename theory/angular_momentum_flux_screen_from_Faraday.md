# Screen-level angular-momentum fluxes from the Faraday tensor

## Purpose

This note is an implementation specification for computing the one-sided
spectral flux densities of the third component of angular momentum through a
detector screen parallel to the $x_1x_2$ plane. The only field input is the
complex Fourier-transformed Faraday tensor
$F^{\alpha\beta}(\omega,x_1,x_2,x_3)$ already evaluated at the screen.

The required outputs at every screen point and every positive angular
frequency are

$$
\frac{d\mathcal M^{\mathrm{int}}_{33,+}}{d\omega},\qquad
\frac{d\mathcal M^{\mathrm{spin}}_{33,+}}{d\omega},\qquad
\frac{d\mathcal M^{\mathrm{orb,int}}_{33,+}}{d\omega}.
$$

All formulas are in SI units. They are exact and make no assumption about a
preferred propagation direction. The energy-flux centroid is assumed to be
at the origin:

$$
X_c(\omega)=Y_c(\omega)=0.
$$

Consequently, the intrinsic and total $J_3$ flux densities coincide.

## Fourier and tensor conventions

The forward transform used by the TeX document is

$$
F^{\alpha\beta}(\omega,{\bf x})
=\frac{1}{2\pi}\int_{-\infty}^{\infty}dt\,
e^{+i\omega t}F^{\alpha\beta}(t,{\bf x}),
$$

with inverse transform

$$
F^{\alpha\beta}(t,{\bf x})
=\int_{-\infty}^{\infty}d\omega\,
e^{-i\omega t}F^{\alpha\beta}(\omega,{\bf x}).
$$

The formulas below are one-sided and apply only for $\omega>0$. The factor
of $4\pi$ must not be removed: $2\pi$ comes from Parseval's identity for
this transform normalization, and a further factor of two converts the
two-sided spectrum of a real field to the one-sided spectrum.

Use the component conventions

$$
E_i=cF^{i0},\qquad
B_1=-F^{23},\qquad B_2=F^{13},\qquad B_3=-F^{12}.
$$

The tensor is antisymmetric, $F^{\beta\alpha}=-F^{\alpha\beta}$, and its six
independent complex components are

$$
F^{10},\ F^{20},\ F^{30},\ F^{12},\ F^{13},\ F^{23}.
$$

In all expressions below, `conj(z)` means complex conjugation,
`real(z)` means the real part, and `imag(z)` means the imaginary part.

## 1. Required momentum-flux components

Only $d\mathcal T_{13,+}/d\omega$ and
$d\mathcal T_{23,+}/d\omega$ are needed. Written explicitly in terms of the
six independent Faraday components, they are

$$
\boxed{
\frac{d\mathcal T_{13,+}}{d\omega}
=-4\pi\epsilon_0c^2\operatorname{Re}\!\left[
F^{10}(F^{30})^*+F^{23}(F^{12})^*
\right]
}
$$

and

$$
\boxed{
\frac{d\mathcal T_{23,+}}{d\omega}
=4\pi\epsilon_0c^2\operatorname{Re}\!\left[
-F^{20}(F^{30})^*+F^{13}(F^{12})^*
\right].
}
$$

Equivalent code-like expressions are

```text
T13 = -4*pi*epsilon0*c^2 * real(F10*conj(F30) + F23*conj(F12))
T23 =  4*pi*epsilon0*c^2 * real(-F20*conj(F30) + F13*conj(F12))
```

Here `T13` and `T23` denote the one-sided spectral quantities
$d\mathcal T_{13,+}/d\omega$ and
$d\mathcal T_{23,+}/d\omega$, not instantaneous stress-tensor components.

## 2. Intrinsic angular-momentum flux density

At a screen point $(x_1,x_2)$, with $X_c=Y_c=0$, compute

$$
\boxed{
\frac{d\mathcal M^{\mathrm{int}}_{33,+}}{d\omega}
=x_1\frac{d\mathcal T_{23,+}}{d\omega}
-x_2\frac{d\mathcal T_{13,+}}{d\omega}.
}
$$

The completely expanded Faraday-tensor expression is

$$
\boxed{
\begin{aligned}
\frac{d\mathcal M^{\mathrm{int}}_{33,+}}{d\omega}
=4\pi\epsilon_0c^2\Big\{&
x_1\operatorname{Re}\!\left[
-F^{20}(F^{30})^*+F^{13}(F^{12})^*
\right]\\
&+x_2\operatorname{Re}\!\left[
F^{10}(F^{30})^*+F^{23}(F^{12})^*
\right]\Big\}.
\end{aligned}
}
$$

Code-like form:

```text
M33_int = x1*T23 - x2*T13
```

Because the centroid is the origin, this also equals the total screen-level
flux density $d\mathcal M_{33,+}/d\omega$.

## 3. Spin angular-momentum flux density

Compute the spin contribution directly from the Faraday components as

$$
\boxed{
\frac{d\mathcal M^{\mathrm{spin}}_{33,+}}{d\omega}
=\frac{4\pi\epsilon_0c^3}{\omega}\operatorname{Im}\!\left[
-(F^{23})^*F^{10}
+(F^{13})^*F^{20}
+(F^{12})^*F^{30}
\right],\qquad \omega>0.
}
$$

Code-like form:

```text
spin_complex = -conj(F23)*F10 + conj(F13)*F20 + conj(F12)*F30
M33_spin = (4*pi*epsilon0*c^3/omega) * imag(spin_complex)
```

Do not evaluate this expression at $\omega=0$, where the explicit
$1/\omega$ factor requires a separate limiting analysis.

## 4. Intrinsic orbital angular-momentum flux density

The intrinsic orbital flux density is defined as the residual after
subtracting the spin flux density from the intrinsic total flux density:

$$
\boxed{
\frac{d\mathcal M^{\mathrm{orb,int}}_{33,+}}{d\omega}
=\frac{d\mathcal M^{\mathrm{int}}_{33,+}}{d\omega}
-\frac{d\mathcal M^{\mathrm{spin}}_{33,+}}{d\omega}.
}
$$

Code-like form:

```text
M33_orb_int = M33_int - M33_spin
```

This residual formula is the definition to use for the spatial screen map.
Do not replace it by $c$ times the azimuthal-derivative formula for the
orbital angular-momentum density. For a general field, the exact canonical
flux separation contains additional mixed electric--magnetic and spin--orbit
terms.

## 5. Integration over the screen

If the detector has pixels with areas $\Delta A_p$, integrate each local
flux density using the same pixel mask and quadrature weights:

$$
\frac{d\mathcal J^q_3}{d\omega}
=\int_A dA\,\frac{d\mathcal M^q_{33,+}}{d\omega}
\simeq\sum_p \Delta A_p
\left.\frac{d\mathcal M^q_{33,+}}{d\omega}\right|_p,
\qquad
q\in\{\mathrm{int},\mathrm{spin},\mathrm{orb,int}\}.
$$

The pointwise and integrated consistency checks are

```text
M33_int == M33_spin + M33_orb_int
J3_int   == J3_spin   + J3_orb_int
```

up to floating-point and quadrature error.

## 6. Minimal implementation order

For each positive frequency and every screen pixel:

1. Load `F10`, `F20`, `F30`, `F12`, `F13`, and `F23` as complex numbers.
2. Compute `T13` and `T23` from the explicit formulas in Section 1.
3. Compute `M33_int = x1*T23 - x2*T13`.
4. Compute `M33_spin` from the imaginary part in Section 3.
5. Compute `M33_orb_int = M33_int - M33_spin`.
6. Store the three real-valued spatial maps.
7. Integrate them over the identical screen grid if spectra
   $d\mathcal J_3/d\omega$ are required.

Useful numerical checks:

- Every final flux-density map should be real to numerical precision.
- Reconstructing `M33_int` from the fully expanded formula must agree with
  `x1*T23 - x2*T13`.
- The residual identity must hold pointwise before screen integration.
- Preserve the physical screen coordinates, including their length units;
  pixel indices must not be substituted for $x_1$ and $x_2$.
- Apply the same spatial mask to all three maps. An asymmetric mask can spoil
  the assumed condition $X_c=Y_c=0$.

