# Screen-level angular-momentum fluxes and densities from the Faraday tensor

## 1. Purpose

This note is an implementation specification for computing the one-sided spectral angular momentum densities and fluxes (specifically the third component along the $z$-axis, $J_3$) across a detector screen parallel to the $x_1x_2$ ($xy$) plane at fixed $x_3 = z$.

The input field is the complex Fourier-transformed Faraday tensor $F^{\alpha\beta}(\omega, x_1, x_2, x_3)$ evaluated at the detector plane.

All formulas are in SI units. We assume the origin of coordinates in the screen plane satisfies $X_c(\omega) = Y_c(\omega) = 0$, so intrinsic and extrinsic components around the $z$-axis coincide.

---

## 2. Fourier and Tensor Conventions

### 2.1 Fourier Transform Convention
The forward and inverse temporal Fourier transforms are defined as:
$$
F^{\alpha\beta}(\omega, \mathbf{x}) = \frac{1}{2\pi}\int_{-\infty}^{\infty} dt\, e^{i\omega t} F^{\alpha\beta}(t, \mathbf{x}),
$$
$$
F^{\alpha\beta}(t, \mathbf{x}) = \int_{-\infty}^{\infty} d\omega\, e^{-i\omega t} F^{\alpha\beta}(\omega, \mathbf{x}).
$$

### 2.2 Parseval / One-Sided Spectral Normalization
For two real physical time-dependent fields $f(t)$ and $g(t)$, Parseval's identity under this definition gives:
$$
\int_{-\infty}^{\infty} f(t) g(t)\, dt = 2\pi \int_{-\infty}^{\infty} f(\omega) g^*(\omega)\, d\omega = 4\pi \operatorname{Re} \int_{0}^{\infty} f(\omega) g^*(\omega)\, d\omega.
$$
Therefore, for any bilinear product in time domain $A(t) = f(t)g(t)$, the one-sided spectral density $\frac{dA_+}{d\omega}$ (such that $\int_{-\infty}^\infty A(t)dt = \int_0^\infty \frac{dA_+}{d\omega} d\omega$) carries a global prefactor of:
$$
C_{\text{Parseval}} = 4\pi.
$$

### 2.3 Field Tensor Identifications (SI)
With metric signature $(+,-,-,-)$ or standard 4-vector notation:
$$
E_i = c F^{i0}, \qquad i \in \{1, 2, 3\},
$$
$$
B_1 = -F^{23}, \qquad B_2 = F^{13}, \qquad B_3 = -F^{12}.
$$
The tensor is antisymmetric ($F^{\beta\alpha} = -F^{\alpha\beta}$), with six independent complex components:
$$
F^{10},\ F^{20},\ F^{30},\ F^{12},\ F^{13},\ F^{23}.
$$

---

## 3. Angular Momentum Densities along $Oz$

In the time domain, the linear momentum volume density is the Poynting momentum $\mathbf{p} = \varepsilon_0 (\mathbf{E} \times \mathbf{B})$. The corresponding physical volume angular momentum density around the $z$-axis is:
$$
j_3^{\text{direct}} = (\mathbf{r} \times \mathbf{p})_3 = x p_y - y p_x = \varepsilon_0 \left[ x (\mathbf{E} \times \mathbf{B})_y - y (\mathbf{E} \times \mathbf{B})_x \right].
$$

In the frequency domain (one-sided spectrum for $\omega > 0$):

### 3.1 Direct Density ($j_3^{\text{direct}}$)
$$
\frac{d j_3^{\text{direct}}}{d\omega} = 4\pi\varepsilon_0 \operatorname{Re}\left[ x \left( E_3 B_1^* - E_1 B_3^* \right) - y \left( E_2 B_3^* - E_3 B_2^* \right) \right].
$$

---

### 3.2 Canonical Spin and Orbital Densities
In monochromatic radiation fields, the linear Poynting momentum decomposes into canonical (orbital) momentum and spin momentum:
$$
\mathbf{p} = \mathbf{p}^{\text{orb}} + \mathbf{p}^{\text{spin}}, \qquad \text{where } \mathbf{p}^{\text{spin}} = \frac{1}{2} \nabla \times \mathbf{s}.
$$
Here $\mathbf{s} = (s_x, s_y, s_z)$ is the physical spin density vector:
$$
\mathbf{s} = \frac{4\pi\varepsilon_0}{4\omega} \operatorname{Im}\left( \mathbf{E}^* \times \mathbf{E} + c^2 \mathbf{B}^* \times \mathbf{B} \right).
$$

Componentwise, the spin density along the $z$-axis is:
$$
\frac{d j_3^{\text{spin}}}{d\omega} = s_z = \frac{4\pi\varepsilon_0}{4\omega} \operatorname{Im}\left[ \left(E_1^* E_2 - E_2^* E_1\right) + c^2 \left(B_1^* B_2 - B_2^* B_1\right) \right].
$$
*(Note: Since $u^* v - v^* u = 2i \operatorname{Im}(u^* v)$, the term $\operatorname{Im}(E_1^* E_2 - E_2^* E_1) = 2 \operatorname{Re}(E_1^* E_2)$ depending on how the cross product is written; retaining the anti-symmetrized form with prefactor $\frac{4\pi\varepsilon_0}{4\omega}$ ensures exact consistency).*

The orbital angular momentum density is:
$$
\frac{d j_3^{\text{orbital}}}{d\omega} = [\mathbf{r} \times \mathbf{p}^{\text{orb}}]_3 = \frac{4\pi\varepsilon_0}{4\omega} \operatorname{Im}\left[ \sum_{k=1}^3 \left( E_k^* \partial_\phi E_k + c^2 B_k^* \partial_\phi B_k \right) \right],
$$
where $\partial_\phi = x \partial_y - y \partial_x$.

---

## 4. Exact Pointwise Identity and the Divergence Correction Term

The relation between the direct density and the canonical spin + orbital split is **not** a local equality; it differs by a curl/divergence:
$$
\mathbf{r} \times \mathbf{p} = \mathbf{r} \times \mathbf{p}^{\text{orb}} + \frac{1}{2} \mathbf{r} \times (\nabla \times \mathbf{s}).
$$

Using the vector identity for any vector field $\mathbf{s}$:
$$
\left[ \mathbf{r} \times (\nabla \times \mathbf{s}) \right]_z = 2 s_z + \frac{\partial}{\partial z}\left( x s_x + y s_y \right) - \left[ \frac{\partial}{\partial x}(x s_z) + \frac{\partial}{\partial y}(y s_z) \right],
$$
we obtain the **exact pointwise identity**:
$$
\frac{d j_3^{\text{direct}}}{d\omega} = \frac{d j_3^{\text{orbital}}}{d\omega} + \frac{d j_3^{\text{spin}}}{d\omega} + \frac{d\mathcal{D}_z}{dz} + \nabla_\perp \cdot \mathbf{M}_\perp,
$$
where:

### 4.1 The $z$-Divergence Correction Term $\frac{d\mathcal{D}_z}{dz}$
$$
\frac{d\mathcal{D}_z}{dz} = \frac{1}{2} \frac{\partial}{\partial z}\left( x s_x + y s_y \right) = \frac{1}{2}\left[ x \frac{\partial s_x}{\partial z} + y \frac{\partial s_y}{\partial z} \right].
$$
The transverse components of spin density $(s_x, s_y)$ are:
$$
s_x = \frac{4\pi\varepsilon_0}{4\omega} \operatorname{Im}\left[ (E_2^* E_3 - E_3^* E_2) + c^2 (B_2^* B_3 - B_3^* B_2) \right],
$$
$$
s_y = \frac{4\pi\varepsilon_0}{4\omega} \operatorname{Im}\left[ (E_3^* E_1 - E_1^* E_3) + c^2 (B_3^* B_1 - B_1^* B_3) \right].
$$

### 4.2 The Transverse Boundary Divergence Term $\nabla_\perp \cdot \mathbf{M}_\perp$
$$
\nabla_\perp \cdot \mathbf{M}_\perp = -\frac{1}{2} \left[ \frac{\partial}{\partial x}(x s_z) + \frac{\partial}{\partial y}(y s_z) \right].
$$

---

## 5. Screen Surface Integration

When integrating across the 2D transverse screen $\iint_{\text{screen}} dx\, dy$:

1. By the 2D divergence theorem, the transverse term $\nabla_\perp \cdot \mathbf{M}_\perp$ becomes a boundary line integral along the edge of the detector screen:
   $$
   \iint_{\text{screen}} \nabla_\perp \cdot \mathbf{M}_\perp\, dx\,dy = -\frac{1}{2} \oint_{\partial \text{screen}} s_z (\mathbf{r} \cdot \mathbf{n}_{\text{edge}})\, d\ell.
   $$
   For a circular aperture of radius $R$:
   $$
   -\frac{R^2}{2} \int_0^{2\pi} s_z(R, \phi)\, d\phi.
   $$
   If the beam is sufficiently localized within the screen such that the fields vanish at the screen boundary, this contour integral is zero.

2. In contrast, the $z$-derivative term **does not vanish** upon integrating over $dx\,dy$:
   $$
   \iint_{\text{screen}} \frac{d j_3^{\text{direct}}}{d\omega}\, dx\,dy = \iint_{\text{screen}} \left( \frac{d j_3^{\text{spin}}}{d\omega} + \frac{d j_3^{\text{orbital}}}{d\omega} \right) dx\,dy + \frac{d}{dz} \iint_{\text{screen}} \frac{1}{2}\left( x s_x + y s_y \right) dx\,dy.
   $$

Therefore, to verify the exact equality on a detector screen at position $z$, you must compute the $z$-derivative of the transverse spin moment:
$$
\Delta J_3(z) = \frac{d}{dz} \iint_{\text{screen}} \frac{1}{2}\left[ x s_x(x, y, z) + y s_y(x, y, z) \right] dx\,dy.
$$
Evaluating the numerical derivative $\frac{\partial}{\partial z}$ requires fields evaluated at two neighboring planes $z \pm \Delta z$:
$$
\frac{\partial s_{x,y}}{\partial z} \approx \frac{s_{x,y}(z + \Delta z) - s_{x,y}(z - \Delta z)}{2\Delta z}.
$$

---

## 6. Implementation Workflow with Divergence Check

For each frequency $\omega > 0$:

1. **Evaluate Fields at $z - \Delta z$, $z$, and $z + \Delta z$:**
   Extract $(E_i, B_i)$ from $F^{\alpha\beta}$.

2. **Compute Field Derivatives at the Screen ($z$):**
   * Transverse: $\partial_\phi f = x \partial_y f - y \partial_x f$.
   * Longitudinal: $\partial_z s_x \approx \frac{s_x(z+\Delta z) - s_x(z-\Delta z)}{2\Delta z}$ and $\partial_z s_y \approx \frac{s_y(z+\Delta z) - s_y(z-\Delta z)}{2\Delta z}$.

3. **Compute Pointwise Quantities on Grid:**
   * $j_3^{\text{direct}}(x, y)$
   * $j_3^{\text{spin}}(x, y) = s_z(x, y)$
   * $j_3^{\text{orbital}}(x, y)$
   * $\frac{d\mathcal{D}_z}{dz}(x, y) = \frac{1}{2}\left( x \frac{\partial s_x}{\partial z} + y \frac{\partial s_y}{\partial z} \right)$
   * $\nabla_\perp \cdot \mathbf{M}_\perp(x, y) = -\frac{1}{2} \left[ \partial_x(x s_z) + \partial_y(y s_z) \right]$

4. **Pointwise Verification:**
   Verify at every pixel:
   $$
   j_3^{\text{direct}} - \left( j_3^{\text{spin}} + j_3^{\text{orbital}} + \frac{d\mathcal{D}_z}{dz} + \nabla_\perp \cdot \mathbf{M}_\perp \right) = 0 \quad (\pm\text{ numerical precision}).
   $$

5. **Integrated Screen Verification:**
   Verify over the screen:
   $$
   \iint j_3^{\text{direct}}\, dx\,dy = \iint \left( j_3^{\text{spin}} + j_3^{\text{orbital}} + \frac{d\mathcal{D}_z}{dz} \right) dx\,dy - \frac{1}{2}\oint s_z (\mathbf{r}\cdot\mathbf{n})\, d\ell.
   $$