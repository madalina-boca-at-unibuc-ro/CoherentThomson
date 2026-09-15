# Numerical calculation of the angular momentum density and flux

## Calculation of the temporal average of a bi-linear combination of vector fields

If real vector quantities ${\bf A}(t)={\bf A}^*(t)$ and ${\bf B}(t) = {\bf B}^*(t)$ are represented by their Fourier transforms:

$${\bf A}(t) = \int\limits_{-\infty}^{\infty}d\omega\, e^{-i\omega t}\tilde{\bf A}(\omega),\qquad {\bf B}(t) = \int\limits_{-\infty}^{\infty}d\omega\, e^{-i\omega t}\tilde{\bf B}(\omega)$$

their frequency components satisfy:
$$\tilde{\bf A}(\omega)={\bf A}^*(-\omega),\qquad \tilde{\bf B}(\omega) = {\bf B}^*(-\omega)$$

We express the time-integrated average of the bilinear form $C_{ij}(t) = A_i(t)B_j(t)$:
$$\langle C_{ij}\rangle = \int\limits_{-\infty}^{\infty}dt\, A_i(t) B_j(t) = \int\limits_{-\infty}^{\infty}dt\int\limits_{-\infty}^{\infty}d\omega\int\limits_{-\infty}^{\infty}d\omega' e^{-it(\omega-\omega')}\tilde A_i(\omega)\tilde B_j^*(\omega')$$

Evaluating the temporal integral gives $2\pi\delta(\omega-\omega')$, yielding:
$$\langle C_{ij}\rangle = \int\limits_{-\infty}^{\infty}d\omega\, \tilde A_i(\omega)\tilde B_j^*(\omega) = 2\int\limits_{0}^{\infty}d\omega\, \operatorname{Re}\left[\tilde A_i(\omega)\tilde B_j^*(\omega)\right]$$

The one-sided Fourier spectral density for $\omega \in (0,\infty)$ is therefore:
<a id="eq-bilinear-spectral-density"></a>
$$\frac{d\langle C_{ij}\rangle}{d\omega} = 2\operatorname{Re}\left[\tilde A_i(\omega)\tilde B_j^*(\omega)\right] = \tilde A_i(\omega)\tilde B_j^*(\omega) + \tilde A_i^*(\omega)\tilde B_j(\omega) \tag{bilinear-spectral-density}$$

---

## Applications for the angular momentum density

### Notations and Conventions
* All field quantities below are evaluated in the frequency domain for a positive frequency $\omega > 0$. The explicit argument $(\omega)$ and spatial coordinate ${\bf x}$ are suppressed for brevity.
* Fourier amplitudes are denoted with a tilde ($\tilde{\bf E}, \tilde{\bf B}, \tilde{\bf A}$).
* We consider the third ($z$) component of the angular momentum along the beam propagation / detector normal axis.

The complex Faraday tensor calculated on the detector grid is:
$$F^{\mu\nu}({\bf x},\omega)=\left(\begin{array}{cccc}0&-\tilde E_x/c& -\tilde E_y/c& -\tilde E_z/c\\\tilde E_x/c&0&-\tilde B_z&\tilde B_y\\\tilde E_y/c&\tilde B_z&0&-\tilde B_x\\\tilde E_z/c&-\tilde B_y&\tilde B_x&0\end{array}\right)$$

From which we extract the electric and magnetic field components:
<a id="eq-definition-of-E-and-B"></a>
$$\tilde E_x = c F^{10},\quad \tilde E_y = c F^{20},\quad \tilde E_z = c F^{30},\quad \tilde B_x = F^{32},\quad \tilde B_y = F^{13},\quad \tilde B_z = F^{21} \tag{definition-of-E-and-B}$$

In the temporal gauge / Coulomb-radiation gauge, the vector potential $\tilde{\bf A}$ and auxiliary vector $\tilde{\bf C}$ satisfy:
$$\tilde {\bf A} = \frac{1}{i\omega}\tilde {\bf E},\qquad \tilde {\bf C}=\frac{c}{i\omega}\tilde {\bf B}$$

In terms of the tensor components:
<a id="eq-definition-of-A"></a>
$$\tilde A_x = \frac{c}{i\omega}F^{10},\quad \tilde A_y = \frac{c}{i\omega}F^{20},\quad \tilde A_z = \frac{c}{i\omega}F^{30} \tag{definition-of-A}$$

<a id="eq-definition-of-C"></a>
$$\tilde C_x = \frac{c}{i\omega}F^{32},\quad \tilde C_y = \frac{c}{i\omega}F^{13},\quad \tilde C_z = \frac{c}{i\omega}F^{21} \tag{definition-of-C}$$

---

### 1. Spin Angular Momentum (SAM) Density ($Oz$ component)

The fundamental definition of the spin angular momentum density vector in the time domain is:
$${\boldsymbol{\cal S}} = \epsilon_0{\bf E}\times{\bf A}$$

Its $z$-component is ${\cal S}_z = \epsilon_0(E_x A_y - E_y A_x)$. Applying [Eq. (bilinear-spectral-density)](#eq-bilinear-spectral-density) gives:
$$\frac{d {\cal S}_z}{d\omega} = 2\epsilon_0\operatorname{Re}\left[\tilde E_x \tilde A_y^* - \tilde E_y \tilde A_x^*\right]$$

Substituting $\tilde A_i^* = \frac{i}{\omega}\tilde E_i^*$:
$$\tilde E_x \tilde A_y^* - \tilde E_y \tilde A_x^* = \frac{i}{\omega}\left(\tilde E_x \tilde E_y^* - \tilde E_y \tilde E_x^*\right) = -\frac{2}{\omega}\operatorname{Im}\left[\tilde E_x \tilde E_y^*\right] = \frac{2}{\omega}\operatorname{Im}\left[\tilde E_x^* \tilde E_y\right]$$

Thus, the SAM spectral density is:
<a id="eq-primary-definition-of-Sz"></a>
$$\frac{d {\cal S}_z}{d\omega} = \frac{4\epsilon_0}{\omega}\operatorname{Im}\left[\tilde E_x^* \tilde E_y\right] = -\frac{4\epsilon_0}{\omega}\operatorname{Im}\left[\tilde E_x \tilde E_y^*\right] \tag{primary-definition-of-Sz}$$

---

### 2. Orbital Angular Momentum (OAM) Density ($Oz$ component)

The canonical orbital angular momentum density in the time domain is:
$${\boldsymbol{\cal L}} = \epsilon_0 \sum_{i=x,y,z} E_i\, ({\bf r}\times{\boldsymbol\nabla})\, A_i$$

Its $z$-component involves the transverse differential operator $\hat{L}_z = ({\bf r}\times{\boldsymbol\nabla})_z = (x\partial_y - y\partial_x)$:
$${\cal L}_z = \epsilon_0 \sum_{i=x,y,z} E_i\, (x\partial_y - y\partial_x) A_i$$

Applying [Eq. (bilinear-spectral-density)](#eq-bilinear-spectral-density):
$$\frac{d {\cal L}_z}{d\omega} = 2\epsilon_0 \sum_{i=x,y,z} \operatorname{Re}\left[\tilde E_i\, \Big(\hat{L}_z \tilde A_i\Big)^*\right]$$

Using $\tilde A_i^* = \frac{i}{\omega}\tilde E_i^*$ and noting that coordinates and spatial derivative operators are real:
$$\Big(\hat{L}_z \tilde A_i\Big)^* = \frac{i}{\omega}\hat{L}_z \tilde E_i^*$$
$$\operatorname{Re}\left[\tilde E_i \left(\frac{i}{\omega}\hat{L}_z \tilde E_i^*\right)\right] = -\frac{1}{\omega}\operatorname{Im}\left[\tilde E_i \hat{L}_z \tilde E_i^*\right] = \frac{1}{\omega}\operatorname{Im}\left[\tilde E_i^*\, \hat{L}_z \tilde E_i\right]$$

#### A. Rectangular Screen (Cartesian Grid)
On a Cartesian grid $(x, y)$, the operator is $\hat{L}_z = x \partial_y - y \partial_x$:
<a id="eq-primary-definition-of-Lz-rectangular"></a>
$$\frac{d {\cal L}_z}{d\omega} = \frac{2\epsilon_0}{\omega}\sum_{i=x,y,z}\operatorname{Im}\left[\tilde E_i^* \left(x \frac{\partial \tilde E_i}{\partial y} - y \frac{\partial \tilde E_i}{\partial x}\right)\right] \tag{primary-definition-of-Lz-rectangular}$$

#### B. Circular Screen (Polar Grid)
On a polar grid $(\rho, \phi)$ where $x = \rho\cos\phi$, $y = \rho\sin\phi$, we use the identity $x\partial_y - y\partial_x = \partial_\phi$:
<a id="eq-primary-definition-of-Lz-circular"></a>
$$\frac{d {\cal L}_z}{d\omega} = \frac{2\epsilon_0}{\omega}\sum_{i=x,y,z}\operatorname{Im}\left[\tilde E_i^* \frac{\partial \tilde E_i}{\partial \phi}\right] \tag{primary-definition-of-Lz-circular}$$

---

## Indications for numerical implementation in Python

1. **Target Screens:** Implement the angular momentum density calculations exclusively for the rectangular screen (Cartesian grid) and the circular screen (polar grid).
2. **Data Ingestion:** For each run and for both exported field states (the emitted radiation and the incident LG beam), load the exported Faraday tensor array $F^{\mu\nu}({\bf x},\omega)$ from disk.
3. **Field Reconstruction:** At every spatial grid node, reconstruct the complex electric field vector $\tilde{\bf E} = (\tilde E_x, \tilde E_y, \tilde E_z)$ and magnetic field vector $\tilde{\bf B}$ using [Eq. (definition-of-E-and-B)](#eq-definition-of-E-and-B).
4. **SAM Density Calculation:**
   - Compute the local SAM spectral density $\frac{d{\cal S}_z}{d\omega}$ using [Eq. (primary-definition-of-Sz)](#eq-primary-definition-of-Sz).
5. **OAM Density Calculation:**
   - **Rectangular Screen:** Compute numerical spatial derivatives $\frac{\partial \tilde E_i}{\partial x}$ and $\frac{\partial \tilde E_i}{\partial y}$ (using central finite differences `numpy.gradient` along axis 1 and axis 0, taking grid spacings $\Delta x, \Delta y$ into account). Evaluate $\frac{d{\cal L}_z}{d\omega}$ using [Eq. (primary-definition-of-Lz-rectangular)](#eq-primary-definition-of-Lz-rectangular).
   - **Circular Screen:** Compute the azimuthal derivative $\frac{\partial \tilde E_i}{\partial \phi}$ along the periodic angular coordinate grid $\phi$ (with central differences and periodic boundary conditions `numpy.gradient(..., axis=phi_axis)`). Evaluate $\frac{d{\cal L}_z}{d\omega}$ using [Eq. (primary-definition-of-Lz-circular)](#eq-primary-definition-of-Lz-circular).
6. **Total Angular Momentum (TAM) Density:**
   - Compute the local total angular momentum spectral density as the direct sum:
     $$\frac{d{\cal J}_z}{d\omega} = \frac{d{\cal L}_z}{d\omega} + \frac{d{\cal S}_z}{d\omega}$$
7. **Surface Integration (Integrated Flux / Power):**
   - Integrate each spectral density ($\frac{d{\cal S}_z}{d\omega}$, $\frac{d{\cal L}_z}{d\omega}$, and $\frac{d{\cal J}_z}{d\omega}$) across the screen using the composite 2D trapezoidal rule:
     - **Rectangular screen:** $\int\int \dots\, dx\, dy$ via `scipy.integrate.trapezoid` or nested `numpy.trapz`.
     - **Circular screen:** $\int\int \dots\, \rho\, d\rho\, d\phi$, explicitly including the radial Jacobian weight $\rho$.
   - Append the integrated values to the `run_log` file with timestamps, run parameters, and clear column labels.
8. **Visualization:**
   - Generate heatmaps for $\frac{d{\cal S}_z}{d\omega}$, $\frac{d{\cal L}_z}{d\omega}$, and $\frac{d{\cal J}_z}{d\omega}$.
   - Format colorscales, colormaps, aspect ratios, and spatial extent ticks consistently with the existing Faraday tensor component heatmaps.