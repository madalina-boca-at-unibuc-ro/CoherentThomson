# Angular Momentum Flux Density Along $Oz$

This document summarizes the local flux density (differential flux per unit area) for the $z$-component of electromagnetic angular momentum, in both the time domain and the spectral (Fourier) domain.

---

## 1. Time Domain

The linear momentum density of the electromagnetic field is:
$$\mathbf{g}(\mathbf{r}, t) = \varepsilon_0 \big( \mathbf{E}(\mathbf{r}, t) \times \mathbf{B}(\mathbf{r}, t) \big) = \frac{1}{c^2} \mathbf{S}(\mathbf{r}, t)$$

The angular momentum density vector with respect to the coordinate origin is:
$$\mathbf{j}(\mathbf{r}, t) = \mathbf{r} \times \mathbf{g}(\mathbf{r}, t) = \varepsilon_0 \Big( \mathbf{r} \times \big( \mathbf{E}(\mathbf{r}, t) \times \mathbf{B}(\mathbf{r}, t) \big) \Big)$$

The instantaneous local **angular momentum flux density** along the $z$-axis (rate of $J_z$ passing through an infinitesimal surface element $dx\,dy$ at position $(x, y, z_0)$) is given by:

$$\mathcal{F}_{J_z}(\mathbf{r}, t) = \left[ \mathbf{r} \times \mathbf{g}(\mathbf{r}, t) \right]_z = x g_y(\mathbf{r}, t) - y g_x(\mathbf{r}, t)$$

Expanding the cross products in terms of field components:

$$\mathcal{F}_{J_z}(\mathbf{r}, t) = \varepsilon_0 \Big[ x \big( E_z(\mathbf{r}, t) B_x(\mathbf{r}, t) - E_x(\mathbf{r}, t) B_z(\mathbf{r}, t) \big) - y \big( E_y(\mathbf{r}, t) B_z(\mathbf{r}, t) - E_z(\mathbf{r}, t) B_y(\mathbf{r}, t) \big) \Big]$$

---

## 2. Spectral / Fourier Domain

Defining the temporal Fourier transform according to:
$$\mathbf{E}(\mathbf{r}, t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{\mathbf{E}}(\mathbf{r}, \omega) e^{-i\omega t} d\omega, \qquad \mathbf{B}(\mathbf{r}, t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{\mathbf{B}}(\mathbf{r}, \omega) e^{-i\omega t} d\omega$$

By Parseval's theorem for cross-products, the time-integrated local angular momentum flux density is:
$$\int_{-\infty}^{\infty} \mathcal{F}_{J_z}(\mathbf{r}, t) \, dt = \int_0^\infty \frac{d\mathcal{F}_{J_z}(\mathbf{r}, \omega)}{d\omega} \, d\omega$$

where the **spectral angular momentum flux density per unit positive angular frequency $d\omega$** is:

$$\frac{d\mathcal{F}_{J_z}(\mathbf{r}, \omega)}{d\omega} = \frac{\varepsilon_0}{\pi} \operatorname{Re}\Big[ x \left( \tilde{E}_z(\mathbf{r}, \omega) \tilde{B}_x^*(\mathbf{r}, \omega) - \tilde{E}_x(\mathbf{r}, \omega) \tilde{B}_z^*(\mathbf{r}, \omega) \right) - y \left( \tilde{E}_y(\mathbf{r}, \omega) \tilde{B}_z^*(\mathbf{r}, \omega) - \tilde{E}_z(\mathbf{r}, \omega) \tilde{B}_y^*(\mathbf{r}, \omega) \right) \Big]$$

> **Note on Complex Conjugation:** 
> One of the Fourier transforms must be complex-conjugated (here $\tilde{\mathbf{B}}^*$). Equivalently, because $\operatorname{Re}[\tilde{\mathbf{E}} \times \tilde{\mathbf{B}}^*] = \operatorname{Re}[\tilde{\mathbf{E}}^* \times \tilde{\mathbf{B}}]$, conjugating $\tilde{\mathbf{E}}$ instead yields the exact same real physical quantity.


> **Note on the splitting in long/short components of the fields:**
Both the electric and magnetic field ar the sum between a long range component and a short range componenent, calculated distinctly in the numerical code. In principle the short range  component is essential for the angular momentum; we will test this behaviour. 
we calculate the contribution of the long lange component only 

$$\frac{d\mathcal{F}_{J_z}^{(ll)}(\mathbf{r}, \omega)}{d\omega} = \frac{\varepsilon_0}{\pi} \operatorname{Re}\Big[ x \left( \tilde{E}_{l,z}(\mathbf{r}, \omega) \tilde{B}_{l,x}^*(\mathbf{r}, \omega) - \tilde{E}_{l,x}(\mathbf{r}, \omega) \tilde{B}_{l,z}^*(\mathbf{r}, \omega) \right) - y \left( \tilde{E}_{l,y}(\mathbf{r}, \omega) \tilde{B}_{l,z}^*(\mathbf{r}, \omega) - \tilde{E}_{l,z}(\mathbf{r}, \omega) \tilde{B}_{l,y}^*(\mathbf{r}, \omega) \right) \Big]$$

and the two terms combining the long and short components
$$\frac{d\mathcal{F}_{J_z}^{(ls)}(\mathbf{r}, \omega)}{d\omega} = \frac{\varepsilon_0}{\pi} \operatorname{Re}\Big[ x \left( \tilde{E}_{l,z}(\mathbf{r}, \omega) \tilde{B}_{s,x}^*(\mathbf{r}, \omega) - \tilde{E}_{l,x}(\mathbf{r}, \omega) \tilde{B}_{s,z}^*(\mathbf{r}, \omega) \right) - y \left( \tilde{E}_{l,y}(\mathbf{r}, \omega) \tilde{B}_{s,z}^*(\mathbf{r}, \omega) - \tilde{E}_{l,z}(\mathbf{r}, \omega) \tilde{B}_{s,y}^*(\mathbf{r}, \omega) \right) \Big]$$

$$\frac{d\mathcal{F}_{J_z}^{(sl)}(\mathbf{r}, \omega)}{d\omega} = \frac{\varepsilon_0}{\pi} \operatorname{Re}\Big[ x \left( \tilde{E}_{s,z}(\mathbf{r}, \omega) \tilde{B}_{l,x}^*(\mathbf{r}, \omega) - \tilde{E}_{s,x}(\mathbf{r}, \omega) \tilde{B}_{l,z}^*(\mathbf{r}, \omega) \right) - y \left( \tilde{E}_{s,y}(\mathbf{r}, \omega) \tilde{B}_{l,z}^*(\mathbf{r}, \omega) - \tilde{E}_{s,z}(\mathbf{r}, \omega) \tilde{B}_{l,y}^*(\mathbf{r}, \omega) \right) \Big]$$
