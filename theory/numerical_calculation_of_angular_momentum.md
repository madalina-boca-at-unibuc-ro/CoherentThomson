
# Numerical calculation of the angular momentum density and flux

## Calculation of the temporal average of a bi-linear combinations of vector fields

If vector quantities are known in the Fourier transform, 

$${\bf A}(t) = \int\limits_{-\infty}^{\infty}d\omega e^{-i\omega t}\tilde{\bf A}(\omega)$$
$${\bf B}(t) = \int\limits_{-\infty}^{\infty}d\omega e^{-i\omega t}\tilde{\bf B}(\omega)$$
and we know that they are real ${\bf A}(t)={\bf A}^*(t)$, ${\bf B}(t) = {\bf B}^*(t)$, we can easily derive the consequence for the Fourier transforms as 
$$\tilde {\bf A}(\omega)={\bf A}^*(-\omega),\qquad \tilde{\bf B}(\omega) = {\bf B}^*(-\omega)$$ 
We want to express the temporal average of a bi-linear form of the type 

$$C_{ij}(t)=A_i(t)B_j(t)$$
$$\langle C_{ij}\rangle = \int\limits_{0}^T dt A_i(t) B_j(t)$$
by direct calculation we obtain 
$$\langle C_{ij}\rangle = \int\limits_{-\infty}^{\infty}dt\int\limits_{-\infty}^{\infty}d\omega\int\limits_{-\infty}^{\infty}d\omega'e^{-it(\omega-\omega')}\tilde A_i(\omega)\tilde B_j^*(\omega')$$
We can calculate the temporal integral which gives $2\pi\delta(\omega-\omega')$, so the final result is 
$$\langle C_{ij}\rangle = \int\limits_{-\infty}^{\infty}d\omega \tilde A_i(\omega)\tilde B_j^*(\omega)=\int\limits_{-\infty}^{\infty}d\omega \tilde A_i(\omega)\tilde B_j(-\omega)=2\int\limits_{0}^{\infty}d\omega \tilde A_i(\omega)\tilde B_j^*(\omega)$$
The Fourier spectrum of the total (time integrated) $\langle C_{ij}\rangle$ is then (for positive frequencies only):
$$\frac{d\langle C_{ij}\rangle}{d\omega} = 2\tilde A_i(\omega)\tilde B_j^*(\omega),\quad \omega\in(0,\infty)$$
Since $\langle C_{ij}\rangle$ is real, the previous result can be put into the form 
$$\frac{dC_{ij}}{d\omega}=\tilde A_i(\omega)\tilde B_j^*(\omega)+\tilde A_i^*(\omega)\tilde B_j(\omega)=2\operatorname{Re}(\tilde A_i(\omega)\tilde B_j^*(\omega))$$

## Applications for the angular momentum density

### Notations
**In the following we only discuss Fourier transforms of different quantities.**  
**We also consider now only the third component of the angular momentum.**  
**The argument $(\omega)$ will usually be omitted; Fourier components are denoted by a tilde. Spatial arguments will also be omitted; all quantities depend on the spatial coordinate ${\bf x}$.**

The basic quantity already calculated on a grid is the Faraday tensor:
$$F^{\mu\nu}({\bf x},\omega)=\left(\begin{array}{cccc}0&-\tilde E_x/c& -\tilde E_y/c& -\tilde E_z/c\\\tilde E_x/c&0&-\tilde B_z&\tilde B_y\\\tilde E_y/c&\tilde B_z&0&-\tilde B_x\\\tilde E_z/c&-\tilde B_y&\tilde B_x&0\end{array}\right)$$

From which we obtain the field components:
<a id="eq-definition-of-E-and-B"></a>
$$\tilde E_x = c F^{10},\quad \tilde E_y = c F^{20},\quad \tilde E_z = c F^{30},\quad \tilde B_x = F^{32},\quad \tilde B_y = F^{13},\quad \tilde B_z = F^{21} \tag{definition-of-E-and-B}$$

From the components we calculate the auxiliary vectors:
$$\tilde {\bf A} = \frac{1}{i\omega}\tilde {\bf E},\qquad \tilde {\bf C}=\frac{c}{i\omega}\tilde {\bf B}$$
which in terms of the tensor components are:
<a id="eq-definition-of-A"></a>
$$A_x = \frac{c}{i\omega}F^{10},\quad A_y = \frac{c}{i\omega}F^{20},\quad A_z = \frac{c}{i\omega}F^{30} \tag{definition-of-A}$$

<a id="eq-definition-of-C"></a>
$$C_x = \frac{c}{i\omega}F^{32},\quad C_y = \frac{c}{i\omega}F^{13},\quad C_z = \frac{c}{i\omega}F^{21} \tag{definition-of-C}$$

### 1. Primary definition of OAM and SAM ($Oz$ components) densities

The fundamental definition in the temporal domain is:
$${\boldsymbol{\cal S}} = {\bf E}\times{\bf A}$$

which gives for the Fourier spectrum:
$$\frac{d {\cal S}_z}{d\omega} = 2\operatorname{Re}\left[\tilde E_x \tilde A_y^* - \tilde E_y \tilde A_x^*\right] = -\frac{2}{\omega}\operatorname{Re}\left[\frac{1}{i}\tilde E_x \tilde E_y^* - \frac{1}{i}\tilde E_y \tilde E_x^*\right]$$

and the final result is:
<a id="eq-primary-definition-of-Sz"></a>
$$\frac{d {\cal S}_z}{d\omega} = -\frac{4}{\omega}\operatorname{Im}\left[\tilde E_x \tilde E_y^*\right] = \frac{4}{\omega}\operatorname{Im}\left[\tilde E_x^* \tilde E_y\right] \tag{primary-definition-of-Sz}$$

## Indications for numerical implementation in Python

1. From the file with the exported values of the Faraday tensor on the screen, read the values of $F^{\mu\nu}$.
2. At every point on the grid, calculate the Cartesian components of the vectors $\tilde{\bf E}$, $\tilde{\bf B}$, $\tilde{\bf A}$, and $\tilde{\bf C}$ using formulas [Eq. (definition-of-E-and-B)](#eq-definition-of-E-and-B), [Eq. (definition-of-A)](#eq-definition-of-A), and [Eq. (definition-of-C)](#eq-definition-of-C).
3. Calculate $\frac{d{\cal S}_z}{d\omega}$ according to [Eq. (primary-definition-of-Sz)](#eq-primary-definition-of-Sz).
4. Represent $\frac{d{\cal S}_z}{d\omega}$ in a heatmap (similar to the existing plots of the Faraday tensor components).
5. Implement the previous operations for both exported Faraday tensors (the emitted radiation and the incident LG beam).