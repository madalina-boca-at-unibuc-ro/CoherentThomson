# Implementing the long-distance direct and simplified Fourier formulas

This guide specifies an extension to a code that already evaluates the exact
Fourier-transformed fields. Preserve the existing exact calculations and add
long-distance direct and simplified calculations using the same trajectory
input, units, Fourier convention, and output representation.

The formulas follow
[FT-long_range_direct-simplified.tex](../../shared_src/elm_el_sc/frequency-domain/FT-long_range_direct-simplified.tex).
The two-vector method follows
[many_electron_far_field_two_projections.md](many_electron_far_field_two_projections.md),
expressed entirely in the notation of the LaTeX source. No abbreviations for
phases, contractions, or prefactors are introduced here.

## 1. Conventions and existing inputs

Use

$$
g_{\mu\nu}=\operatorname{diag}(1,-1,-1,-1),\qquad
\widetilde f(\omega)=\frac1{2\pi}\int dt\,e^{i\omega t}f(t),\qquad
k=\frac{\omega}{c},
$$

$$
E_i=cF^{i0},\qquad B_i=-\frac12\epsilon_{ijk}F^{jk}.
$$

The symbol $e$ is the signed source charge, consistent with the existing
exact implementation. Bold-vector dot products are Euclidean; four-vector
dot products use the metric above.

Reuse the laboratory trajectory and its proper-time derivatives:

$$
r_0^0(\tau)=ct'(\tau),\qquad
u^\mu=\frac{dr_0^\mu}{d\tau},\qquad
w^\mu=\frac{du^\mu}{d\tau}.
$$

Here $t'$ is laboratory emission time, $t$ is observation time, and $\tau$
is proper time. In particular, $\mathbf u$ is the spatial part of the
four-velocity, not the ordinary laboratory three-velocity; $\mathbf w$ is
the spatial part of the four-acceleration. Do not pass ordinary velocity
or acceleration into these formulas without the appropriate conversion.

For each trajectory and observation point, fix

$$
\boldsymbol{\mathfrak R}_0=\mathbf r(\tau_m),\qquad
\mathbf x_0=\mathbf x-\boldsymbol{\mathfrak R}_0,\qquad
\mathbf r_0(\tau)=\mathbf r(\tau)-\boldsymbol{\mathfrak R}_0,
$$

$$
\mathbf R_0(\tau)=\mathbf x_0-\mathbf r_0(\tau),\qquad
\mathbf n_0=\frac{\mathbf x_0}{|\mathbf x_0|},\qquad
n_0=(1,\mathbf n_0).
$$

Require $|\mathbf x_0|>0$. Keep this reference point and $\mathbf n_0$
fixed throughout the integral, including when splitting it into numerical
subintervals. The spatial translation does not reset laboratory time.

The contractions and transverse vectors are

$$
n_0\cdot u=u^0-\mathbf n_0\cdot\mathbf u,\qquad
n_0\cdot w=w^0-\mathbf n_0\cdot\mathbf w,
$$

$$
u_\perp=(0,\mathbf u-\mathbf n_0(\mathbf n_0\cdot\mathbf u)),\qquad
w_\perp=(0,\mathbf w-\mathbf n_0(\mathbf n_0\cdot\mathbf w)).
$$

Because $\mathbf n_0$ is constant, $du_\perp/d\tau=w_\perp$.
For a future-directed timelike trajectory, $n_0\cdot u>0$, although it
can be small for strongly beamed motion.

## 2. Approximation and integration interval

The long-distance calculation uses

$$
\frac1{|\mathbf R_0|}\simeq\frac1{|\mathbf x_0|},\qquad
n_{R_0}\simeq n_0,\qquad
|\mathbf R_0|\simeq|\mathbf x_0|-\mathbf n_0\cdot\mathbf r_0
\quad\text{in the phase}.
$$

It omits the time-domain velocity field proportional to
$1/|\mathbf R_0|^2$. Require a source segment small compared with
$|\mathbf x_0|$, and

$$
\frac{|k|\,|\mathbf r_0-\mathbf n_0(\mathbf n_0\cdot\mathbf r_0)|^2}
{2|\mathbf x_0|}\ll1
$$

throughout the contributing segment. For strongly beamed motion also
check that the direction replacement gives a small relative change in
$n_0\cdot u$; a small geometric angle alone may not suffice.

All formulas below use retarded proper-time limits $\tau_m,\tau_M$.
If the API specifies an observation-time window, determine its endpoints
using

$$
ct(\tau)\simeq r_0^0(\tau)+|\mathbf x_0|
-\mathbf n_0\cdot\mathbf r_0(\tau).
$$

Use this same approximate prescription for both long-distance methods.
The exact calculation uses the exact retarded relation. A comparison
must specify whether it fixes the observation window or the source
proper-time interval; these are different choices at finite distance.
For an observation-window API, choose the fixed spatial reference before
solving for its retarded endpoints; the choice
$\boldsymbol{\mathfrak R}_0=\mathbf r(\tau_m)$ is a convenient convention
when the source interval is already specified, not a reason to move the
reference during the endpoint solve.

These are proper-time integrals. If the existing quadrature uses $t'$,
apply $d\tau=c\,dt'/u^0$ to each integral. The observation-time Jacobian
$dt/d\tau=(n_0\cdot u)/c$ is already included in the formulas below;
do not multiply by it a second time.

## 3. Direct tensor formula

The direct long-distance result is

$$
\boxed{
\begin{aligned}
F_a^{\alpha\beta}(ck,\mathbf x)
={}&\frac1{2\pi}\frac{e}{4\pi\epsilon_0c^2}
\frac{e^{ik|\mathbf x_0|}}{|\mathbf x_0|}
\int_{\tau_m}^{\tau_M}d\tau\,
e^{ik(r_0^0-\mathbf n_0\cdot\mathbf r_0)}\\
&\times\frac{(n_0\cdot u)(n_0^\alpha w_\perp^\beta-n_0^\beta w_\perp^\alpha)
-(n_0\cdot w)(n_0^\alpha u_\perp^\beta-n_0^\beta u_\perp^\alpha)}
{(n_0\cdot u)^2}.
\end{aligned}
}
$$

Replacing $u_\perp,w_\perp$ by $u,w$ in this numerator gives the same
tensor. The longitudinal contributions cancel. This is useful for
comparison with an existing direct tensor implementation.

No separate endpoint term is added to the direct formula, even on a
finite interval. It needs reliable four-acceleration data. Prefer the
existing equations of motion or acceleration evaluator to numerical
differentiation of noisy velocity samples.

## 4. Simplified tensor formula

Integration by parts gives

$$
\boxed{F_a^{\alpha\beta}=F_l^{\alpha\beta}+F_b^{\alpha\beta},}
$$

$$
\boxed{
F_l^{\alpha\beta}(ck,\mathbf x)
=-ik\frac1{2\pi}\frac{e}{4\pi\epsilon_0c^2}
\frac{e^{ik|\mathbf x_0|}}{|\mathbf x_0|}
\int_{\tau_m}^{\tau_M}d\tau\,
e^{ik(r_0^0-\mathbf n_0\cdot\mathbf r_0)}
(n_0^\alpha u_\perp^\beta-n_0^\beta u_\perp^\alpha),
}
$$

$$
\boxed{
F_b^{\alpha\beta}(ck,\mathbf x)
=\frac1{2\pi}\frac{e}{4\pi\epsilon_0c^2}
\frac{e^{ik|\mathbf x_0|}}{|\mathbf x_0|}
\left[
e^{ik(r_0^0-\mathbf n_0\cdot\mathbf r_0)}
\frac{n_0^\alpha u_\perp^\beta-n_0^\beta u_\perp^\alpha}{n_0\cdot u}
\right]_{\tau_m}^{\tau_M}.
}
$$

Square brackets mean upper endpoint minus lower endpoint, including the
phase and every factor inside the brackets. Here $F_l$ denotes the bulk
term of this simplified long-distance formula. Do not equate it alone
with an exact acceleration-field contribution bearing a similar label
in the existing code. The physical long-distance result is $F_a$.

The simplified method needs four-velocity, but no four-acceleration.
Retain $F_b$ for finite windows. Its omission is justified only by an
explicit full-time endpoint/convergence prescription. The derivation uses

$$
\frac{d}{d\tau}e^{ik(r_0^0-\mathbf n_0\cdot\mathbf r_0)}
=ik(n_0\cdot u)e^{ik(r_0^0-\mathbf n_0\cdot\mathbf r_0)}.
$$

## 5. Two fixed transverse vectors for either method

Choose a right-handed orthonormal basis of the transverse plane:

$$
\mathbf e_a\cdot\mathbf n_0=0,\qquad
\mathbf e_a\cdot\mathbf e_b=\delta_{ab},\qquad
\mathbf e_1\times\mathbf e_2=\mathbf n_0,\qquad a,b=1,2,
$$

$$
e_a^\mu=(0,\mathbf e_a).
$$

The basis index $a=1,2$ in $e_a,p_a$ is a summation index; the subscript
on $F_a$ retains its existing meaning as the long-distance field.

For a robust construction, choose the Cartesian unit vector $\mathbf a$
with the smallest $|\mathbf a\cdot\mathbf n_0|$, then calculate

$$
\mathbf e_1=
\frac{\mathbf a-(\mathbf a\cdot\mathbf n_0)\mathbf n_0}
{|\mathbf a-(\mathbf a\cdot\mathbf n_0)\mathbf n_0|},\qquad
\mathbf e_2=\mathbf n_0\times\mathbf e_1.
$$

Calculate this basis once per trajectory/observation-point geometry.
Do not rotate it with the instantaneous velocity or retarded direction.

The decompositions are

$$
u_\perp^\mu=\sum_{a=1}^2(\mathbf e_a\cdot\mathbf u)e_a^\mu,\qquad
w_\perp^\mu=\sum_{a=1}^2(\mathbf e_a\cdot\mathbf w)e_a^\mu.
$$

Thus the full transverse vectors need not be built at every quadrature
sample. Project the original laboratory vectors directly onto the two
basis vectors. The trajectory remains in laboratory coordinates; no
Lorentz boost or new trajectory integration is needed.

### 5.1. Direct evaluation of the two amplitudes

For each $a=1,2$, calculate

$$
\boxed{
p_a=-\int_{\tau_m}^{\tau_M}d\tau\,
e^{ik(r_0^0-\mathbf n_0\cdot\mathbf r_0)}
\frac{(n_0\cdot u)(\mathbf e_a\cdot\mathbf w)
-(n_0\cdot w)(\mathbf e_a\cdot\mathbf u)}{(n_0\cdot u)^2}.
}
$$

At each quadrature sample evaluate the phase, $n_0\cdot u$, $n_0\cdot w$,
and the two projections of each spatial vector. Integrate the two complex
components together using the same nodes. There is no additional boundary
correction. Keep the leading minus sign: the reconstruction below uses
$e_a^\alpha n_0^\beta-n_0^\alpha e_a^\beta$.

### 5.2. Simplified evaluation of the same two amplitudes

For each $a=1,2$, calculate

$$
\boxed{
p_a=ik\int_{\tau_m}^{\tau_M}d\tau\,
(\mathbf e_a\cdot\mathbf u)e^{ik(r_0^0-\mathbf n_0\cdot\mathbf r_0)}
-\left[
\frac{\mathbf e_a\cdot\mathbf u}{n_0\cdot u}
e^{ik(r_0^0-\mathbf n_0\cdot\mathbf r_0)}
\right]_{\tau_m}^{\tau_M}.
}
$$

Integrate the two projected velocities with the phase, multiply by $ik$,
and subtract the endpoint difference. Evaluate endpoint values at the
actual limits, using the trajectory interpolator if necessary; do not
silently replace them with nearby stored samples.

The integral part reconstructs $F_l$, and the negative endpoint part
reconstructs $F_b$, using the reconstruction below. The positive tensor
boundary sign in section 4 and negative scalar boundary sign here agree
because their antisymmetric basis factors have opposite ordering.

Both methods compute the same $p_a$. Indeed,

$$
\frac{d}{d\tau}\left(\frac{\mathbf e_a\cdot\mathbf u}{n_0\cdot u}\right)
=\frac{(n_0\cdot u)(\mathbf e_a\cdot\mathbf w)
-(n_0\cdot w)(\mathbf e_a\cdot\mathbf u)}{(n_0\cdot u)^2}.
$$

Integration by parts in the direct scalar formula therefore yields the
simplified scalar formula exactly within the same long-distance geometry.
The direct method can avoid cancellation between large bulk and endpoint
terms; the simplified method avoids requiring acceleration. Compare them
with adequately resolved quadrature before interpreting a discrepancy as
a physical approximation error.

## 6. Common reconstruction for both methods

After either calculation of $p_1,p_2$, reconstruct

$$
\boxed{
F_a^{\alpha\beta}(ck,\mathbf x)
=\frac1{2\pi}\frac{e}{4\pi\epsilon_0c^2}
\frac{e^{ik|\mathbf x_0|}}{|\mathbf x_0|}
\sum_{a=1}^2p_a(e_a^\alpha n_0^\beta-n_0^\alpha e_a^\beta).
}
$$
