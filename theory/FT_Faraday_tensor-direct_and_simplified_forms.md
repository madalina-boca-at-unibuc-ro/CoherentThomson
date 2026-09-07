# Fourier transform of the Faraday tensor — two equivalent forms

Source: `shared_src/elm_el_sc/frequency-domain/FT-direct.tex` (§ "direct
formula") and `shared_src/elm_el_sc/frequency-domain/FT-simplified.tex`
(§ "simplified form"). Each ends in a highlighted result box; both boxes are
transcribed below.

**Purpose of this note**: these are two of a planned set of *equivalent*
analytic expressions for the Fourier transform of the radiation field
$F^{\alpha\beta}(\omega,{\bf x})$ emitted by a single accelerated charge in
Thomson scattering. The end goal is a numerical cross-check: implement each
form independently and verify they produce the same $F^{\alpha\beta}(\omega,{\bf x})$
for the same input trajectory.

## Common notation (shared by both forms)

- $\tau$: proper time of the emitting charge, the variable of integration.
- ${\bf r}_0(\tau)$: particle trajectory (3-vector); $r_0^0(\tau)=c t(\tau)$ is
  the corresponding lab-time coordinate along the worldline.
- $u^\mu(\tau) = dr_0^\mu/d\tau$: four-velocity.
- $w^\mu(\tau) = du^\mu/d\tau$: four-acceleration.
- ${\bf x}_0$: fixed observation point (3-vector); $x_0^0 = ct$ (observation lab time, the outer Fourier-conjugate variable before the change of variable to $\tau$).
- ${\bf R}_0(\tau) = {\bf x}_0 - {\bf r}_0(\tau)$, $|{\bf R}_0| = |{\bf x}_0-{\bf r}_0(\tau)|$.
- ${\bf n}_{R_0}(\tau) = {\bf R}_0(\tau)/|{\bf R}_0(\tau)|$ (unit 3-vector).
- $n_{R_0} = (1, {\bf n}_{R_0})$: null four-vector, $R_0 = |{\bf R}_0|\, n_{R_0}$.
- $k = \omega/c$.
- Dot products are Minkowski four-products, e.g. $u\cdot n_{R_0} = u^0 - {\bf u}\cdot{\bf n}_{R_0}$ (metric signature $(+,-,-,-)$, consistent with $u\cdot u = c^2$).
- The phase in every integral is $k\,(r_0^0(\tau) + |{\bf R}_0(\tau)|)$, i.e. $\omega$ times the retarded-time relation $ct = r_0^0(\tau)+|{\bf x}_0-{\bf r}_0(\tau)|$ (the Jacobian of the change of variable $t\to\tau$ is $d(ct)/d\tau = n_{R_0}\cdot u$, already absorbed into the derivation, not left as a prefactor).
- Both forms decompose the result as $F^{\alpha\beta} = F_l^{\alpha\beta} + F_s^{\alpha\beta}$ ("long"/short-range split), but **the split is defined differently in each derivation** — see the note at the end. Only the **total** $F^{\alpha\beta}$ is guaranteed to agree between the two forms; do not expect $F_l$ (direct) to equal $F_l$ (simplified) term by term.

All prefactors use $e$ = charge, $\epsilon_0$ = vacuum permittivity, $c$ = speed of light.

---

## Form 1 — "direct" (from `FT-direct.tex`)

Starting point: direct Fourier transform of the closed-form primary expression
for $F^{\alpha\beta}(x)$ (Liénard–Wiechert field written without an explicit
$\tau$-derivative).

$$
F^{\alpha\beta}(ck,{\bf x}) = F^{\alpha\beta}_{l}(ck,{\bf x}) + F^{\alpha\beta}_{s}(ck,{\bf x})
$$

$$
F_l^{\alpha\beta}(ck,{\bf x}) =\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c}
\int_{-\infty}^{\infty} d\tau\; e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}\;
\frac{(u\cdot n_0)(n_0^{\alpha}w^{\beta}-n_0^{\beta}w^{\alpha})-(w\cdot n_0)(n_0^{\alpha}u^{\beta}-n_0^{\beta}u^{\alpha})}
{|{\bf R}_0|\,(u\cdot n_0)^3}
$$

$$
F_s^{\alpha\beta}(ck,{\bf x}) =\frac{1}{2\pi}\frac{e\,c^2}{4\pi\epsilon_0c}
\int_{-\infty}^{\infty} d\tau\; e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}\;
\frac{n_0^{\alpha}u^{\beta}-n_0^{\beta}u^{\alpha}}
{|{\bf R}_0|^2\,(u\cdot n_0)^3}
$$

(here $n_0 \equiv n_{R_0}(\tau)$, same object as above; the source file uses the
shorthand $n_0$ inside the highlightbox).

Notes for implementation:
- $F_l$ scales as $1/|{\bf R}_0|$ (radiation/far-field term); $F_s$ scales as
  $1/|{\bf R}_0|^2$ (near-field term) and carries an extra factor $c^2$.
- Both integrands are antisymmetric in $\alpha\leftrightarrow\beta$ by
  construction (each is built from $n_0^{\alpha}(\cdot)^\beta - n_0^\beta(\cdot)^\alpha$-type combinations), so only the 6 independent components of the antisymmetric tensor need to be evaluated.
- Denominator power: $(u\cdot n_0)^3$ in both terms.

---

## Form 2 — "simplified" (from `FT-simplified.tex`)

Starting point: Fourier transform of Jackson's compact form
$F^{\alpha\beta}(x) = \frac{e}{4\pi\epsilon_0 c}\frac{1}{u\cdot R_0}\frac{d}{d\tau}\!\left[\frac{R_0^\alpha u^\beta - R_0^\beta u^\alpha}{R_0\cdot u}\right]_{\tau=\tau_r}$,
followed by an integration by parts in $\tau$ that moves a derivative onto the
phase factor (boundary terms dropped; they only contribute as $\omega\to0$).

$$
F^{\alpha\beta}(ck,{\bf x}) = F^{\alpha\beta}_{l}(ck,{\bf x}) + F^{\alpha\beta}_{s}(ck,{\bf x})
$$

$$
F^{\alpha\beta}_{l}(ck,{\bf x}) =
\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c}\int_{-\infty}^{\infty} d\tau\;
e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}
\left(-ik\,\frac{u\cdot n_{R_0}}{|{\bf R}_0|}\right)
\frac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}{n_{R_0}\cdot u}
$$

$$
F^{\alpha\beta}_{s}(ck,{\bf x}) =
\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c}\int_{-\infty}^{\infty} d\tau\;
e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}
\left(\frac{{\bf n}_{R_0}\cdot{\bf u}}{|{\bf R}_0|^2}\right)
\frac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}{n_{R_0}\cdot u}
$$

Notes for implementation:
- Both terms share the same tensor factor $\dfrac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}{n_{R_0}\cdot u}$; they differ only in the scalar prefactor multiplying it: $-ik\,(u\cdot n_{R_0})/|{\bf R}_0|$ for $F_l$ (linear in $k$, radiation term) vs. $({\bf n}_{R_0}\cdot{\bf u})/|{\bf R}_0|^2$ for $F_s$ ($k$-independent, near-field term). Note ${\bf n}_{R_0}\cdot{\bf u}$ in $F_s$ is a **3-vector** dot product (spatial part of $u$), unlike the four-dot $u\cdot n_{R_0}$ used in $F_l$.
- This form is algebraically simpler to implement (one shared tensor factor, two scalar weights) and is a good candidate for a first numerical implementation.

---

## Shared notation block (reproduced from both source highlightboxes)

$$
k = \frac{\omega}{c},\qquad
R_0 = (|{\bf R}_0|,{\bf R}_0) = |{\bf R}_0|(1,{\bf n}_{R_0}) = |{\bf R}_0|\,n_{R_0},\qquad
{\bf R}_0 = {\bf x}_0-{\bf r}_0(\tau),\qquad
{\bf n}_{R_0} = \frac{{\bf R}_0}{|{\bf R}_0|}.
$$

In all integrals, every quantity built from $u$, $w$, ${\bf r}_0$, ${\bf R}_0$, $n_{R_0}$ is evaluated **as a function of the integration variable $\tau$** (not at a fixed retarded time) — the retarded-time constraint has already been absorbed into the change of variable $t\to\tau$ and appears only through the phase $e^{ik(r_0^0(\tau)+|{\bf R}_0(\tau)|)}$.

## Cross-check target

For a given trajectory ${\bf r}_0(\tau)$ (hence $u,w$ known analytically or
numerically) and observation point ${\bf x}_0$, both forms should be
numerically integrated over $\tau$ to give the same total
$F^{\alpha\beta}(ck,{\bf x})$ for every $\alpha,\beta,k$. Recommended check:
compare the **sum** $F_l+F_s$ from each form (not the individual $F_l$, $F_s$
pieces, which are not required to agree between the two derivations — see the
note on differing decompositions above).
