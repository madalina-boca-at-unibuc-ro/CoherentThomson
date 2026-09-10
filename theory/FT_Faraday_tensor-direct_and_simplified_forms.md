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
- The phase in every integral is $k\,(r_0^0(\tau) + |{\bf R}_0(\tau)|)$, i.e. $\omega$ times the retarded-time relation $ct = r_0^0(\tau)+|{\bf x}_0-{\bf r}_0(\tau)|$.
- **Jacobian of the change of variable $t\to\tau$.** The retarded-time relation gives $d(ct)/d\tau = u\cdot n_{R_0}$ directly, but the outer integral is over $t$, not $ct$, so the actual (dimensionless) Jacobian is $dt/d\tau = (u\cdot n_{R_0})/c$. This carries an explicit extra factor of $1/c$ relative to $d(ct)/d\tau$, and it must appear as an explicit multiplicative factor in the integrand once the variable is changed — it is not something that cancels or gets absorbed invisibly. Both source derivations originally mishandled this (see the bug note under Form 1 and Form 2 below); both are now fixed and consistently carry an overall $e/(4\pi\epsilon_0 c^2)$.
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
F_l^{\alpha\beta}(ck,{\bf x}) =\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c^2}
\int_{-\infty}^{\infty} d\tau\; \frac{1}{|{\bf R}_0(\tau)|}\,e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}\;
\frac{(u\cdot n_{R_0})(n_{R_0}^{\alpha}w^{\beta}-n_{R_0}^{\beta}w^{\alpha})-(w\cdot n_{R_0})(n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha})}
{(u\cdot n_{R_0})^2}
$$

$$
F_s^{\alpha\beta}(ck,{\bf x}) =\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0}
\int_{-\infty}^{\infty} d\tau\; \frac{1}{|{\bf R}_0(\tau)|^2}\,e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}\;
\frac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}
{(u\cdot n_{R_0})^2}
$$

Notes for implementation:
- $F_l$ scales as $1/|{\bf R}_0|$ (radiation/far-field term); $F_s$ scales as
  $1/|{\bf R}_0|^2$ (near-field term). Neither carries an explicit factor of
  $c^2$ any more — see the Jacobian bug note below.
- Both integrands are antisymmetric in $\alpha\leftrightarrow\beta$ by
  construction (each is built from $n_{R_0}^{\alpha}(\cdot)^\beta - n_{R_0}^\beta(\cdot)^\alpha$-type combinations), so only the 6 independent components of the antisymmetric tensor need to be evaluated.
- Denominator power: $(u\cdot n_{R_0})^2$ in both terms (**not** $^3$ — see below).
- **Do not confuse $n_{R_0}(\tau)$ with $n_0$.** $n_{R_0}(\tau)=(1,{\bf n}_{R_0}(\tau))$ is the *exact*, $\tau$-dependent unit direction from the emitting charge to the observation point, defined above and used throughout this note. $n_0=(1,{\bf n}_0)$, ${\bf n}_0={\bf x}_0/|{\bf x}_0|$, is a *different*, constant vector used only in the separate long-distance/far-field expansion (`shared_src/elm_el_sc/time-domain/c2.tex`, `c3.tex`); it never appears in the exact Fourier-transform formulas here. An earlier version of the `FT-direct.tex` source used `n_0` as shorthand for `n_{R_0}` in this section — that was a notation bug in the LaTeX source and has since been fixed to write `n_{R_0}` explicitly; this note now matches the corrected source.
- **Jacobian bug (fixed).** An earlier version of `FT-direct.tex` derived $d(ct)/d\tau=u\cdot n_{R_0}$ but then never multiplied the integrand by the corresponding $dt/d\tau=(u\cdot n_{R_0})/c$ when switching the outer integral from $t$ to $\tau$ — the Jacobian was dropped entirely. The source has been corrected to insert it explicitly, which (after using $u\cdot n_{R_0}=(u\cdot R_0)/|{\bf R}_0|$ to cancel one power of $u\cdot R_0$) changes the overall prefactor from $e/(4\pi\epsilon_0c)$ to $e/(4\pi\epsilon_0c^2)$ and drops the denominator power from $(u\cdot n_{R_0})^3$ to $(u\cdot n_{R_0})^2$ in both $F_l$ and $F_s$; $F_s$ also loses its separate $c^2$ prefactor (it cancels against the new $1/c^2$). The formulas above already reflect the corrected source.
- **Typesetting fix (this doc only, not a source/physics bug).** $|{\bf R}_0|$ and $|{\bf R}_0|^2$ were previously
  written in the constant prefactor in front of the $\int d\tau$, as if independent of the integration variable —
  but $|{\bf R}_0|=|{\bf R}_0(\tau)|$ is $\tau$-dependent (the electron-to-screen distance at that trajectory
  point), so it cannot be pulled outside the integral. Moved inside, multiplying the integrand, matching how Form
  2 below already has it. `Radiation::compute_radiation` (`radiation.cpp`) was never affected by this — its
  direct-form prefactors (`long_range_prefactor_direct`/`short_range_prefactor_direct`) are computed once per
  trajectory point inside the `i_tau` loop, using that point's own `R`, i.e. the code already treated $|{\bf
  R}_0|$ as $\tau$-dependent; only this `.md` transcription had it factored out.

---

## Form 2 — "simplified" (from `FT-simplified.tex`)

Starting point: Fourier transform of Jackson's compact form
$F^{\alpha\beta}(x) = \frac{e}{4\pi\epsilon_0 c}\frac{1}{u\cdot R_0}\frac{d}{d\tau}\!\left[\frac{R_0^\alpha u^\beta - R_0^\beta u^\alpha}{R_0\cdot u}\right]_{\tau=\tau_r}$,
followed by an integration by parts in $\tau$ that moves a derivative onto the
phase factor.

**Boundary term — only negligible for an infinite integration range.** The
original derivation drops the boundary term produced by this integration by
parts, on the grounds that it only contributes as $\omega\to0$ — true *only*
when $\tau$ ranges over all of $(-\infty,\infty)$ and the bracketed quantity
decays at both ends. `compute_radiation`'s actual $\tau$-integral is over a
**finite** window (`tau_0` to `tau_max`, i.e. one electron's finite recorded
trajectory), so this justification does not apply, and the boundary term
$F_b^{\alpha\beta}$ below must in general be kept. This was confirmed in a
related derivation (see "On-axis $F^{03}$ cancellation" below): dropping
$F_b$ leaves $F_l^{03}$ carrying a spurious, individually large artifact of
the integration by parts that does not appear in Form 1's direct
(non-integrated-by-parts) derivation — a plausible root cause of the
single-electron Thomson-dipole mismatch flagged as an open validation gap in
`CLAUDE.md`.

$$
F^{\alpha\beta}(ck,{\bf x}) = F^{\alpha\beta}_{l}(ck,{\bf x}) + F^{\alpha\beta}_{s}(ck,{\bf x})
$$

$$
F^{\alpha\beta}_{l}(ck,{\bf x}) =
\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c^2}\int_{-\infty}^{\infty} d\tau\;
e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}
\left(-ik\,\frac{u\cdot n_{R_0}}{|{\bf R}_0|}\right)
\frac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}{n_{R_0}\cdot u}
$$

$$
F^{\alpha\beta}_{s}(ck,{\bf x}) =
-\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c^2}\int_{-\infty}^{\infty} d\tau\;
e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}
\left(\frac{{\bf n}_{R_0}\cdot{\bf u}}{|{\bf R}_0|^2}\right)
\frac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}{n_{R_0}\cdot u}
$$

Notes for implementation:
- Both terms share the same tensor factor $\dfrac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}{n_{R_0}\cdot u}$; they differ only in the scalar weight multiplying it: $-ik\,(u\cdot n_{R_0})/|{\bf R}_0|$ for $F_l$ (linear in $k$, radiation term) vs. $-({\bf n}_{R_0}\cdot{\bf u})/|{\bf R}_0|^2$ for $F_s$ ($k$-independent, near-field term). Note ${\bf n}_{R_0}\cdot{\bf u}$ in $F_s$ is a **3-vector** dot product (spatial part of $u$), unlike the four-dot $u\cdot n_{R_0}$ used in $F_l$.
- This form is algebraically simpler to implement (one shared tensor factor, two scalar weights) and is a good candidate for a first numerical implementation.
- **Jacobian bug (fixed).** An earlier version of `FT-simplified.tex` combined the Jacobian $dt/d\tau=(u\cdot n_{R_0})/c$ with Jackson's own $1/(u\cdot R_0)$ prefactor into $\dfrac{u\cdot n_{R_0}}{u\cdot R_0}=\dfrac1{|{\bf R}_0|}$ — dropping the $1/c$, so the combined factor had units of $1/\text{length}$ instead of being dimensionless as required. The correct combination is $\dfrac{u\cdot n_{R_0}}{c\,(u\cdot R_0)}=\dfrac1{c\,|{\bf R}_0|}$. The source has been corrected: every prefactor from this point on carries $c^2$ instead of $c$ (the shape of the integration-by-parts result is otherwise unchanged). The formulas above already reflect the corrected source.
- **Sign bug (fixed).** An earlier version of `FT-simplified.tex` computed $\dfrac{d}{d\tau}\!\left[\dfrac{e^{ik(r_0^0+|{\bf R}_0|)}}{|{\bf R}_0|}\right] = ik(u\cdot n_{R_0})\dfrac{e^{ik(\cdots)}}{|{\bf R}_0|} - \dfrac{e^{ik(\cdots)}\,{\bf n}_{R_0}\cdot{\bf u}}{|{\bf R}_0|^2}$ — i.e. with a **minus** sign on the $d/d\tau(1/|{\bf R}_0|)$ term. Direct differentiation gives the opposite: from $|{\bf R}_0|^2={\bf R}_0\cdot{\bf R}_0$ and $d{\bf R}_0/d\tau=-{\bf u}$, $\;d|{\bf R}_0|/d\tau = -{\bf n}_{R_0}\cdot{\bf u}$, so $d(1/|{\bf R}_0|)/d\tau = +({\bf n}_{R_0}\cdot{\bf u})/|{\bf R}_0|^2$ — **plus**, not minus. That sign error propagated through the subsequent integration by parts and flipped the sign of the entire $F_s$ term. The source has been corrected (the $F_s$ formula above now carries the extra minus sign this fix introduces); $F_l$ is unaffected, since it comes from the *other* term in that same derivative (the one built from $\phi'=n_{R_0}\cdot u$), which had the correct sign already. This was found by comparing this "simplified" total field against the independently-derived "direct" form (Form 1): the two agreed to good precision on 5 of the 6 independent Faraday-tensor components, but disagreed by an order of magnitude (and had unrelated phase/sign) specifically on $F^{03}$ — the one component small enough, by transversality, for a near-field ($F_s$) sign error to dominate its total instead of being swamped by the (correct, ~7-orders-larger) $F_l$ contribution that masks this bug everywhere else.

### Form 2's boundary term $F_b$

Writing $\Phi(\tau) = r_0^0(\tau) + |{\bf R}_0(\tau)|$ (the same phase argument shared by $F_l$/$F_s$) and
$T^{\alpha\beta}(\tau) = \dfrac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}{n_{R_0}\cdot u}$ (the tensor factor
shared by $F_l$/$F_s$), the integration by parts that produces $F_l+F_s$ from Jackson's compact form also produces,
unavoidably, a boundary term evaluated at the two ends of the $\tau$-integration range:

$$
F_b^{\alpha\beta}(ck,{\bf x};\tau_{\min},\tau_{\max}) = \frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c^2}
\left[\frac{e^{ik\Phi(\tau)}}{|{\bf R}_0(\tau)|}\,T^{\alpha\beta}(\tau)\right]_{\tau=\tau_{\min}}^{\tau=\tau_{\max}}
$$

so that $F^{\alpha\beta} = F_l^{\alpha\beta} + F_s^{\alpha\beta} + F_b^{\alpha\beta}$ *exactly*, for any $\tau_{\min},
\tau_{\max}$ (finite or infinite) — this piece falls out directly from evaluating $[fg]$ at the integration limits in
$\int f\,g' = [fg] - \int f'g$ and carries no risk of a sign ambiguity from the differentiation step itself. Only in
the $\tau_{\min}\to-\infty,\ \tau_{\max}\to\infty$ limit (with the bracketed quantity decaying at both ends) does
$F_b\to0$, recovering the textbook two-term result. For `compute_radiation`'s finite trajectory window this term is
not automatically small — see the on-axis case below, where it is proven to be individually as large as $F_l$
itself.

### On-axis $F^{03}$ cancellation (worked special case)

The following identity was proven independently in a different project's analysis of this same simplified-form
derivation, and reproduced here since it is directly relevant to `compute_radiation`'s open validation gap — the
implementation it originally referenced (Python module paths, plotting code) is specific to that other project and
is **not** part of this repository; only the mathematical result is transcribed.

For a screen point lying on the particle's own axis of motion (canonical $z$-axis, ${\bf x}=(0,0,\pm|{\bf x}|)$, so
${\bf n}_{R_0}=(0,0,\pm1)$ is constant along the whole trajectory) and the longitudinal component $\alpha=0,\beta=3$:

$$
T^{03}(\tau) = \frac{u^z(\tau) - n_{R_0}^z\,u^0(\tau)}{u^0(\tau) - n_{R_0}^z\,u^z(\tau)}
$$

- **Forward screen** ($n_{R_0}^z=+1$): $T^{03}(\tau) = \dfrac{u^z-u^0}{u^0-u^z} \equiv -1$ — exactly, for every
  $\tau$ and every velocity.
- **Backward screen** ($n_{R_0}^z=-1$): $T^{03}(\tau) = \dfrac{u^z+u^0}{u^0+u^z} \equiv +1$ — exactly.

Because $T^{03}$ is a $\tau$-independent constant $c_\pm=\mp1$ on axis, it factors out of both the $F_l^{03}$
integral and the $F_b^{03}$ boundary evaluation above:

$$
F_l^{03} = \frac{e}{4\pi\epsilon_0c^2}\,c_\pm\!\int_{\tau_{\min}}^{\tau_{\max}}\!\left(\frac{-ik}{|{\bf R}_0(\tau)|}\right)e^{ik\Phi(\tau)}\,d\tau,
\qquad
F_b^{03} = \frac{e}{4\pi\epsilon_0c^2}\,c_\pm\left[\frac{e^{ik\Phi(\tau)}}{|{\bf R}_0(\tau)|}\right]_{\tau_{\min}}^{\tau_{\max}}
$$

and, since $\dfrac{d}{d\tau}e^{ik\Phi(\tau)} = ik\,(u\cdot n_{R_0})\,e^{ik\Phi(\tau)}$, a further integration by parts
on $F_l^{03}$ alone shows $F_l^{03}+F_b^{03}\approx0$ up to the (generally nonzero, $1/|{\bf R}_0|^2$-suppressed)
piece from $d(1/|{\bf R}_0|)/d\tau$ — i.e. **$F_l^{03}$ and $F_b^{03}$ are each individually large on-axis, as an
artifact of the integration by parts, and they cancel each other almost completely.** Form 1 ("direct",
no integration by parts) never produces this artifact: its long-range term is driven by the four-acceleration $w$,
and the corresponding on-axis combination $w^z\mp w^0$ vanishes identically. This is consistent with the two forms'
totals agreeing (per the "Cross-check target" below) while their individual $F_l$ pieces do not, and it means any
implementation of Form 2 that drops $F_b$ — as `compute_radiation` currently does for `radiation_formula=simplified`
(`radiation.cpp`, `long_range_prefactor`/`short_range_prefactor`) — will see a large, spurious $F^{03}$ contribution
near the beam axis that a correct (direct, or boundary-term-complete simplified) implementation would not.

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
compare the **sum** $F_l+F_s$ (Form 1) against $F_l+F_s+F_b$ (Form 2 — see
"Form 2's boundary term $F_b$" above; $F_b$ is *not* optional for a finite
$\tau$-integration range), not the individual $F_l$, $F_s$ pieces, which are
not required to agree between the two derivations — see the note on differing
decompositions above.
