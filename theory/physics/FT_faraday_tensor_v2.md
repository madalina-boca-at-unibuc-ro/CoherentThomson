# Fourier transform of the Faraday tensor — two equivalent forms

Sources: `shared_src/elm_el_sc/frequency-domain/FT-direct.tex` (direct
formula), `shared_src/elm_el_sc/frequency-domain/FT-simplified-v2.tex`
(corrected potential-based derivation), and
`md_helpers/fourier_faraday_clean_derivation.md`, section 7.
The previous simplified short-range formula was incorrect, independently
of the earlier Jacobian correction. The corrected formulas below supersede it.
Existing notation and the form labels are retained for the coding agent;
the numbered tags are legacy identifiers, not current TeX equation numbers.

Related implementation notes: units.md defines the atomic-unit convention,
equations_of_motion.md defines the trajectory solver, and
emitted_field_time.md defines the corresponding time-domain field.

**Task for the coding agent:** correct the formula for $F^{\alpha\beta}$
in the existing simplified implementation using Form 2 below. Keep the
existing notation, interfaces, trajectory data, integration bounds, phase
calculation, and numerical integration machinery unchanged. The direct
implementation is the reference for comparison and requires no changes.

The required changes are the simplified short-range kernel and the endpoint
contribution to its total. The simplified long-range kernel is algebraically
unchanged.

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
- **Jacobian of the change of variable $t\to\tau$.** The required Jacobian is $dt/d\tau=(u\cdot n_{R_0})/c$. The displayed formulas already include this Jacobian; do not multiply them by it again.
- The direct form has $F=F_l+F_s$. On a finite interval, the simplified
  form has $F=F_l+F_s+F_b$, including the endpoint tensor below.
  **The long/short-range split differs between the derivations.** Only
  the total tensors must agree; do not compare the separate pieces as
  though they were identical.

All quantities follow units.md: $q$ is the signed charge and $c=1/\alpha$.
Do not reuse $q$ for $n_{R_0}\cdot u$; do not set $c=1$.
The Fourier and tensor conventions are

$$
F^{\alpha\beta}(\omega,{\bf x})=
\frac1{2\pi}\int_{t_m}^{t_M}dt\,e^{i\omega t}F^{\alpha\beta}(t,{\bf x}),
\qquad F^{\alpha\beta}=\partial^\alpha A^\beta-\partial^\beta A^\alpha,
$$

with $E_i=cF^{i0}$ and $B_i=-\epsilon_{ijk}F^{jk}/2$. In these
formulas ${\bf x}$ denotes the same observation point as ${\bf x}_0$.
The bounds $\tau_m,\tau_M$ in the formulas below are the bounds already
implemented in the code. No new endpoint calculation is required. The
displayed kernels already include the change-of-variable Jacobian.

---

## Form 1 — "direct" (from `FT-direct.tex`)

Starting point: direct Fourier transform of the closed-form primary expression
for $F^{\alpha\beta}(x)$ (Liénard–Wiechert field written without an explicit
$\tau$-derivative).

$$
F^{\alpha\beta}(ck,{\bf x}) = F^{\alpha\beta}_{l}(ck,{\bf x}) + F^{\alpha\beta}_{s}(ck,{\bf x})
\tag{IV.1.1.14}
$$

$$
F_l^{\alpha\beta}(ck,{\bf x}) =\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c^2}
\int_{\tau_m}^{\tau_M} d\tau\; e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}\;
\frac{1}{|{\bf R}_0(\tau)|}
\frac{(u\cdot n_{R_0})(n_{R_0}^{\alpha}w^{\beta}-n_{R_0}^{\beta}w^{\alpha})-(w\cdot n_{R_0})(n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha})}
{(u\cdot n_{R_0})^2}
\tag{IV.1.1.15}
$$

$$
F_s^{\alpha\beta}(ck,{\bf x}) =\frac{e}{8\pi^2\epsilon_0}
\int_{\tau_m}^{\tau_M} d\tau\; e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}\;
\frac{1}{|{\bf R}_0(\tau)|^2}
\frac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}
{(u\cdot n_{R_0})^2}
\tag{IV.1.1.16}
$$

Notes for implementation:

- $F_l$ scales as $1/|{\bf R}_0|$ and carries the atomic-unit prefactor
  $q/c^2$; $F_s$ scales as $1/|{\bf R}_0|^2$ and carries $q$.
- Both integrands are antisymmetric in $\alpha\leftrightarrow\beta$ by
  construction (each is built from $n_{R_0}^{\alpha}(\cdot)^\beta - n_{R_0}^\beta(\cdot)^\alpha$-type combinations), so only the 6 independent components of the antisymmetric tensor need to be evaluated.
- Denominator power: $(u\cdot n_{R_0})^2$ in both terms (**not** $^3$ — see below).
- **Do not confuse $n_{R_0}(\tau)$ with $n_0$.** $n_{R_0}(\tau)=(1,{\bf n}_{R_0}(\tau))$ is the *exact*, $\tau$-dependent unit direction from the emitting charge to the observation point, defined above and used throughout this note. $n_0=(1,{\bf n}_0)$, ${\bf n}_0={\bf x}_0/|{\bf x}_0|$, is a *different*, constant vector used only in separate long-distance/far-field expansions; it never appears in the exact Fourier-transform formulas here. An earlier version of the `FT-direct.tex` source used `n_0` as shorthand for `n_{R_0}` in this section — that was a notation bug in the LaTeX source and has since been fixed to write `n_{R_0}` explicitly; this note now matches the corrected source.
- **Jacobian bug (fixed).** Restoring the previously dropped $dt/d\tau$ changed the atomic-unit prefactor from $q/c$ to $q/c^2$ and reduced the denominator power from $(u\cdot n_{R_0})^3$ to $(u\cdot n_{R_0})^2$. The displayed formulas contain the correction.

---

## Form 2 — "simplified" (corrected from `FT-simplified-v2.tex`)

Starting point: Fourier-transform the Liénard–Wiechert potentials, keeping
both the temporal integration-by-parts endpoint and the spatial dependence
of the retarded proper-time limits. Alternatively, Jackson's fixed-event
proper-time derivative can be used only with the explicit chain-rule
correction. Replacing it directly by a derivative along the retarded
observation path loses finite-distance terms.

To express the corrected short-range tensor, define one additional vector:

$$
s^\mu(\tau)=(0,{\bf n}_{R_0}(\tau)).
$$

Thus $s^0=0$, $s^i=n_{R_0}^i$ for $i=1,2,3$, and
$s_\mu=(0,-{\bf n}_{R_0})$. These are contravariant components in
the integrands. The distinction from $n_{R_0}^0=1$ is essential:
only spatial differentiation of the potential amplitude produces $1/|{\bf R}_0|^2$.
No existing trajectory or direction variable needs to be renamed.

$$
F^{\alpha\beta}(ck,{\bf x}) = F^{\alpha\beta}_{l}(ck,{\bf x})
+ F^{\alpha\beta}_{s}(ck,{\bf x})+F^{\alpha\beta}_{b}(ck,{\bf x})
\tag{IV.1.2.13}
$$

$$
F^{\alpha\beta}_{l}(ck,{\bf x}) =
\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c^2}\int_{\tau_m}^{\tau_M}d\tau\;
e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}
\left(-\frac{ik}{|{\bf R}_0|}\right)
\left(n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}\right)
\tag{IV.1.2.14}
$$

$$
F^{\alpha\beta}_{s}(ck,{\bf x}) =
\frac{1}{2\pi}\frac{e}{4\pi\epsilon_0c^2}\int_{\tau_m}^{\tau_M}d\tau\;
e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}
\frac{s^{\alpha}u^{\beta}-s^{\beta}u^{\alpha}}{|{\bf R}_0|^2}
\tag{IV.1.2.15}
$$

$$
\boxed{
F_b^{\alpha\beta}(ck,{\bf x})=
\frac{e}{8\pi^2\epsilon_0 c^2}
\left[
e^{ik\left(r_0^0(\tau)+|{\bf R}_0(\tau)|\right)}
\frac{n_{R_0}^{\alpha}u^{\beta}-n_{R_0}^{\beta}u^{\alpha}}
{|{\bf R}_0|(n_{R_0}\cdot u)}
\right]_{\tau_m}^{\tau_M}.
}
$$

Notes for implementation:

- Preserve the existing $F_l$, $F_s$ notation. Include $F_b$ exactly once in the
  simplified total (correct an existing endpoint term if present); do not add it to the direct total, which already
  integrates the actual field without integration by parts.
- The simplified $F_l$ is algebraically unchanged: its previous factors
  $(u\cdot n_{R_0})$ and $1/(n_{R_0}\cdot u)$ cancel exactly.
- **Replace the simplified $F_s$ kernel.** Remove the old scalar factor
  $-({\bf n}_{R_0}\cdot{\bf u})/(n_{R_0}\cdot u)$ and use
  $s^\alpha u^\beta-s^\beta u^\alpha$ instead of
  $n_{R_0}^\alpha u^\beta-n_{R_0}^\beta u^\alpha$.
  There is no $n_{R_0}\cdot u$ denominator in this corrected bulk term.
- The two simplified bulk terms no longer share a single tensor factor.
  For spatial indices $i,j$, their numerators coincide. For an electric
  component with indices $i0$, they are respectively
  $n_{R_0}^i u^0-u^i$ and $n_{R_0}^i u^0$; the $0i$ components have
  the opposite signs. Respect the code's stored index order.
- The endpoint sign is **positive**, with the bracket meaning upper
  endpoint minus lower endpoint. Use the same endpoint trajectory data,
  phase, and observation geometry as in the quadrature. The endpoint
  term requires no acceleration and no extra Jacobian.
- All three simplified contributions carry $q/(2\pi c^2)$ in atomic
  units. The earlier Jacobian fix alone did not repair the short-range
  kernel or eliminate the need for the endpoint term.
- Antisymmetry still allows evaluation of only six independent components.
  The simplified form needs $r_0,u$; the direct form also needs $w$.
- These formulas do not divide by $k$ and are valid at $k=0$ on a finite
  interval. Do not implement a spurious zero-frequency singularity.

---

## Shared notation block

$$
k = \frac{\omega}{c},\qquad
R_0 = (|{\bf R}_0|,{\bf R}_0) = |{\bf R}_0|(1,{\bf n}_{R_0}) = |{\bf R}_0|\,n_{R_0},\qquad
{\bf R}_0 = {\bf x}_0-{\bf r}_0(\tau),\qquad
{\bf n}_{R_0} = \frac{{\bf R}_0}{|{\bf R}_0|}.
\tag{IV.1.1.17; IV.1.2.16}
$$

In all integrals, every quantity built from $u$, $w$, ${\bf r}_0$, ${\bf R}_0$, $n_{R_0}$ is evaluated **as a function of the integration variable $\tau$** (not at a fixed retarded time) — the retarded-time constraint has already been absorbed into the change of variable $t\to\tau$ and determines the phase $e^{ik(r_0^0(\tau)+|{\bf R}_0(\tau)|)}$, the integration limits, and the endpoint evaluations.

## Cross-check target

For the same trajectory ${\bf r}_0(\tau)$, observation point ${\bf x}_0$,
frequency, and endpoints, verify

$$
\boxed{
(F_l+F_s)_{\mathrm{direct}}
=(F_l+F_s+F_b)_{\mathrm{simplified}}.
}
$$

Compare the **total** tensors, not the individual long/short pieces.
The formulas are exact at finite distance away from the source; do not
replace $n_{R_0}(\tau)$ by a constant direction in either implementation
for this cross-check.

Use the existing trajectory, endpoints, and numerical settings for this
comparison, with absolute and relative tolerances. This check concerns the
complete tensors displayed above.
