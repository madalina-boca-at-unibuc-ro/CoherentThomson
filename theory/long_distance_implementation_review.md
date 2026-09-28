# Long-distance approximation: reading of the coding guide

Review of [long_distance_direct_simplified_coding_guide.md](long_distance_direct_simplified_coding_guide.md),
checked against the current exact implementation in
[radiation.cpp](../src/core/radiation/radiation.cpp). No inconsistencies were found.

## 1. What the approximation does

For each **(electron, screen point)** pair, fix the reference point at that electron's own first
trajectory sample,

$$
\boldsymbol{\mathfrak R}_0=\mathbf r(\tau_m),\qquad
\mathbf x_0=\mathbf x-\boldsymbol{\mathfrak R}_0,\qquad
\mathbf r_0(\tau)=\mathbf r(\tau)-\boldsymbol{\mathfrak R}_0,\qquad
\mathbf n_0=\frac{\mathbf x_0}{|\mathbf x_0|},\qquad n_0=(1,\mathbf n_0),
$$

and replace the exact per-$\tau$ geometry by a fixed one:

| | exact (current code) | long-distance |
|---|---|---|
| amplitude | $1/\lvert\mathbf x-\mathbf r(\tau)\rvert$, recomputed every $\tau$ | $1/\lvert\mathbf x_0\rvert$, constant |
| direction | $n_{R_0}(\tau)$, recomputed every $\tau$ | $n_0$, constant |
| phase | $k\,\big(x^0+\lvert\mathbf x-\mathbf r(\tau)\rvert\big)$ | $k\lvert\mathbf x_0\rvert+k\,\big(r_0^0-\mathbf n_0\cdot\mathbf r_0\big)$ |
| short-range $F_s\propto 1/R^2$ | kept | **dropped** |

The two phases agree to first order, since
$|\mathbf x_0-\mathbf r_0|\simeq|\mathbf x_0|-\mathbf n_0\cdot\mathbf r_0$.

Because $\boldsymbol{\mathfrak R}_0$ is chosen per electron, the validity condition

$$
\frac{k\,\big|\mathbf r_0-\mathbf n_0(\mathbf n_0\cdot\mathbf r_0)\big|^2}{2|\mathbf x_0|}\ll 1
$$

involves only each electron's own excursion, not the beam size. For `config/config.cfg` (screen at
144000 λ, angles of about 3 mrad) this is far below 1.

## 2. The two formulas

Common prefactor, identical to the current `general_factor` ($e$ = signed charge `q_0`):

$$
C=\frac1{2\pi}\frac{e}{4\pi\epsilon_0c^2},\qquad
\Phi(\tau)=r_0^0-\mathbf n_0\cdot\mathbf r_0 .
$$

**Direct:** the current direct long-range integrand, with $n_{R_0}\to n_0$ and $R\to|\mathbf x_0|$.
There is no $F_s$ and no boundary term:

$$
F_a^{\alpha\beta}=C\,\frac{e^{ik|\mathbf x_0|}}{|\mathbf x_0|}\int_{\tau_m}^{\tau_M}d\tau\,e^{ik\Phi}\,
\frac{(n_0\cdot u)(n_0^\alpha w^\beta-n_0^\beta w^\alpha)-(n_0\cdot w)(n_0^\alpha u^\beta-n_0^\beta u^\alpha)}{(n_0\cdot u)^2}.
$$

**Simplified:** bulk plus boundary, no $F_s$:

$$
F_l^{\alpha\beta}=-ik\,C\,\frac{e^{ik|\mathbf x_0|}}{|\mathbf x_0|}\int_{\tau_m}^{\tau_M}d\tau\,e^{ik\Phi}\,
(n_0^\alpha u_\perp^\beta-n_0^\beta u_\perp^\alpha),
$$

$$
F_b^{\alpha\beta}=C\,\frac{e^{ik|\mathbf x_0|}}{|\mathbf x_0|}
\left[e^{ik\Phi}\,\frac{n_0^\alpha u_\perp^\beta-n_0^\beta u_\perp^\alpha}{n_0\cdot u}\right]_{\tau_m}^{\tau_M},
\qquad u_\perp=\big(0,\ \mathbf u-\mathbf n_0(\mathbf n_0\cdot\mathbf u)\big).
$$

### Why $u_\perp$ instead of $u$, and what it changes

Write $u=u_\perp+\ell$, with $\ell=\big(u^0,\ \mathbf n_0(\mathbf n_0\cdot\mathbf u)\big)$. Then

$$
\ell = (\mathbf n_0\cdot\mathbf u)\,n_0 + (n_0\cdot u)\,t,\qquad t=(1,\mathbf 0),
$$

so $n_0\wedge\ell=(n_0\cdot u)\,(n_0\wedge t)$: using $u$ instead of $u_\perp$ adds a term
proportional to $(n_0\cdot u)$ to both $F_l$ and $F_b$. Since

$$
\frac{d}{d\tau}e^{ik\Phi}=ik\,(n_0\cdot u)\,e^{ik\Phi},
$$

the extra piece in the bulk term is

$$
-ik\int d\tau\,e^{ik\Phi}(n_0\cdot u)\,(n_0\wedge t)=-\big[e^{ik\Phi}\big]_{\tau_m}^{\tau_M}(n_0\wedge t),
$$

and the extra piece in the boundary term is
$+\big[e^{ik\Phi}(n_0\cdot u)/(n_0\cdot u)\big](n_0\wedge t)=+\big[e^{ik\Phi}\big](n_0\wedge t)$.
**They cancel exactly**, so the physical sum $F_l+F_b$ is the same either way.

Consequences:

- The individual `LR_`/`BR_` columns will **not** match the exact code's columns (the guide warns
  about this). Only the sum is comparable.
- With $u_\perp$, each term is purely transverse to $\mathbf n_0$ on its own. The ~1000x on-axis
  $F^{03}$ cancellation between $F_l$ and $F_b$ seen in the exact simplified formula should not occur
  here.

## 3. Two-projection method (guide sections 5–6)

Build a fixed right-handed transverse basis once per (electron, screen point):
$\mathbf e_a\cdot\mathbf n_0=0$, $\mathbf e_1\times\mathbf e_2=\mathbf n_0$, $e_a^\mu=(0,\mathbf e_a)$.
Then $u_\perp=\sum_a(\mathbf e_a\cdot\mathbf u)\,e_a$, and instead of 6 bivector components per
$\tau$ only **two complex scalars** $p_1,p_2$ are accumulated.

Direct:

$$
p_a=-\int_{\tau_m}^{\tau_M}d\tau\,e^{ik\Phi}\,
\frac{(n_0\cdot u)(\mathbf e_a\cdot\mathbf w)-(n_0\cdot w)(\mathbf e_a\cdot\mathbf u)}{(n_0\cdot u)^2}.
$$

Simplified:

$$
p_a=\underbrace{ik\int_{\tau_m}^{\tau_M}d\tau\,(\mathbf e_a\cdot\mathbf u)\,e^{ik\Phi}}_{\to F_l}
\;\underbrace{-\left[\frac{\mathbf e_a\cdot\mathbf u}{n_0\cdot u}\,e^{ik\Phi}\right]_{\tau_m}^{\tau_M}}_{\to F_b}.
$$

Reconstruction (both methods):

$$
F_a^{\alpha\beta}=C\,\frac{e^{ik|\mathbf x_0|}}{|\mathbf x_0|}\sum_{a=1}^2p_a\,
\big(e_a^\alpha n_0^\beta-n_0^\alpha e_a^\beta\big).
$$

**Sign check.** Both leading minus signs come from $n_0\wedge e_a=-\,e_a\wedge n_0$:

- Simplified bulk: $n_0\wedge u_\perp=-\sum_a(\mathbf e_a\cdot\mathbf u)\,(e_a\wedge n_0)$, so
  $-ik\,(n_0\wedge u_\perp)=+ik\sum_a(\mathbf e_a\cdot\mathbf u)(e_a\wedge n_0)$, which matches the $+ik$ in $p_a$.
- Simplified boundary: $+[\dots n_0\wedge u_\perp]$ becomes $-[\dots](e_a\wedge n_0)$, which matches the $-[\,\cdot\,]$ in $p_a$.
- Direct: the numerator $(n_0\cdot u)\,n_0\wedge w_\perp-(n_0\cdot w)\,n_0\wedge u_\perp$ becomes
  $-\sum_a[\dots](e_a\wedge n_0)$, which matches the leading $-$ in $p_a$.

The $(0,j)$ components of the reconstruction basis are
$(e_a\wedge n_0)^{0j}=e_a^0n_0^j-n_0^0e_a^j=-e_a^j$, so $E_j=cF^{j0}=c\,C\,\frac{e^{ik|\mathbf x_0|}}{|\mathbf x_0|}\sum_a p_a e_a^j$
is manifestly transverse.

## 4. How it maps onto the existing code

- **Phase:** $|\mathbf x_0|+r_0^0-\mathbf n_0\cdot\mathbf r_0$ is independent of frequency, so
  `radiation_phase_argument`'s contract and the `cexp *= cexp_step` recurrence still hold. It needs
  no square root per $\tau$ (only a dot product), so this should be faster than the exact calculation.
- **Units:** $r_0^0=ct'$ is `position[0]`, $u$ is `momentum` ($m=1$), $w$ is `acceleration`, and
  $k$ is the `frequencies_list` entry (already $\omega/c$). This is the same as the exact code.
- **Integration variable:** the integral is over proper time $\tau$ with the uniform step `d_tau`,
  so the existing trapezoid weights apply unchanged. $dt/d\tau$ is already inside the formulas, so
  no extra Jacobian is needed.
- **Endpoints:** $\tau_m,\tau_M$ are the first and last stored samples, so the "evaluate at the
  actual limits" rule holds exactly, with no interpolation.
- **Missing references:** the `.tex` and `many_electron_far_field_two_projections.md` files the guide
  links to aren't in this repo. The guide seems self-contained without them.

## 5. Open decisions

1. **How to select it.** Proposed: two new `radiation_formula` values, `long_distance_direct` and
   `long_distance_simplified`. The alternative is a separate on/off key.
2. **Output columns.** Proposed: for simplified, the bulk goes into `LR_`, the endpoint into `BR_`,
   and `SR_` = 0; for direct, everything goes into `LR_`. Every plotting script keeps working
   unchanged.
3. **Implementation route.** Proposed: the two-projection method (2 scalars instead of 6 per $\tau$).
   Validation: the two long-distance formulas should agree with each other, and with the exact code
   wherever the approximation is valid.
4. **Debug mode.** `debug_radiation.cpp` mirrors only the exact simplified formula. Either extend it
   to the long-distance version, or keep debug mode exact-only for now.

## 6. Decisions taken and implementation

- Selected by `radiation_formula = long_distance_simplified | long_distance_direct`
  (`Radiation::RadiationFormula`). The exact `simplified`/`direct` are unchanged.
- Output: `SR_` = 0; simplified puts the bulk in `LR_` and the endpoint in `BR_`; direct puts everything in `LR_`.
- Two-projection method, in `compute_radiation_long_distance` ([radiation.cpp](../src/core/radiation/radiation.cpp)).
- Debug mode was left exact-only at first, and later removed from the code entirely (2026-09-28).

## 7. Validation

Setup: `config/config.cfg` with a rectangular 9x9 screen (±400 λ, backward), 8 electrons, 3 harmonics. All numbers
are max |ΔF| / max |F| on the total field $F_l+F_s+F_b$, per harmonic h1, h2, h3.

**Independent Python evaluation** of $p_a$ from `electron.dat` (1 electron, 3x3 screen) agrees with the C++ at h1
to $3\times10^{-6}$, once a global phase from the 7-digit `fundamental_frequency` in `run_log.txt` is removed.

**Approximation error vs. distance** (`trajectory_NT = 100`):

| $D$ [λ] | LD direct vs exact direct | LD simplified vs LD direct | exact direct vs exact simplified |
|---|---|---|---|
| 144000 | 3.6e-5, 5.9e-5, 6.2e-3 | 1.0e-6, 6.9e-4, 2.6e-1 | 1.7e-4, 1.0e-1, 1.0 |
| 14400 | 3.6e-4, 7.3e-4, 8.2e-3 | 1.0e-5, 1.7e-3, 3.4e-1 | 1.6e-4, 2.5e-2, 1.0 |
| 1440 | 1.0e-2, 1.8e-2, 2.9e-2 | 8.7e-5, 4.3e-3, 4.0e-1 | 1.2e-4, 7.2e-3, 5.4e-1 |
| 144 | 1.2e-1, 2.6e-1, 4.8e-1 | 7.8e-4, 6.0e-2, 1.8 | 8.5e-4, 5.8e-2, 1.0 |

The LD error at h1 scales as $1/D$ (from the amplitude and direction freeze), steepening at small $D$ where the
screen angles grow.

**Convergence in `trajectory_NT`** (at $D = 144000$ λ):

| NT | LD simplified vs LD direct | exact simplified vs exact direct | LD direct vs exact direct |
|---|---|---|---|
| 100 | 1.0e-6, 6.9e-4, 2.6e-1 | 1.7e-4, 1.0e-1, 47 | 3.6e-5, 5.9e-5, 6.2e-3 |
| 200 | 4.9e-7, 1.8e-4, 6.5e-2 | 4.2e-5, 2.9e-2, 12 | 3.6e-5, 5.9e-5, 6.1e-3 |
| 400 | 7.5e-7, 4.7e-5, 1.6e-2 | 1.0e-5, 7.6e-3, 2.9 | 3.6e-5, 5.9e-5, 6.0e-3 |

- Simplified vs direct converges as $O(d\tau^2)$ (4x per doubling) in both geometries, so the weak-harmonic
  discrepancies are quadrature error, not a bug.
- LD simplified converges ~200x better than exact simplified at h3, which is the $u_\perp$ benefit expected in
  section 2.
- LD vs exact does not change with NT: that is the true approximation error.

**Speed:** ~0.017 s/electron (LD simplified) vs ~0.032 (exact simplified) on the 9x9 screen, about 2x faster.

**Compiler issue:** GCC 15.2's SLP vectorizer miscompiled this function under `-O3 -march=native`, giving wrong
harmonics after the first. It is now disabled for this function only; see `theory/dev_notes.md` (2026-09-28).

## 8. Stepped phase: `long_distance_simplified_approx` / `long_distance_direct_approx`

**Why.** In both the exact and the plain long-distance code, the cost per (screen point, $\tau$) is dominated by
the complex exponential $e^{ik\Phi}$. Its argument is $k\Phi\sim10^{7}$ rad, and one `std::polar` costs ~19 ns.
So the long-distance approximation alone gave only ~1.6x.

**Method.** With $\mathbf n_0$ fixed, the phase is linear along the trajectory:

$$
\Phi_{j+1}-\Phi_j=\Delta\Phi_j=\Delta r^0_j-\mathbf n_0\cdot\Delta\mathbf r_j ,\qquad
e^{ik\Phi_{j+1}}=e^{ik\Phi_j}\,e^{ik\Delta\Phi_j}.
$$

$\Delta\Phi_j$ comes from differences of stored samples (small numbers, so accurate). The factor
$e^{ik\Delta\Phi_j}$ is a Taylor polynomial (cos to order 10, sin to order 11), accurate to $10^{-16}$ for
$|k\Delta\Phi|\le 0.25$. Larger angles are halved $m$ times and the result squared $m$ times ($m\le6$).
The phasor is re-anchored from the exact phase every 64 steps. Eight screen points are advanced together so the
CPU overlaps their independent chains; a single chain is latency-bound and barely beats `std::polar`. For several
frequencies, both the base phasor $e^{ik_0\Phi}$ and the step phasor $e^{i\,\Delta k\,\Phi}$ are advanced, and
the usual frequency recurrence is applied at each $\tau$.

Because $k_1\propto 1/(n_0\cdot u)$ (Doppler), the per-step angle is about $2\pi N/\text{trajectory\_NT}$ for
harmonic $N$, independent of geometry.

**Accuracy.** Tested against an extended-precision (80-bit) evaluation on identical inputs:

| | direct phase (LD) | stepped phase (approx) |
|---|---|---|
| h1 | ~5e-10 | ~7e-10 |
| h2, h3 | 4e-8 – 8e-7 | 5e-8 – 8e-7 |
| h5 – h7 | 6e-8 – 1.4e-5 | 2e-8 – 1.7e-6 |

Both are limited by float64 rounding of $k\Phi$, amplified by the cancellation in the $\tau$ sum
($\sum|\text{terms}|/|\text{sum}|$ up to $10^4$ at weak harmonics). So approx and LD differ by 1e-5 – 1e-4 at
weak harmonics in full runs, while agreeing to ~1e-8 – 1e-10 at h1. Release and `-O0` builds agree.

**Speed** (`config.cfg`: circular 64x64, 20 threads; ns per (screen point, $\tau$)):

| | exact | long_distance | long_distance_approx |
|---|---|---|---|
| simplified, h1 | 58 | 36 | **11** |
| direct, h1 | 65 | 41 | **18** |
| simplified, h5–h7 | 112 | 74 | **24** |

## 9. Constant phase moved out of the $\tau$ loop

Both long-distance paths used to put $|\mathbf x_0|$ inside the per-$\tau$ phase argument, i.e. evaluate
$e^{ik(|\mathbf x_0|+r^0-\mathbf n_0\cdot\mathbf r_0)}$ at each step. That is mathematically fine, but
$|\mathbf x_0|\approx2\times10^9$ a.u. dominates the argument. A double then resolves the varying part only to
~$5\times10^{-7}$ a.u., i.e. ~$10^{-9}$ rad of random phase per step, and the up-to-$10^4$ cancellation of the
$\tau$ sum at weak harmonics amplifies that.

Now the phase is split as

$$
e^{ik\Phi}=\underbrace{e^{ik\,(|\mathbf x_0|+r^0(\tau_m))}}_{\text{once per frequency, at reconstruction}}\;
\underbrace{e^{ik\,\left(r^0-r^0(\tau_m)-\mathbf n_0\cdot\mathbf r_0\right)}}_{\text{inside the }\tau\text{ loop}} .
$$

The constant factor still carries ~$10^{-9}$ rad of rounding, but as one overall phase per (electron, screen
point) after the sum, where nothing amplifies it. The time is not reset: $r^0(\tau_m)$ sits in the constant
factor, so relative phases between electrons are unchanged.

float64 error vs. the extended-precision reference (same inputs):

| | constant inside (before) | constant outside (now) |
|---|---|---|
| h1 | 5e-10 – 6e-10 | 6e-10 – 1e-9 |
| h3 | 4e-8 – 8e-7 | ~2e-9 |
| h6 | 3e-6 – 1.4e-5 | 3e-9 – 5e-9 |

Approx vs. unstepped LD in full runs went from ≤3e-4 to ≤1e-7 (Release). Results are otherwise unchanged
(LD vs exact direct: 3.6e-5 at h1), and so is the speed.
