# Implementation details: radiation computation modes and optimizations

This document describes how `Radiation::compute_radiation`
([radiation.cpp](../src/core/radiation/radiation.cpp)) evaluates the Fourier-transformed radiated field, the six
computation modes selected by the config key `radiation_formula`, and the numerical and performance
optimizations behind them. The derivations are in
[FT_Faraday_tensor-direct_and_simplified_forms-v2.md](FT_Faraday_tensor-direct_and_simplified_forms-v2.md)
(exact formulas) and [long_distance_direct_simplified_coding_guide.md](long_distance_direct_simplified_coding_guide.md)
(long-distance formulas). Validation data is in
[long_distance_implementation_review.md](long_distance_implementation_review.md).

## 1. What is computed

For every electron, frequency $k=\omega/c$ and screen point $\mathbf x$, the code evaluates the Faraday tensor
$F^{\alpha\beta}(ck,\mathbf x)$ of the radiated field as an integral over the electron's proper time $\tau$,
and sums it coherently over all electrons:

$$
F^{\alpha\beta}_{\text{total}}(ck,\mathbf x)=\sum_{\text{electrons}}F^{\alpha\beta}(ck,\mathbf x).
$$

Common ingredients, shared by all modes:

| quantity | meaning | where it comes from |
|---|---|---|
| $r^\mu(\tau)$ | position, $r^0=ct'$ | `Electron::State::position` (RK4 trajectory) |
| $u^\mu=dr^\mu/d\tau$ | four-velocity ($m=1$) | `State::momentum` |
| $w^\mu=du^\mu/d\tau$ | four-acceleration | `State::acceleration` (only the direct modes use it) |
| $\tau_j$, $d\tau$ | uniform proper-time grid, $N_\tau$ samples | `Electron::get_d_tau`, `get_N_tau` |
| $k$ | $\omega/c$ | `frequencies_list` (harmonics, or a dense linear scan) |
| $C=\dfrac{q_0}{2\pi\,4\pi\epsilon_0c^2}$ | overall prefactor, $q_0=-1$ signed charge | `general_factor`, applied once in `run_simulation` |

**Quadrature.** Integrals over $\tau$ use the trapezoidal rule: weight $d\tau$, halved at the two endpoints.
Boundary (endpoint) terms are exact evaluations and get no weight.

**Output.** Every mode fills three packed tensors per (frequency, screen point), exported as the
`LR_`/`SR_`/`BR_` column groups of `radiation_field.dat`:

- long range $F_l$ (`LR_`),
- short range $F_s$ (`SR_`),
- boundary $F_b$ (`BR_`).

**Only the sum $F_l+F_s+F_b$ is physical.** How the field is split between the three differs between modes.

## 2. The six modes

| `radiation_formula` | geometry | integrand needs | `LR_` | `SR_` | `BR_` |
|---|---|---|---|---|---|
| `simplified` (default) | exact | $u$ | bulk | $1/R^2$ term | endpoints |
| `direct` | exact | $u$, $w$ | $1/R$ term | $1/R^2$ term | 0 |
| `long_distance_simplified` | long distance | $u$ | bulk | 0 | endpoints |
| `long_distance_direct` | long distance | $u$, $w$ | everything | 0 | 0 |
| `long_distance_simplified_approx` | long distance, stepped phase | $u$ | bulk | 0 | endpoints |
| `long_distance_direct_approx` | long distance, stepped phase | $u$, $w$ | everything | 0 | 0 |

An unknown value stops the run with an error listing the valid names.

### 2.1 Exact modes: `simplified` and `direct`

For each $\tau$ the exact distance and direction to the screen point are recomputed:

$$
\mathbf R(\tau)=\mathbf x-\mathbf r(\tau),\qquad R=|\mathbf R|,\qquad
n=(1,\ \mathbf R/R),\qquad \text{phase } e^{ik\,(r^0+R)} .
$$

With $n\wedge u$ denoting the bivector $n^\alpha u^\beta-n^\beta u^\alpha$:

**`simplified`** (integration by parts; velocity only):

$$
F_l=-ik\,C\int d\tau\,e^{ik(r^0+R)}\,\frac{n\wedge u}{R},\qquad
F_s=C\int d\tau\,e^{ik(r^0+R)}\,\frac{s\wedge u}{R^2},\ \ s=(0,\mathbf n),\qquad
F_b=C\left[e^{ik(r^0+R)}\,\frac{n\wedge u}{R\,(n\cdot u)}\right]_{\tau_m}^{\tau_M}.
$$

**`direct`** (Fourier transform of the Liénard–Wiechert field; needs the acceleration):

$$
F_l=C\int d\tau\,e^{ik(r^0+R)}\,
\frac{(n\cdot u)\,n\wedge w-(n\cdot w)\,n\wedge u}{R\,(n\cdot u)^2},\qquad
F_s=C\int d\tau\,e^{ik(r^0+R)}\,\frac{c^2\,n\wedge u}{R^2\,(n\cdot u)^2}.
$$

The explicit $c^2$ in the direct $F_s$ compensates for the $1/c^2$ contained in $C$.

The two exact modes agree to the trapezoid error $O(d\tau^2)$. Comparing them is the main correctness check.

### 2.2 Long-distance modes: `long_distance_simplified` and `long_distance_direct`

For each **(electron, screen point)** pair, the geometry is frozen at the electron's own first sample:

$$
\boldsymbol{\mathfrak R}_0=\mathbf r(\tau_m),\quad
\mathbf x_0=\mathbf x-\boldsymbol{\mathfrak R}_0,\quad
\mathbf r_0(\tau)=\mathbf r(\tau)-\boldsymbol{\mathfrak R}_0,\quad
\mathbf n_0=\frac{\mathbf x_0}{|\mathbf x_0|},\quad n_0=(1,\mathbf n_0),
$$

with the amplitude $1/|\mathbf x_0|$, the linearized phase
$\Phi=|\mathbf x_0|+r^0-\mathbf n_0\cdot\mathbf r_0$, and the $1/R^2$ term dropped.
Only space is shifted; time is not reset, so relative phases between electrons are kept.

**Validity.** $k\,|\mathbf r_{0\perp}|^2/(2|\mathbf x_0|)\ll1$ over each electron's own excursion. The beam size
does not enter, because the reference point is per electron. Measured on `config.cfg`, the error against exact
`direct` scales as $1/D$: $3.6\times10^{-5}$ at the fundamental for a screen at 144000 λ, and $10^{-2}$ at 1440 λ.
This error does not depend on `trajectory_NT`.

**Two-projection method.** A fixed right-handed transverse basis $\mathbf e_1,\mathbf e_2\perp\mathbf n_0$ is built
once per (electron, screen point): the Cartesian axis least aligned with $\mathbf n_0$ is projected out and
normalized, then $\mathbf e_2=\mathbf n_0\times\mathbf e_1$. The tensor is carried by two complex scalars:

$$
F=C\,\frac{e^{ik(|\mathbf x_0|+r^0(\tau_m))}}{|\mathbf x_0|}\sum_{a=1,2}p_a\,(e_a\wedge n_0),\qquad e_a=(0,\mathbf e_a),
$$

$$
\text{simplified: }\ p_a=ik\int d\tau\,(\mathbf e_a\cdot\mathbf u)\,e^{ik\phi}
-\left[\frac{\mathbf e_a\cdot\mathbf u}{n_0\cdot u}\,e^{ik\phi}\right]_{\tau_m}^{\tau_M},
\qquad
\text{direct: }\ p_a=-\int d\tau\,e^{ik\phi}\,
\frac{(n_0\cdot u)(\mathbf e_a\cdot\mathbf w)-(n_0\cdot w)(\mathbf e_a\cdot\mathbf u)}{(n_0\cdot u)^2},
$$

where $\phi=r^0-r^0(\tau_m)-\mathbf n_0\cdot\mathbf r_0$ is the in-loop phase (see §4.3). For simplified, the
integral part goes to `LR_` and the endpoint part to `BR_`.

Because only the transverse part of $u$ appears ($u_\perp$), each term is transverse to $\mathbf n_0$ on its own.
The large on-axis $F^{03}$ cancellation between `LR_` and `BR_` of the exact simplified mode does not occur. As a
consequence, the long-distance simplified mode converges much faster in `trajectory_NT` at weak harmonics.

### 2.3 Stepped-phase modes: `*_approx`

These are the same long-distance formulas with the same outputs; only the evaluation of $e^{ik\phi}$ changes
(§4.4). They agree with the unstepped long-distance modes to $\lesssim10^{-7}$ and are 3–4× faster.

### 2.4 Choosing a mode

- **Production runs at a distant screen:** `long_distance_simplified_approx`. It is the fastest mode, needs no
  acceleration, and is the most accurate at weak harmonics.
- **Checking a long-distance result:** rerun with `long_distance_direct_approx`. The two should agree to
  $O(d\tau^2)$.
- **Screens close to the beam, or checking the approximation itself:** `direct` and `simplified`. For
  harmonics beyond the fundamental, prefer `direct` or raise `trajectory_NT`, because the exact simplified mode
  converges slowly there (§5).

## 3. Parallel structure and loop order

- **Threads.** `run_simulation` splits the beam into contiguous electron chunks, one per thread
  (`num_threads=0` uses every hardware thread). Each thread accumulates into its own private field, and the
  fields are summed at the end, so there is no locking. The whole beam is generated upfront.
- **Loop order: screen point outer, $\tau$ inner.** Each screen point's $\tau$ sum accumulates in small local
  buffers that stay in cache, and the large `[N_freq][N_screen]` output array is written once per
  (screen point, frequency), not once per $\tau$.
- **Packed tensors.** Only the 6 independent components ($01,02,03,12,13,23$) of each antisymmetric tensor are
  accumulated (`PackedFaraday`). The full $4\times4$ tensor is rebuilt once, after the thread reduction.
- **Two scalars instead of six.** In the long-distance modes the $\tau$ loop accumulates only $p_1,p_2$ per
  frequency. The six tensor components are formed once per (screen point, frequency) at reconstruction.

## 4. Optimizations of the phase factor

For each (screen point, $\tau$) every mode needs the complex exponential $e^{ik\phi}$, and it dominates the run
time. At a distant screen its argument is large ($k\Phi\sim10^{7}$ rad), and one `std::polar` (a sine and a
cosine) costs about 19 ns, more than all the other per-step arithmetic of the long-distance modes together.

### 4.1 Frequency recurrence (all modes)

`frequencies_list` is always an arithmetic progression $k_i=k_0+i\,\Delta k$ (harmonics, or a linear scan). So

$$
e^{ik_i\phi}=e^{ik_0\phi}\,\big(e^{i\Delta k\,\phi}\big)^i ,
$$

and each $\tau$ needs two exponentials whatever the number of frequencies; every further frequency costs one
complex multiplication. This requires the phase argument to be independent of frequency, which is why the
prefactors are kept separate from the phase (`radiation_phase_argument`). A per-electron check falls back to
one exponential per frequency if the list is not evenly spaced.

### 4.2 Why the long-distance approximation alone gives only ~1.6×

The long-distance modes remove the square root, the division and the six-component bivector algebra of the exact
modes, but they still evaluate $e^{ik\phi}$ at every (screen point, $\tau$). Since that exponential is the
dominant cost, the gain is limited to about 1.6× (36 ns vs 58 ns per step on `config.cfg`).

### 4.3 Constant phase outside the $\tau$ loop (long-distance modes)

The full phase is split into a constant and a varying part:

$$
e^{ik\Phi}=\underbrace{e^{ik\,(|\mathbf x_0|+r^0(\tau_m))}}_{\text{once per frequency and screen point}}\ \cdot\
\underbrace{e^{ik\,(r^0-r^0(\tau_m)-\mathbf n_0\cdot\mathbf r_0)}}_{\text{inside the }\tau\text{ loop}} .
$$

This changes no physics and costs no speed. It is purely a precision fix:

- $|\mathbf x_0|\approx2\times10^9$ a.u. would otherwise dominate the in-loop argument, and a double resolves the
  varying part only to about $5\times10^{-7}$ a.u., i.e. about $10^{-9}$ rad of random phase error per step.
- At weak harmonics the $\tau$ sum cancels by up to $10^4$ (the sum of the magnitudes of the terms is up to $10^4$
  times the result), and that cancellation amplified the per-step noise to $\sim10^{-5}$ relative.
- With the constant outside, its rounding is a single overall phase per (electron, screen point) that nothing
  amplifies. The float64 error at h6 fell from $1.4\times10^{-5}$ to $5\times10^{-9}$ against an 80-bit reference.

The exact modes still carry $R\approx|\mathbf x_0|$ inside the loop, so their weak-harmonic results are limited
to about $10^{-5}$ by this rounding.

### 4.4 Stepped phase (the `*_approx` modes)

With $\mathbf n_0$ fixed, the long-distance phase is linear along the trajectory, so between consecutive samples

$$
\phi_{j+1}-\phi_j=\Delta\phi_j=\Delta r^0_j-\mathbf n_0\cdot\Delta\mathbf r_j,\qquad
e^{ik\phi_{j+1}}=e^{ik\phi_j}\,e^{ik\Delta\phi_j}.
$$

$\Delta r^0_j,\Delta\mathbf r_j$ are differences of neighbouring stored samples (small, accurate numbers),
precomputed once per electron. The pieces:

1. **Small-angle phasor.** $e^{ia}$ for $|a|\le0.25$ from a Taylor polynomial: cos to order 10 and sin to
   order 11, accurate to $10^{-16}$ (`small_angle_phasor`).
2. **Angle halving.** Because the Doppler-shifted fundamental satisfies $k_1\propto1/(n_0\cdot u)$, the step angle
   of harmonic $N$ is about $2\pi N/\text{trajectory\_NT}$ (0.063 for $N=1$, NT = 100), whatever the geometry.
   From $N\approx4$ the angle exceeds 0.25. The polynomial is then evaluated at $a/2^m$ and squared $m$ times
   ($m\le6$, i.e. up to 16 rad per step); beyond that the step falls back to `std::polar`.
3. **Re-anchoring.** Every 64 steps the phasor is reset from the exact in-loop phase with `std::polar`, which
   bounds the accumulated rounding error.
4. **Blocks of 8 screen points.** A single phasor chain is limited by the latency of its dependent
   multiplications: about 12 ns per step, barely better than `std::polar`. Eight screen points are advanced
   together as independent chains, which the CPU overlaps and the compiler turns into SIMD (AVX2) loops, bringing
   the cost to about 3.5 ns per step. The loops are written branch-free for this; the halving count is common to
   the block.
5. **Several frequencies.** Both the base phasor $e^{ik_0\phi}$ and the step phasor $e^{i\Delta k\,\phi}$ are
   advanced, and the frequency recurrence of §4.1 is applied at each $\tau$. With a single frequency the step
   phasor is skipped.
6. **Work skipped per step.** Endpoint terms (with their divisions) are computed only at the two endpoints.

**Accuracy.** The stepped phase is at least as accurate as the directly evaluated one; at high harmonics it is
more accurate, because its increments are computed from small numbers. The approx and unstepped long-distance
modes agree to $\lesssim10^{-7}$, mostly to all 7 digits printed in `radiation_field.dat`.

## 5. Precision limits worth knowing

- **Weak harmonics cancel strongly.** In the $\tau$ sum at h2 and above, the result is up to $10^4$ times smaller
  than the sum of the magnitudes of its terms, so any per-step error is amplified by that factor. This affects
  both the rounding of the phase (§4.3) and the quadrature error.
- **Exact `simplified` at weak harmonics.** Its $F_l$ and $F_b$ are individually large and cancel, so its
  trapezoid error is amplified. On `config.cfg` with `trajectory_NT=100`, it differs from `direct` by 10% at h2
  and about 100% at h3. The difference shrinks 4× per doubling of `trajectory_NT`, so it is quadrature error,
  not a bug.
- **Output precision.** `radiation_field.dat` prints 7 significant digits, and `run_log.txt` prints
  `fundamental_frequency` to 7 digits. Reconstructing phases from these files, e.g. in an independent Python
  check, carries a global phase error of order $k\Phi\times10^{-7}$.

## 6. Compiler caveat

GCC 15.2's SLP vectorizer, combined with FMA contraction under `-O3 -march=native`, miscompiled the complex
accumulation of `compute_radiation_long_distance`: every frequency after the first came out wrong (h2 about 100×,
h3 about $10^4$×). Sanitizers found no undefined behaviour, and the logged per-step inputs were correct. SLP is
disabled for the two long-distance functions only, with `__attribute__((optimize("no-tree-slp-vectorize")))`;
loop vectorization, which the stepped mode relies on, stays enabled. The exact modes were not affected.

When changing these functions, or moving to another compiler or machine, rerun a long-distance mode with
`N_harmonics >= 2` and compare against a `-O0` build. All build types write the same
`bin/coherent_thomson_solver`, and a build with nothing to recompile does not relink, so copy each binary out
after a `--clean-first` build before comparing.

GCC 13.3 (on a server) was checked and does not show the miscompile. The workaround is kept on every machine
anyway, deliberately, rather than made compiler-conditional.

The stepped mode must stay vectorizable: after editing it, check with `-fopt-info-vec` that the stepping loops
(`advance_phasors`) still vectorize.

## 7. Performance summary

`config.cfg` with a circular 64×64 screen, 20 threads (i7-12700H), Release build. The table gives CPU time per
(electron, screen point, $\tau$ step), summed over threads:

| mode | h1 | h5–h7 (3 harmonics) |
|---|---|---|
| `simplified` | 58 ns | 112 ns |
| `direct` | 65 ns | 119 ns |
| `long_distance_simplified` | 36 ns | 74 ns |
| `long_distance_direct` | 41 ns | 79 ns |
| `long_distance_simplified_approx` | **11 ns** | **24 ns** |
| `long_distance_direct_approx` | **18 ns** | **31 ns** |

The direct approx mode is slower than the simplified one because its integrand needs a division by
$(n_0\cdot u)^2$ at every step. The trajectory integration (RK4) is negligible next to the radiation sum: a few ms
per electron, against about 0.4 s for a 4096-point screen.

`run_simulation` prints per-thread timings and CPU-seconds per electron (per screen point), computed from the
**sum of per-thread times**, not wall-clock time. The hot loops are deliberately left uninstrumented. On a matched
config the Release build is about 4× faster than the Python reference; an older "C++ is 3× slower" result came
from a Debug build, which is 12–13× slower.
