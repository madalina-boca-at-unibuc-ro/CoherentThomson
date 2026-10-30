# Ponderomotive lens: the dressing as a dented source

Written 2026-10-30. A physical picture of the radial phase that ponderomotive dressing puts on the emission
(`lg_decomposition_implementation.md`, Verifications 4 and 6), a parameter-free model of it, and the tests
that would confirm it. **Test 1 is done and confirms the model** (Results below): the dressing factor is in
`lg_decomposition.py`. Tests 2–4 are still open.

## Picture (electrons at rest)

- **The disk dents.** Inside the pulse, each electron drifts along the laser direction by an amount ∝ its local
  `⟨a²⟩`. With the laser phase `ξ` and `A(r, ξ) = ⟨a²⟩` at the electron, the displacement is
  `s(r, ξ) = (λ_L/2π) ∫ A/2 dξ'`. The flat electron disk becomes a dent shaped like the intensity profile:
  `∝ e^{−2r²/w0²}` for a Gaussian, parabolic near the axis and Gaussian further out. It deepens linearly in time
  while the pulse passes.
- **Each electron's emission is phase-locked to the laser.** It emits N harmonic cycles per laser cycle, so every
  radius emits the same number of cycles.
- **For a screen facing the laser (our backscatter geometry at rest),** observer time runs as
  `dt_obs/dξ = (1 + A)/ω_L`. The electron's wave train is stretched, which is the redshift
  `ω = Nω_L/(1 + A)`, and delayed, by the extra path of the drift counted twice: the pulse reaches the electron later,
  and the wave has farther to travel back. Frequency shift and phase are one effect: the shift is the time
  derivative of the displacement.
- **At a fixed frequency, the delay is a phase.** The Fourier transform of a fixed-cycle wave train at ω has
  phase `ω·t_c(r) + const`, with `t_c` the centre of the train. So the relative phase is about
  `φ_N(r) ≈ 2 N k_L δz_mid(r)`, with `δz_mid` the drift at mid-pulse. It scales as `∝ N · a0² · n_cyc · e^{−2r²/w0²}`.
- **The axis lags most, so the lens focuses,** like glass that is thicker in the middle. It is aberrated, Gaussian
  rather than parabolic, which explains why there is no single-beam fit and where the `p ≥ 1` content comes from.
- **It does not vanish as `Z → ∞`.** The dent's depth is about a tenth of a wavelength. What matters is the phase,
  ~1 rad, not the geometric curvature compared with `Z`. The far field is the Fourier transform of `A(r)e^{iφ(r)}`,
  so the angular pattern changes, in amplitude and in phase. An apparent focus shifted by `z0` shows as the angular
  phase `k z0 θ²/2` (`θ = r/Z`). That phase doesn't depend on `Z` at fixed angle, so it stays measurable as long as
  `Z` is known. (Corrected 2026-10-30: an earlier version said `z0` becomes unmeasurable.)
- **What it changes on the screen (lattice run, Results):**
  - N = 1: the intensity core is almost unchanged, with a weak wing added. The wavefront carries `k z0 θ²/2` with
    `z0 ≈ 900 λ`, 1.7 rad at 25 mrad.
  - N = 2: the ring moves from 6.5 to 13 mrad.

  The LG fit's "waist" (17.4 λ for N = 1) is the best single-mode compromise. It is not the beam's divergence: a
  17.4 λ Gaussian at 883 λ would be 17.5 mrad wide on the screen, while the simulated core is about 13.5 mrad wide
  (12.8 mrad without dressing).
- **Cross-check:** observing along the laser direction gives `dt_obs/dξ = 1`, so neither redshift nor lens.

## Model: a per-electron dressing factor

For an electron at radius `r`, weak field, with the run's evaluation frequency `ω_e`:

$$d_N(r) = \frac{\int d\xi\, g(\xi)^N \exp\!\Big\{i\int_{\xi_0}^{\xi}\big[\tfrac{\omega_e}{\omega_L}(1 + A(r,\xi')) - N\big]\,d\xi'\Big\}}
{\int d\xi\, g(\xi)^N},$$

with `g` the envelope (flat top plus wings `e^{−Δφ²/σ²}`, as in `laser_field.cpp`) and `ξ_0` before the pulse. `ξ_0`
is the same for every electron, since they all sit in the `z = 0` plane at the waist. The coherent-source prediction
(`coherent_source_prediction` in `lg_decomposition.py`) then uses `u^N · d_N(r)` instead of `u^N`. The angular factor
and the Fresnel propagation are unchanged.

`d_N` contains both effects:

- **Phase:** the dent above.
- **Amplitude:** electrons whose `A(r)` differs from the resonant value lose amplitude in a sinc-like way. The
  width is set by the line, `~1/(N n_cyc)`.

Both grow as `N · n_cyc · a0²`.

**Which `A` is resonant matters.** `ω_e` is the harmonic dressed with `⟨a²⟩ = a0²/2`. The LG laser's on-axis field,
however, is `√2 ×` the plane wave's (`pulse_envelope_and_lg_normalization.md` §2), so the axis sees `A = a0²`. The
resonant ring is at `r ≈ 0.59 w0`, and the axis itself is detuned. The model must take `A(r)` from the code's
actual laser field, not from the nominal `a0`. For LG `p, l ≠ 0`, `A(r)` is ring-shaped and so is the dent.

**Limits:** `d_N → 1` as `a0 → 0`. For a plane wave (constant `A`), `d_N` doesn't depend on `r`.

## Open question this addresses: why N = 2 focuses more strongly

A pure phase lens has `f = k_N w0²/(4 φ_N)` near the axis. With `φ_N ∝ N` and `k_N ∝ N`, `f` is the same for every
harmonic. Applied to N = 2 (source waist `w0/√2`, same `z_R`), it predicts 12.6 λ at 980 λ. The fit gives
**7.6 λ at 542 λ** (Verification 6).

The model's amplitude part (stronger detuning apodization for N = 2) and the non-quadratic phase are the
candidates. The model includes both, so it either explains N = 2 or rules dressing out as the cause.

**Settled by test 1:** with the dressing factor, the prediction matches N = 2 to 99.8%, so the measured waist is
explained. Phase alone gives 95%; the detuning amplitude gives the rest. The thin-lens argument fails because the
phase is not parabolic: `d_2` has a phase of −3.0 rad at `w0` relative to the axis.

**The magnitude, before test 1:** the naive flat-top estimate is `φ_N(0) − φ_N(∞) ≈ N π n_cyc A_axis`. For
50 cycles and `A_axis = a0² = 0.01`, that's ~1.6 rad. Extrapolated quadratically to `r = w0`, it gives ~3 rad. The
fitted thin lens (`f ≈ 2000 λ`) corresponds to 0.98 rad at `w0`. A thin-lens fit averages a Gaussian-shaped phase
under a Gaussian amplitude, so a factor of a few is plausible, but only the full model can say. (The model gives −1.52 rad at `w0` for N = 1, and it fits.)

## Results (2026-10-30, test 1 on lattice runs)

The 2026-10-27 runs are not on the server, so the Verification 6 setup was rerun with the new polar-lattice beam:
- disk of 2 w0 = 50 λ, Gaussian LG(0,0) laser, `w0 = 25 λ`;
- backscatter screen at 2e4 λ, `R_max = 40 w0`;
- 256 × 64 screen points, `trajectory_NT 100`.

Configs are in `tests/ponderomotive_lens/`; output is in `~/output_tests/ponderomotive_lens/`. All runs used one
exclusive node (job 30462).

| Run | Beam | `a0` | Time (256 CPUs) |
|---|---|---|---|
| `20261030_092541` | lattice, spacing 2.44 λ, 1321 electrons | 0.1 | 39 s |
| `20261030_092940` | lattice, spacing 1.22 λ, 5153 electrons | 0.1 | 132 s |
| `20261030_092744` | lattice, spacing 2.44 λ, 1321 electrons | 0.01 | 39 s |
| `20261030_093318` | random, 65536 electrons | 0.1 | 1601 s |

Coherent-source prediction fidelity for `E_−`, without → with `d_N`:

| Run | N = 1 | N = 2 |
|---|---|---|
| lattice, `a0 = 0.1` | 0.7756 → **0.9999** | 0.2199 → **0.9981** |
| lattice (fine), `a0 = 0.1` | 0.7765 → 0.9999 | 0.2187 → 0.9985 |
| random, `a0 = 0.1` | 0.7781 → 0.9988 | 0.2140 → 0.9933 |
| lattice, `a0 = 0.01` | 0.9999 → 0.9999 | 0.9983 → 0.9985 |

- **The dressing factor explains the field with no free parameter.** Without it, the random run reproduces
  Verification 6 (0.778 / 0.214). With it, the lattice run reaches 0.9999 for N = 1 and 0.998 for N = 2. The random
  run is capped by its speckle: 0.07% (N = 1) and 0.56% (N = 2) of the power at other `l`, against 1e-7 and 1e-4 on
  the lattice.
- **The weak spin-flipped `E_+`** (1e-7 of the power) goes from 0.187 to 0.997 at N = 1 and from 0.035 to 0.906 at
  N = 2.
- **Values of `d_N`** (`A_axis = a0² = 0.01`, evaluated at the run's frequency):
  - `|d_1(0)| = 0.87`, `d_1(w0)/d_1(0) = 1.06 e^{−1.52i}`;
  - `|d_2(0)| = 0.56`, `d_2(w0)/d_2(0) = 1.33 e^{−3.04i}`.

  The axis is detuned because the run's frequency is resonant at `A = a0²/2`, i.e. at `r ≈ 0.59 w0`.
- **Phase vs amplitude** (lattice, `a0 = 0.1`, fidelity):

  | Factor | N = 1 | N = 2 |
  |---|---|---|
  | full `d_N` | 0.9999 | 0.9981 |
  | phase only, `d/\|d\|` | 0.9988 | 0.9491 |
  | amplitude only, `\|d\|` | 0.7840 | 0.3227 |

  The lens is the delay phase, as pictured. The detuning amplitude adds 5% for N = 2.
- **Fitted beams** (lattice, `a0 = 0.1`, `--scan-z0`, series waist):
  - N = 1: 17.4 λ at 883 λ;
  - N = 2: 7.9 λ at 579 λ.

  These match Verification 6 (17.9 λ at 914 λ; 7.6 λ at 542 λ).
- **N = 2's "coherent tail" was sampling.** At `a0 = 0.01` the lattice gives N = 2 fidelity 0.998, with 0.29% of the
  power outside the prediction's 99.9% radius. The random runs gave 0.971 and 2.7%
  (`lg_decomposition_implementation.md`, Verification 7). Speckle is ~4× larger at N = 2, since its transverse
  `q`-area on the same screen scales as `N²`.

## Tests

### 1. Dressing factor on the coherent-source prediction (done; Results above)

Add `d_N(r)` to `coherent_source_prediction` behind an option. Inputs, all from the run folder:
- `A(r)`: the code's laser normalization;
- the envelope: `laser_flat_duration`, `laser_wing_sigma`, `laser_wing_sigma_cutoff`;
- `ω_e`: `run_log.txt`.

Then compare the fidelities with the ones already measured:

| Run | N = 1 now | N = 2 now | Expected with `d_N` |
|---|---|---|---|
| Gaussian, `a0 = 0.01` (`20261027_194300`) | 0.9991 | 0.971 | unchanged (`d_N ≈ 1`) |
| Gaussian, `a0 = 0.1` (`20261027_204817`) | 0.778 | 0.214 | close to the `a0 = 0.01` values |
| LG(2,2), `a0 = 0.1` (`20261027_213655`) | 0.991 | 0.942 | not worse; ring-shaped dent |

Also run the waist fit (`--scan-z0`) on the predicted field: it should give 17.9 λ at ~914 λ for N = 1 and 7.6 λ
at ~542 λ for N = 2.

The picture is confirmed if the fidelity recovers with no free parameter. If it doesn't, the picture is
incomplete. Candidates then: the transverse ponderomotive push, or the source's thickness along z.

### 2. Single-electron calibration of `d_N(r)` (seconds per run)

Like `tests/one_electron_on_axis.cfg`:
- one electron at `r = 0, 0.25, …, 1.5 w0`;
- one screen point (on axis for N = 1);
- `a0 = 0.1` and `a0 = 0.01`.

Take the complex `F` at the run's frequency. Divide by the `a0 = 0.01` value and remove the geometric phase
`k(|X − x_e| − |X|)`. That gives `d_1(r)` from the full code, to compare with the formula. For N = 2 the dominant
helicity vanishes on axis, so use one fixed off-axis point whose angle is much larger than `r/Z`. This is the
fastest check of the formula itself, before test 1 puts it into the prediction.

With the lattice beam (`tests/ponderomotive_lens/`), each run below takes under a minute on 256 CPUs.

### 3. Evaluate at the undressed frequency (small C++ change, one lattice run)

Add a config key, for example `frequency_dressing` (default `true`). When it is `false`, `frequencies_list` is built
from `p` rather than the dressed `q` (`simulation.cpp`, `init_simulation_parameters`). Then rerun
`config_expansion_rest_a0_01.cfg` with it off.

- **Prediction:** the lens phase changes only by `ω_u/ω_e ≈ 1 + a0²/2`. `|d_N|` changes, though: the resonance moves
  to the edge (`A → 0`), so every electron is detuned and the axis most. The apparent waist and the `p` content
  shift, and the lens stays. Test 1's model gives the numbers beforehand, with `ω_e = Nω_L`.
- **What it separates:**
  - A persisting lens means a real emission delay (this picture).
  - A vanishing lens means an artefact of evaluating at a frequency some electrons don't emit at.

### 4. Scalings (lattice runs, under a minute each)

- **Pulse length:** `laser_flat_duration 25` vs `50` at `a0 = 0.1`. A dent that deepens through the pulse
  predicts `1/f ∝ n_cyc`; a static frequency mismatch would not scale with length.
- **Field strength:** `a0 = 0.05` (and 0.07) at rest, for `1/f ∝ a0²`. Only 0.01 and 0.1 exist so far.
- **Optional:** a screen at 10× the distance. `w0` and the mode content should stay; `z0` becomes poorly
  determined.

**Breakdown expected** when the phase reaches a few rad, i.e. several times the pulse length or `a0 ≫ 0.1`. Early and
late emission then come from different depths: a thick source with a spread of foci. In addition, the transverse
ponderomotive push, which grows as `t²`, starts to move electrons radially. That push is a fraction of λ for the
current runs.

## Not covered

**A moving beam** (runs A–D of the scientific report, `γ ≈ 13.5`). The same dent exists there: in the lab the
electrons slow down rather than drift. The model would need the Doppler mapping of the source, the same item still
open for the coherent-source prediction.
