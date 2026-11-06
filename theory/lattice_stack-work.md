# Lattice stack: current work

Working notes for the longitudinal lattice stack and microbunching project. Three documents track it:

| Document | Content |
|---|---|
| [CoherentThomson-lattice-stack-implementation.md](CoherentThomson-lattice-stack-implementation.md) | The initial proposal (the brief). |
| This file | The current work: analysis, decisions, open questions. |
| [lattice-stack-implemented.md](lattice-stack-implemented.md) | What was actually done, and how it was checked. |

## Status (2026-11-06)

- Working by discussion: each part of the brief is discussed before it is implemented.
- Done: the trajectory memory fix, commit `3052872` on branch `lattice-stack`.
- Phase 1 is implemented and its checkpoint passes (`lattice-stack-implemented.md`). It waits for the user's review and is
  not committed yet. Phase 2 does not start until the user says so.

## Decisions so far

| # | Question | Decision |
|---|---|---|
| 1 | `a0` for validation | **0.1** (run A's value) |
| 2 | Sheet normalization | **`w = 1/M`** per sheet, so a stack carries the same total charge as the flat disk |
| 3 | Memory fix first, as its own commit | **Yes**, done |
| 4 | Which `q` sets the modulation period | **The analytic `q = k_s - k_L`**, see below |
| 5 | `q` in the form factor `b` | **The exact `k_s − (K_s·P / K_L·P) k_L`** (2026-11-06) |
| 6 | Correct report Sec. 6 | **Yes**, Markdown done; PDF still to rebuild, as this laptop has no pandoc |

### The analytic `q`

The "diagnostic `q`" and the analytic `q` are the same quantity. "Diagnostic" only names the routine that prints it:
[form_factor.cpp:33-39](../src/core/simulation/form_factor.cpp#L33-L39) evaluates `q = k_s - k_L` directly from the
formula, with nothing measured. Here `k_s` is the dressed fundamental along the detector normal (`<a²> = a0²/2`), and
`k_L` is the laser wavevector. So:

- `Λ = 2π/|q| = 2π/|k_s - k_L|`, and the modulation period follows from it directly.
- For radiation harmonic `N`, on axis and head-on, `|q_N| = N (k_1 + k_L)`, because harmonic `N` absorbs `N` laser
  photons. **Report Sec. 6.1 has `N k_1 + k_L`, which is wrong for `N ≥ 2`** (see "The exact `q`" below).
- The brief's sheet-translation test is the only thing that *measures* a slope. It becomes a one-off check that the
  solver follows the analytic value, not an input. The prediction is `2k/(1+β₀)`, which differs from the analytic `q`
  by 7e-6 at `a0 = 0.1` (9e-4 rad across 280 nm), far too little to matter.
- The formula should exist once in the code, shared by the log and the stack generator, rather than being copied.

Still to settle for Phase 3: whether the config gives `Λ_mod` as a length or in units of the computed `Λ`. Units of
`Λ` would keep the phase matching exact when `a0` changes (`Λ` = 13.77 nm at `a0 = 0.1`, 13.97 nm at `0.2`).

### The exact `q` (measured 2026-11-06)

A recorded open point (`theory/dev_notes.md`, `coherent_incoherent_decomposition.md`) asked whether the phase should use
`k_s − k_L` or `k_s − (K_s·P / K_L·P) k_L`. Translating the beam and fitting the field phase settles it: the field gets
`e^{-i q·r}` with `q_N = N [k_1 − (K_1·P / K_L·P) k_L]`, undressed `P`. On axis and head-on, this is `2 N k_1 / (1+β₀)`.

Head-on, slope over the screen, in `k_L` (`q_formula_head_on.py`):

| `a0` | N | Measured | Dressed form | `N (k_1 + k_L)` | `N k_1 + k_L` |
|---|---|---|---|---|---|
| 0.1 | 1 | 726.531123 | 726.531125 | 726.536100 | 726.536100 |
| 0.1 | 2 | 1453.062208 | 1453.062250 | 1453.072200 | 1452.072200 |
| 0.2 | 1 | 715.846822 | 715.846892 | 715.866500 | 715.866500 |
| 0.2 | 2 | 1431.693636 | 1431.693784 | 1431.733000 | 1430.733000 |

Tilted, small plane-wave disk (R = 0.1 λ_L, 20 electrons), one on-axis point, N = 1 (`q_formula_tilted_small_disk.py`):

| θ | `a0` | Axis | Measured | Dressed form | `k_s − k_L` |
|---|---|---|---|---|---|
| 0.9π | 0.2 | x | −0.302958 | −0.302958 | −0.309017 |
| 0.9π | 0.2 | z | 698.329299 | 698.329309 | 698.347957 |
| 0.9π | 0.1 | x | −0.307480 | −0.307480 | −0.309017 |
| 0.9π | 0.1 | z | 708.751768 | 708.751725 | 708.756457 |
| 0.75π | 0.2 | x | −0.693245 | −0.693245 | −0.707107 |
| 0.75π | 0.2 | z | 611.016322 | 611.016345 | 611.030207 |

Consequences:

- **The stack period:** `Λ = 2π/|k_s − k_L|` (decision 4) stays, as the difference is negligible head-on: 9e-4 rad across
  280 nm at `a0 = 0.1`.
- **Tilted wide disks are unreliable:** a 46.5 λ_L disk at θ = 0.9π misses the x slope by 2–4% (`q_formula_tilted_wide_disk.py`,
  and also the screen fit). The small disk does not, which points to the fixed integration window (`tau_0_traj`, known bug).
  Head-on over a stack this does not apply, and the Phase 1 scan tests that directly.
- **Report numbers:** with `N (k_1 + k_L)` in `compare_bunch_length.py`, the Sec. 6.3 table changes at N = 2 only:

| `q₁σ` | `q₂σ` before → after | ρ̃ before → after | `\|c\|/ρ̃` before → after |
|---|---|---|---|
| 0.25 | 0.4996 → 0.5000 | 0.88265 → 0.88250 | 0.99676 → 0.99694 |
| 0.5 | 0.9993 → 1.0000 | 0.60697 → 0.60654 | 0.98885 → 0.98955 |
| 1 | 1.9986 → 2.0000 | 0.13573 → 0.13535 | 0.93546 → 0.93808 |

## Open items

- **Report PDF** must be rebuilt where pandoc exists (`theory/scientific_report/README.md`, section PDF). The Markdown is corrected.
- **The period in Phase 3:** decision 4 predates decision 5. The log now lists both values, and they differ by 7e-6 at `a0 = 0.1`.
- **CLAUDE.md edits proposed on 2026-11-06, not applied:** the lattice constraint should also name `beam_sigma_z 0`;
  `beam_lattice_spacing` should be listed as optional; the `plot_processed.py` usage needs correcting; and a pointer to
  the brief should be added.
- **Next:** the user's review of Phase 1, then the discussion of Phase 2.
- **Phase 3b plan, as discussed:** one modulation period (about 275 nm, 296,000 trajectories at 1 nm), at 1 nm and 0.5 nm on the
  small screen; a full screen once if the mode decomposition is wanted; more periods only for an explicit finite-train test.
- **Later for Phase 2:** a standalone Python tool for `b` at any `q`, `N` and ζ from the sheet weights, checked against the log.
- **Phase 3 keys:** `H` need not be a whole number of periods (log `H/Λ_mod`, warn); sheets per period should be whole, with `d`
  adjusted and logged.

## Analysis of the brief (as given to the user on 2026-11-06)

Copied as written, apart from link paths, heading levels and one removed paragraph about a file that no longer exists. Questions 1 to 4 at the end are now answered above, and
the memory fix is done.

I've read the brief and the code it names, and for Schaap I searched the paper's cached text for the key equations
rather than reading it in full. Nothing is implemented, and the only files I wrote are two memory notes and a couple of
scratch scripts. The plan holds up, but the memory problem below has to be fixed before any production-size run.

### Blocker: memory

`Electron` keeps its whole trajectory and nothing ever frees it
([electron.cpp:22](../src/core/particle/electron.cpp#L22), line number before the fix). Each state is about 104 bytes,
and a run has 20,000 states per electron (200 cycles × `trajectory_NT` 100), so about 2 MB per electron:

| Run (1,076 sites per sheet) | Electrons | Memory kept |
|---|---|---|
| Flat disk | 1,076 | 2 GB |
| 280 sheets (1 nm spacing) | 301,280 | **~630 GB** |
| 560 sheets (0.5 nm spacing) | 602,560 | **~1.25 TB** |

I'd free each trajectory right after its electron's radiation is computed, as a separate commit checked for
byte-identical output. Scaling from the report's 37 min for 65,536 electrons, a full screen with 301,280 electrons would
take about 2.8 h on 256 cores.

### Part by part

**Shared conventions**

- **The slope test needs no new code:** `beam_center_z` already translates a sheet, so we can run it with the current
  binary.
- **Expected slope:** from the plane-wave translation symmetry, I get a field ratio of `e^{-iq z}` with
  `q = 2k/(1+β₀)`, where β₀ is the initial velocity. This is my derivation, not a measurement.
- **`q` itself is safe:** the slope differs from the diagnostic `q = k_s + k_L` by only 7e-6 (`a0` = 0.1) to 3e-5
  (`a0` = 0.2). That's at most about 0.003 rad across 280 nm, so the diagnostic `q` can set the modulation period.
- **The sign does matter:** `form_factor.cpp` uses `e^{+iq·r}`, so the complex field ratio should match the conjugate of
  `b`. This makes no difference for a symmetric uniform stack, but it does for modulated ones.

**Phase 1: uniform stack**

- **Generator:** it's a small extension, and `M = 1` should stay identical because the z position and random draw order
  don't change. The existing check already rejects non-zero `beam_cylinder_height` and `beam_sigma_z`, and I'd keep that.
- **Weights:** these go in the three accumulation points of [radiation.cpp](../src/core/radiation/radiation.cpp#L396)
  (lines 396, 560 and 844), plus the mean momentum and `A` in `form_factor.cpp`.
- **Normalisation:** I'd use `w = 1/M` per sheet, so a stack carries the same total charge as the flat disk and the
  fields compare directly.
- **Also needs updating:** `plot_beam_density_z.py` assumes a cylinder, and `plot_beam_scatter` should be off for runs
  this size.
- **Main risk:** every electron starts at `t = 0`, so the integration window doesn't shift with the sheet. By my
  estimate, the endpoint terms then pick up phase with slope `k` rather than `q`, and the pulse's sharp 4σ edge is
  clipped slightly differently per sheet. This could set the noise floor at the `H = Λ` minimum, which is what the
  brief's absolute-error reporting would show.

**Phase 2: density model**

- **Formula confirmed:** the direct integral matches `(-1)^n J_n(nζ)` to 6 digits at (n, ζ) = (1, 1), (20, 1) and
  (20, 1.2). The formula also matches Eq. (9) in the paper.
- **Binning error is larger than the brief's estimate:** the density has a cusp at ζ = 1 and folds above it, so the error
  is about twice the 0.8% figure and doesn't always shrink steadily:

| n | ζ | Bins per period | Error in `b` | Error after sinc correction |
|---|---|---|---|---|
| 1 | 1 | 14 / 28 | 1.65% / 0.37% | 0.82% / 0.16% |
| 20 | 1 | 280 / 560 | 1.63% / 0.37% | 0.80% / 0.16% |
| 20 | 1.111 | 280 / 560 / 1120 | 2.1% / 1.2% / 0.38% | 1.3% / 1.4% / 0.43% |
| 1 | 1.5 | 14 / 28 / 56 | 0.62% / 0.66% / 0.15% | 0.22% / 0.45% / 0.10% |

These are top-hat bins from 40 million quadrature points. So a 1 nm versus 0.5 nm comparison won't prove convergence when
ζ > 1, and I'd add a third, finer grid to the diagnostic.

**Phase 3a: `b₁`**

- Bin count must be an integer per period, and the period `Λ` depends on `a0` (13.77 nm at 0.1, 13.97 nm at 0.2). The
  actual spacing therefore comes from the chosen config, as the brief expects.
- The same comparison needs a uniform stack of identical length and charge, which is a cheap control.

**Phase 3b: `b₂₀`**

- The values are |J₂₀(20)| = 0.165 at ζ = 1, with a maximum of 0.243 at ζ = 1.111.
- At 280 nm per period, the cost is the same as the table above for one period, so several periods at full-screen
  resolution will be expensive. This phase is the one that needs a cost decision.

### Questions on the shared part and Phase 1

1. **Which `a0` should validation use?** I recommend run A's 0.1, because its lens is weaker and easier to separate from
   stack effects. Part II used run B's 0.2, so choose that if you want direct continuity with Sec. 6.
2. **Is `w = 1/M` the normalisation you want?** The alternative is one unit of charge per electron, but that makes the
   field grow with `M`, so the fixed-charge comparison would need rescaling.
3. **Should the memory fix go first, as its own commit?** I recommend yes.
4. **Do you agree that the diagnostic `q` sets the modulation period,** provided the measured slope confirms the
   prediction above?

## Supporting facts from the code

- `compute_radiation` runs `compute_trajectory` for every electron, and `main.cpp` has already run it on the first 10
  for `electron.dat`, so those 10 are integrated twice. The radiation loops use the first `N_tau` states, so this looks
  harmless. I checked the exact path only, and it costs a little extra time and memory.
- Lattice settings of the report's run B config (`theory/scientific_report/configs/runs_A-D_lattice/`): `trajectory_NT 100`,
  `radiation_formula long_distance_direct_approx`, `laser_wing_sigma 25 cycles`, `laser_wing_sigma_cutoff 4`,
  `laser_flat_duration 0`, `average_pz 13.483 mc`, zero momentum spread, `beam_lattice_spacing 2.5 lambda`.
- Paper's constants (Schaap et al. 2023): `ζ = L/f`, `f = cγ³β³/[Δγ(ω₁−ω₂)]`, `b_n = (−1)^n J_n(nζ)` (Eq. 9), maximum of
  `J_n` at `ζ ≈ 1 + 0.81 n^(−2/3)` for `n ≥ 4`, proposal `b₂₀ ≈ 10%` at 13.6 nm from a 272 nm modulation.

## Scripts behind the numbers

Run with `.venv/bin/python -I <script>` (needs numpy and scipy).

### Bessel values, phase slope and memory

```python
import numpy as np
from scipy.special import jv
z=np.linspace(.001,2,200001)
for n in (1,20):
    v=abs(jv(n,n*z)); i=v.argmax(); print('n',n,'max |J_n(n z)|=%.4f at zeta=%.4f; at zeta=1: %.4f'%(v[i],z[i],abs(jv(n,n))))
M=2000000;p=2*np.pi*np.arange(M)/M
for n,ze in [(1,1.0),(20,1.0),(20,1.2)]:
    print('n',n,'zeta',ze,'direct b=',np.round(np.mean(np.exp(-1j*n*(p+ze*np.sin(p)))),6),'(-1)^n J=%.6f'%((-1)**n*jv(n,n*ze)))
pz=13.483;g=np.hypot(1,pz);b=pz/g;k0=(1+b)/(1-b);print('undressed w1/wL',k0)
for a0,w in [(0.1,725.5),(0.2,714.9)]:
    qd=w+1;qt=2*w/(1+b);dq=qd-qt
    print('a0',a0,'q_diag',qd,'q_transl',qt,'rel diff %.2e'%(dq/qd),'phase over 280 nm %.2e rad'%(dq*2*np.pi*.28/10),'Lambda nm %.3f'%(1e4/qd))
N=20000;S=13*8
for sheets,sites in [(1,1076),(280,1076),(560,1076),(14*5,1076)]:
    print(sheets,'sheets',sheets*sites,'electrons -> %.0f GB'%(sheets*sites*N*S/1e9))
```

Output:

```
n 1 max |J_n(n z)|=0.5819 at zeta=1.8412; at zeta=1: 0.4401
n 20 max |J_n(n z)|=0.2433 at zeta=1.1110; at zeta=1: 0.1647
n 1 zeta 1.0 direct b= (-0.440051-0j) (-1)^n J=-0.440051
n 20 zeta 1.0 direct b= (0.164748+0j) (-1)^n J=0.164748
n 20 zeta 1.2 direct b= (0.161913+0j) (-1)^n J=0.161913
undressed w1/wL 729.1637845659956
a0 0.1 q_diag 726.5 q_transl 726.4949753613063 rel diff 6.92e-06 phase over 280 nm 8.84e-04 rad Lambda nm 13.765
a0 0.2 q_diag 715.9 q_transl 715.8804381609893 rel diff 2.73e-05 phase over 280 nm 3.44e-03 rad Lambda nm 13.968
1 sheets 1076 electrons -> 2 GB
280 sheets 301280 electrons -> 627 GB
560 sheets 602560 electrons -> 1253 GB
70 sheets 75320 electrons -> 157 GB
```

The frequencies 725.5 and 714.9 are the report's rounded values (Sec. 2.2), so the slope differences are estimates.

### Binning error

```python
import numpy as np
from scipy.special import jv
# Map uniform initial phase through phi=phi0+zeta sin(phi0), top-hat bin onto Nb sheet centres, then
# compare the sheet-weight Fourier sum at harmonic n with (-1)^n J_n(n zeta) and with its sinc-corrected value.
Q=40_000_000
p0=2*np.pi*(np.arange(Q)+.5)/Q
print(' n zeta  bins  |  exact        sheet sum     err       err after sinc correction  min/max bin weight rel. mean')
for n,zeta,bins_list in [(1,1.0,(14,28,56)),(1,0.5,(14,28)),(1,1.5,(14,28,56)),(20,1.0,(280,560,1120)),(20,1.111,(280,560,1120))]:
    phi=np.mod(p0+zeta*np.sin(p0),2*np.pi)
    for Nb in bins_list:
        h=np.bincount(np.minimum((phi/(2*np.pi)*Nb).astype(int),Nb-1),minlength=Nb)/Q
        c=2*np.pi*(np.arange(Nb)+.5)/Nb
        b=np.sum(h*np.exp(-1j*n*c))
        ex=(-1)**n*jv(n,n*zeta)
        s=np.sinc(n/Nb)   # sin(pi x)/(pi x), x=n/Nb: bin-average factor
        print(f'{n:3d} {zeta:5.3f} {Nb:5d} | {ex:+.6f}  {b.real:+.6f}{b.imag:+.1e}j  {abs(b-ex)/abs(ex):7.3%}   {abs(b/s-ex)/abs(ex):7.3%}   {h.min()*Nb:.3f}/{h.max()*Nb:.2f}')
```

Output:

```
 n zeta  bins  |  exact        sheet sum     err       err after sinc correction  min/max bin weight rel. mean
  1 1.000    14 | -0.440051  -0.432808+9.0e-07j   1.646%    0.816%   0.502/3.21
  1 1.000    28 | -0.440051  -0.438420+4.5e-07j   0.371%    0.161%   0.501/5.03
  1 1.000    56 | -0.440051  -0.439679+2.3e-07j   0.084%    0.032%   0.500/7.91
  1 0.500    14 | -0.242268  -0.240223-8.7e-18j   0.844%    0.007%   0.670/1.81
  1 0.500    28 | -0.242268  -0.241760-4.3e-17j   0.210%    0.000%   0.667/1.94
  1 1.500    14 | -0.557937  -0.554482-6.0e-17j   0.619%    0.220%   0.401/4.19
  1 1.500    28 | -0.557937  -0.554262-6.9e-17j   0.659%    0.450%   0.400/4.86
  1 1.500    56 | -0.557937  -0.557076-6.9e-17j   0.154%    0.102%   0.400/6.30
 20 1.000   280 | +0.164748  +0.162057-9.0e-07j   1.633%    0.803%   0.500/22.94
 20 1.000   560 | +0.164748  +0.164141-4.5e-07j   0.368%    0.159%   0.500/36.36
 20 1.000  1120 | +0.164748  +0.164610-2.3e-07j   0.084%    0.031%   0.500/57.66
 20 1.111   280 | +0.243302  +0.238118+5.7e-16j   2.131%    1.304%   0.474/22.27
 20 1.111   560 | +0.243302  +0.246158-7.2e-16j   1.174%    1.387%   0.474/41.56
 20 1.111  1120 | +0.243302  +0.244216+9.7e-17j   0.376%    0.428%   0.474/55.83
```
