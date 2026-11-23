# Mode content across the emission line; purity near the vortex core

Deferred on purpose (2026-10-29). Stage reports state only that all results are at one frequency per harmonic and that
we intend to study how the mode content varies across the line.

## Motivation

Some groups want to absorb one XUV photon carrying angular momentum in an atom (or "something") and drive a transition
forbidden in first order for plane-wave photons. With the 10 μm laser and γ ≈ 13.5, N = 1 is at ~90 eV and N = 2
at ~180 eV.

## Which output is relevant: fixed frequency, not a time slice

In first-order perturbation theory the transition amplitude of an atom at x picks out the Fourier component of the
field at the transition frequency ω₀:

  `c_fi ∝ ∫ dt ⟨f|Ĥ_int(t)|i⟩ e^{iω₀t} ∝ M_fi[Ẽ(x, ω₀), ∇Ẽ(x, ω₀), …]`.

The multipole order enters through which spatial derivatives of Ẽ at the atom appear. For a resonance with line profile
L of width Γ:

  `P ∝ ∫ dω L(ω − ω₀) |M_fi[Ẽ(x, ω)]|²`,

which is quadratic, with no interference between frequencies. A superposition of fields over a band (a time slice)
never enters one-photon absorption. It would matter only for nonlinear or time-resolved processes.

Notes on time slices, for completeness: `E(x, t) = ∫ dω g(ω) F(x, ω) e^{−iωt}`, and by linearity
`a_pl(t) = ∫ dω g(ω) a_pl(ω) e^{−iωt}`. The arbitrary common phase at each ω (time origin, `e^{ikZ}`) only shifts t.
The physical content is in the ratio `a_10(ω)/a_00(ω)`: constant across the band → the same mode content at every t;
varying → spatio-temporal coupling (expected: the dressing follows the envelope). Superpose at retarded time
(multiply by `e^{−ikR}` first), because the pulse front is curved, arriving later off axis by `r²/(2Zc)`.

## Why the ω-dependence matters

- Our line width is ~1/N_cycles ≈ 1–2% (~0.9 eV at N = 1). The axis-to-edge ponderomotive shift (a0²/2 ≈ 0.5%) is
  inside it.
- **Narrow atomic lines** (bound-bound, meV or less): the atom selects one frequency somewhere in our line, which need
  not be our on-axis centre. Within the line, the dressing changes which electrons are resonant, so the mode content
  (waist, lens, p content) probably varies. All current results are one slice at the centre.
- **Core-hole / autoionizing lines** (Γ ~ 0.1–1 eV, comparable to our line): the line-weighted average of the
  quadratic quantities.
- Either way: compute the full screen at a band of frequencies around each harmonic.

## What it needs

- `dense_frequency_spectrum` currently keeps only screen point 0 (`main.cpp`). The phase recurrence needs only evenly
  spaced frequencies, so allowing the full screen is a small change. Memory and cost scale with the number of
  frequencies.
- Line resolution: the step must satisfy `δω < 2π/T_emitted` (the emitted pulse has ~N × N_cycles cycles), i.e.
  0.1–0.2% over ±3%: 30–60 frequencies per harmonic, 15–30× a two-harmonic run.
- Cheap first check: a few frequencies across the N = 1 line on a coarse screen (e.g. 32×32), and see whether
  `a_10/a_00` (or the dominant waist) moves.
- Python: run the existing decomposition per frequency, plot `|a_10/a_00|` and its phase against ω, and line-weighted
  quadratic observables.

## Purity near the vortex core

The selection-rule effects of twisted light are strongest for an atom at the phase singularity: the field vanishes there
but its gradient doesn't (trapped-ion quadrupole transition, Schmiegelow et al., Nat. Commun. 7, 12998 (2016)).
- **Purity near the core matters more than total purity.** A small admixture with a different `l` dominates the field
  where the main ring vanishes and moves the zero off centre. (Example: the 9% `l = 0` in N = 2 of
  `20261029_084449`, which turned out to be the envelope effect, pulse_envelope_and_lg_normalization.md.) Incoherent
  speckle also fills the core.
- Proposed figures of merit (Python only, no new runs): **core contrast** (intensity at the zero ÷ peak) and the
  **displacement of the zero**.
- **The LG decomposition carries through to the experiment:** the atom isn't at 40 m, the beam is refocused, and ideal
  optics map LG modes onto LG modes with the same weights. The decomposition predicts the focal spot. Different Gouy
  phases make a displaced zero rotate around the axis through the focus.
- **Target geometry:** for a localized atom the effects are large near the core. For an extended gas target averaged
  over atom positions, the orbital part largely washes out (theory work by Afanasev, Serbo, Surzhykov and co-workers;
  to be checked). Ask which kind of target the people asking have in mind.
- The equal-area circular grid is coarse near the axis, so core studies need a finer centre (an extra inner
  screen, or another radial grid).
