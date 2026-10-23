# Split of the propagation phase: e^{ik|X|} · e^{ik(R − |X|)} (all detectors, all modes)

Revision of `stabilize_phase_axis_aligned_detector.md`. **Status: implemented 2026-10-23** (summary in
`implementation_details.md` §4.5).

## 1. Why

Each electron's radiation carries the propagation phase k R_j, with R_j = |X − r_j| from the electron to a
screen point X. In double precision, k R_j is rounded at the level k R_j · 10⁻¹⁶, separately for every electron
(and, in the exact modes, every τ step). Those errors are what decorrelate the electrons, so they directly
degrade the coherent sum.

The error grows fast with γ. The emitted k grows like γ², and keeping the Fresnel number N_F = a²/(Dλ_s) fixed
for a fixed source size a requires D ∝ γ², so kD ∝ γ⁴. Scaling from the 2026-10-22 run (γ ≈ 10, D = 2·10⁶ λ,
k₁ = 0.0165 a.u.⁻¹, kD ≈ 5·10⁹ rad):

| γ   | kD (rad) | phase rounding per (electron, screen point) |
|-----|----------|---------------------------------------------|
| 10  | 5·10⁹    | ~10⁻⁶ rad                                   |
| 30  | 4·10¹¹   | ~10⁻⁴ rad                                   |
| 100 | 5·10¹³   | ~10⁻² rad                                   |
| 300 | 4·10¹⁵   | ~1 rad: coherence destroyed                 |

Higher harmonics multiply k further. The distance sweep in §5 confirms the onset.

## 2. The split

For each screen point X, exactly,

$$
e^{ikR_j}=e^{ik|\mathbf X|}\,e^{ik\,\Delta R_j},\qquad
\Delta R_j=R_j-|\mathbf X|=\frac{|\mathbf r_j|^2-2\,\mathbf X\cdot\mathbf r_j}{R_j+|\mathbf X|}\approx-\mathbf n\cdot\mathbf r_j .
$$

- e^{ik|X|} is computed **once per (frequency, screen point)**, in `run_simulation`. Its rounding is shared by
  every electron at that point, so it cancels from the relative phases between electrons.
- ΔR_j is computed without cancellation and depends only on the electron's position, whatever the screen's
  distance, angle or orientation (planar, tilted, backward, spherical: no special cases, no direction guard).

This is an exact identity, not a far-field approximation. **Everything per-electron is unchanged:** each
electron's own n₀,j, its own amplitude 1/R_j, its own transverse basis, and its own τ integral over its own
(LG-driven) trajectory. The output keeps its meaning: `radiation_field.dat` still contains the full e^{ikR_j}.

**Correction to the first version of this proposal.** It subtracted one screen-wide D (the screen distance,
with ΔR = (s² − 2X·r + |r|²)/(R + D), s² the in-plane offset). That ΔR contains the screen-point term
ρ²/(2D) ≈ Dθ²/2, which grows with D at a fixed angular window; at D = 8·10¹² w0 the exact `direct` mode was still
3% off at h3 off-axis. Subtracting each point's own |X| removes that term, and needs no s² at all.

## 3. Code

- `Detector_2D::add_point` stores each lab-frame point together with `point_distances[i] = |X_i|` (computed from
  the stored point, so that |X|² and X·r refer to the same rounded point); `get_point_distance(i)`.
- `radiation.cpp`: `distance_minus_reference(X, r, |X|, R)`.
  - Exact modes: `radiation_phase_argument(x, R − |X|)` per τ (still frequency-independent).
  - Long-distance modes (both unstepped and stepped): reconstruction constant `|x_0| − |X| + r⁰(τ_m)`.
- `simulation.cpp`: after the thread reduction, each (frequency, point) is multiplied by
  `general_factor · e^{ik|X|}`.

## 4. Unchanged

Trapezoid weights, the boundary term, the stepped-phase recurrence and re-anchoring, the two transverse
projections, n₀,j, 1/R_j, the reference sample r(τ_m), the output convention.

## 5. Validation (2026-10-23, laptop, Release)

Base `config/config.cfg` (γ ≈ 10, a0 = 0.2, LG), 20 electrons, N_harmonics = 3. Old = binary built from the code
just before the change.

**Regression, old vs. new**, max |ΔF| / max |F| per harmonic, all 6 modes × 5 geometries (circular forward at
D = 2·10⁶ λ; circular backward; rectangular tilted θ_d = 0.02π, φ_d = 0.3π at 2000 λ; spherical; beam shifted to
x = 3 λ, z = 50 λ with height 10 λ):

- Long-distance modes: ~10⁻⁶ at the far circular screen, ~10⁻⁹ at the 2000 λ screen (scales with kD: the old
  rounding), up to ~3·10⁻⁵ at h3 of the simplified forms (weak-harmonic cancellation).
- Exact modes: up to 2·10⁻² at h3 on the far screen. The new result is the correct one: against
  `long_distance_direct` (whose own error is O(1/D)), exact `direct` at h3 went from 2.1·10⁻² to 8·10⁻⁵ (forward),
  3.6·10⁻³ to 9.9·10⁻⁴ (backward), 2.7·10⁻⁵ to 8·10⁻⁶ (shifted beam). The remaining backward 10⁻³ is the
  long-distance approximation's: it is independent of `trajectory_NT` and falls 100× per 100× in D.
- Exact `simplified` vs. long-distance differs equally before and after (its known quadrature error).

**Distance sweep**, D = 8·10⁴ … 8·10¹² w0 (kD = 5·10⁹ … 5·10¹⁷), angular window fixed (R_max = 2.5·10⁻⁴ D),
one on-axis point (the 20 electrons coherent) and two off-axis (incoherent), D²·Σ|F|²:

- Old: off-axis points drift from 8·10⁸ w0; at 8·10¹⁰ the long-distance result is off by 20–60%, and at
  8·10¹² exact `direct` is ~10⁴× too small everywhere.
- New: every mode converged to 6–7 digits from 8·10⁸ to 8·10¹² w0, including exact `direct` at h3 off-axis
  (1.635843 vs. long-distance 1.635842, ×10⁻²⁸).

**Speed:** exact modes ~40% faster (0.53 → 0.32 s per electron, 256 points, measured twice), presumably because
sin/cos of small arguments avoid the slow argument reduction; long-distance modes unchanged.
