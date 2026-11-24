# 14: Collision timing (implemented)

Spec: `theory/specs/14_collision_timing.md`.

> **⚠ Collision geometry.** The start time, the free-flight pull-back and the collision point are derived for a
> laser phase φ = k(ct − n̂·x) with the focus at the origin, electrons in free flight before the pulse, and the
> head-on baseline. If the collision geometry changes (co-propagating or crossed beams, a moved focus, an envelope
> that is not a function of φ alone), re-derive spec §3 and re-run every check, not only the related ones.
> 1 − n̂·β → 0 for co-propagation makes the start run away.

- **Status:** done 2026-11-24, except the Python cross-check (spec 7.7): the Python reference still has the old
  arrangement (electrons start at t = 0, pulse front at them), so it can only be compared in the old geometry
  (`laser_delay = cutoff·σ`, which `config_cross_check.cfg` keeps). Deferred by the user.
- **Commit:** see `git log` (branch `spec14-collision-timing`).
- **Configs:** the seven reference configs and `config.cfg` have `laser_delay 0` (collision at the focus) and new
  comments for `laser_delay` and the heatmap; `config/reference/README.md` notes it. Report, check and cross-check
  configs keep `laser_delay = cutoff·σ` (old geometry); the reports' LG runs use w0 = 25 λ_L, where the old
  off-focus collision costs 0.06% of the peak intensity (user's decision: no note).

## Code

- `LaserField` (`laser_field.cpp`): the `delay = wing_sigma_cutoff * wing_sigma` overwrite is gone; `delay` is
  `laser_delay`, the phase at which the plateau starts. New `get_delay()`.
- `Particle::free_flight_crossing(laser, phi, x, beta)` (`electron_factory.{hpp,cpp}`): ct = (φ/k + n̂·x)/(1 − n̂·β)
  and the point on x + β ct.
- `place_electron` takes the laser instead of `tau_0`: each electron starts at `free_flight_crossing(φ_min)` with its
  own velocity, x⁰ = c t₀, τ₀ = t₀/γ. `generate_electron` and `generate_cylinder_beam` lose `tau_0`.
- `Electron::get_reference_position()`: the configured position at t = 0 (x⁰ = 0), a new constructor argument
  defaulting to the initial position. Used by `compute_form_factor` (`b(q)`), `axial_positions` and
  `plot_electron_beam_scatter`.
- `simulation_parameters::tau_0_traj` removed.
- `main.cpp` heatmap: t = delay / ω.
- `run_log.txt`: `laser_delay` (phase and cycles) in the Laser section, and a new section
  `Collision (beam centre, free flight)`: the pulse front (start), the peak or plateau start, the plateau end if
  any, each with ct, point, distance from the origin, distance along n and, for LG, that over z_R with the on-axis
  peak intensity; then the beam's earliest and latest start ct.
- CLAUDE.md: the "Known bugs, revisit together" bullet replaced by a collision-timing bullet and the
  collision-geometry warning (the heatmap-window note kept). `AGENTS.md` (untracked copy) not updated.

## Deviations and judgment calls

- **Old-geometry agreement is not to ~1e-12 as the spec expected.** For θ = π the laser's n̂ is
  (sin π, 0, −1) with sin π ≈ 1.2e-16 (`create_unit_light_like_vector` does not snap; the polarization vectors do).
  So n̂·x ≈ 1e-16·x and off-axis electrons start at ct ≈ ±2.6e-10 a.u. (1.4e-15 λ) instead of 0. Each term (LR, BR)
  then changes by ~3e-9 (stepped-phase LD mode), and their cancellation (to 5% at N = 1, 6e-4 at N = 2) magnifies
  it in the total. Not snapped: snapping n̂ would also change the old laser phase.
- The plane-wave check (spec 7.3) is not "a global phase to rounding": moving the collision by 49.93 λ_L moves the
  source relative to the screen. Corrected for the path difference e^{+ik Δz cos θ}, the residual drops to the 1/R
  level (|c| − 1 = 2e-5 = Δz/D).

## Checks (configs, scripts, logs: `~/Dropbox/work/AI-outputs/Claude/2026-11-24/spec14_checks/`; runs in
`~/output_tests/spec14/`)

Old binary: `03141ed` (dirty only by the untracked AGENTS.md), saved before the change.

| # | Check | Result |
|---|---|---|
| 1 | `checks/thomson_dipole`, old vs new | `radiation_field.dat`, `electron.dat` byte-identical; dipole check \|E\|/\|E_dip\| ∈ [1, 1.000058] |
| 1 | 2D LG lattice (a0 = 1, 284 e, 16×16, simplified LD approx), `laser_delay = cutoff·σ` | beam scatter and b(q) identical; each term 3e-9; total 3.8e-6 (N = 1), 4.4e-5 (N = 2), see above |
| 2 | `random_3d` (512 e, 16×16), old vs new | beam scatter and b(q) byte-identical; total 6.3e-8 (N = 1), 6.5e-7 (N = 2) |
| 3 | plane wave, `laser_delay` 100 vs 0 cycles | \|c\| = 1 − 2.0e-5; residual 2.3e-4 → 2.9e-5 (N = 1) and 2.9e-3 → 4.7e-4 (N = 2) after the path correction |
| 4 | LG at the focus (`laser_delay 0`), a0 = 1, direct, 4418 e, 32×32, cutoff 4 vs 6 | N = 1 ratio 1.000037 (was 1.80); collision printed at 0 along n; N = 2: 7.4e-21 → 6.9e-24 |

Check 4 also settles the n2_residue side finding (`~/Dropbox/work/AI-outputs/Claude/2026-11-24/`): N = 1 at ω₁ is
2.89e-17 at the focus vs 1.79e-17 (cutoff 4) and 9.95e-18 (cutoff 6) off focus, the red-edge sensitivity to the
peak a reached.

## Checks 5–6 (same folder; lattice-stack outputs in its `lattice_stack/`)

| # | Check | Result |
|---|---|---|
| 5 | simplified vs direct (exact formulas) at the focus: 2D LG lattice, a0 = 1, 284 e, 16×16 | N = 1: 4.35e-2 (NT 100) → 1.08e-2 (NT 200), 4.0× per doubling. N = 2: 1.0e3 → 2.5e2 relative to direct, also 4×: at a0 = 1 near axis the physical N = 2 is ~0 and simplified's O(dτ²) LR/BR residue dominates (n2_residue finding) |
| 5 | rotation covariance: laser, beam momentum and detector rotated 0.3 rad about y, at the focus, `long_distance_direct_approx` | screen rotated to 2e-7 of D; N = 1: E, B → R E, R B with one factor −1.000000 (the π-rotated laser frame `rotation_matrix_from_direction` picks), residual 8.0e-7; N = 2 (the cutoff-4 leakage, ~1e-4 of N = 1): residual 4.0e-3 |
| 5 | one-electron configs (`checks/one_electron/`), old vs new | `radiation_field.dat`, `electron.dat` byte-identical |
| 5 | `checks/nt_scan` | matches the table in `work/future_dev/numerical_parameters_and_cost.md` §1 (e.g. p_z = 13.48: N = 2 simplified 6.28e-4 at NT 100, direct 5.1e-7) |
| 6 | lattice stack, `~/output_tests/lattice_stack/` deleted first | `phase1_stack_scan`: as `phase1_results.txt` to ~1e-9. `gaussian_stack_checks` (old = pre-spec-14): beam scatter, laser and incident files byte-identical, field and trajectory not (start offsets); checks 2–4 as `gaussian_stack_results.txt` to the 1e-8 noise. `phase3a_checkpoint`: the screen results as before to ~1e-8; **the single on-axis point at N = 2 improves** from |c − pred| = 7e-4 … 1.2e-2 to 5e-5 … 3e-4 (every sheet now covers the same part of the pulse). `form_factor_exact_q`, `q_formula_*`, `phase1_translation_slope`: translation relation unchanged (dressed q to ≤1e-7 head-on) |
| 6 | `rework_checks` | not run: it needs the pre-rework binary (`35b3a1a`, accepts `beam_stack_spacing`); the pre-spec-14 binary refuses that key. Not affected by spec 14 |

The tracked `checks/lattice_stack/phase3a_results.txt` was restored after the run (it is the spec-09 checkpoint's
record); the new output is `phase3a_results_new.txt` next to the other outputs.
