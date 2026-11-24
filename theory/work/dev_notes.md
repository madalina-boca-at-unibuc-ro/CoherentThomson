# Development notes

New entries, newest last. The record up to 2026-11-23 is frozen in `dev_notes_2026-07-28_to_2026-11-23.md` (paths
there predate the reorganization; `theory/archive/reorganization_plan.md` maps them).

## 2026-11-24: collision timing (spec 14)

- LG runs collided off the focus: the pulse peak met the beam at z = β/(1+β)·cutoff·σ·λ_L (49.9 λ_L for cutoff 4,
  74.9 for 6) while the focus is at the origin. With w0 = 10 λ_L that is 2.5–5.4% below the peak intensity, and N = 1
  at ω₁ (the red edge of the line at a0 = 1) changed 1.8× between cutoff 4 and 6. Found while chasing a fake N = 2
  signal; record: `~/Dropbox/work/AI-outputs/Claude/2026-11-24/screen_distance_and_harmonics_findings.md`.
- Spec 14 implemented to its checkpoint: `laser_delay` honoured, each electron starts at the pulse front, b(q) and
  the beam diagnostics use t = 0 positions. At the focus N = 1 no longer depends on the cutoff (ratio 1.000037).
  `theory/implemented/14_collision_timing.md`.
- Side result (n2_residue tests): the near-axis N = 2 of `simplified` at a0 = 1 was an O(dτ²) residue of the LR/BR
  cancellation (J_z = 0); use `long_distance_direct*` and cutoff 6 for weak harmonics.
