# Collision timing: honour `laser_delay`, start each electron at the pulse front

> **Status: APPROVED 2026-11-24.** Written 2026-11-24. Stop at the checkpoint (§8).

> **⚠ Collision geometry warning (carry into CLAUDE.md and the `implemented/` note).** Everything in this spec (the
> start time, the free-flight pull-back, the trajectory length, which positions `b(q)` uses, the collision point
> printed in `run_log.txt`) is derived for the present geometry: a laser whose phase is φ = k(ct − n̂·x) with its
> focus at the origin, electrons in free flight before the pulse, and the head-on baseline (n̂ = −ẑ, beam along +ẑ).
> **If the collision geometry is ever changed** (laser tilt beyond the existing `laser_direction_*` use, a
> co-propagating or crossed beam, a focus away from the origin, a pulse model whose envelope is not a function of
> φ alone, fields present at the start), **re-derive §3 and re-run every check in §7**, not only the ones that look
> related. In particular 1 − n̂·β (§3) becomes small for co-propagation: the pull-back and the lab duration of the
> interaction then grow without bound.

## 1. Problem

Origin: `~/Dropbox/work/AI-outputs/Claude/2026-11-24/screen_distance_and_harmonics_findings.md`, §4 side finding.

- `LaserField`'s constructor overwrites the parsed delay: `delay = wing_sigma_cutoff * wing_sigma`
  (`laser_field.cpp:34`). The `laser_delay` key has no effect.
- `tau_0_traj = 0.0` is hardcoded (`simulation.cpp:48`), and every electron starts at x⁰ = c·τ₀ = 0 at its
  configured position (`electron_factory.cpp:28`).
- Together: at t = 0 the pulse plateau start (the peak for `laser_flat_duration 0`) is `cutoff·σ·λ_L` ahead of
  the beam, which was the intent. But the LG focus is at the origin (`laser_field.cpp:212`), where the electrons
  start, while they meet the peak at

  z_c = β/(1 + β) · cutoff · σ · λ_L   (49.9 λ_L for cutoff 4, σ = 25 cycles; 74.9 λ_L for cutoff 6).

- Effect: with w0 = 10 λ_L (z_R = 314 λ_L) the peak intensity at the collision is 97.5% (cutoff 4) or 94.6%
  (cutoff 6) of the focus value, while ω_N is dressed with the focus a0. At a0 = 1, N = 1 at ω₁ sits on the red edge
  of the line and changed 1.8× between cutoff 4 and 6 with both radiation formulas (runs in
  `~/output_tests/n2_residue/`).
- No config-only fix: moving the beam back (`beam_center_z = −β·cutoff·σ·λ_L`) puts the collision at the focus, but
  the trajectory length is fixed at `2·cutoff·σ + flat` cycles of phase from the start, so the integration then ends
  at the peak.

## 2. New meaning of the keys

- **`laser_delay`** (cycles; may be 0 or negative): the time t_d at which the plateau start crosses the origin. The
  envelope function is unchanged (plateau over φ ∈ [delay, delay + flat]), so for `laser_flat_duration 0` this is
  the pulse peak. Window: φ ∈ [φ_min, φ_max] = [delay − cutoff·σ, delay + flat + cutoff·σ], as the existing
  `phi_min`/`phi_max` formulas already give once the overwrite is removed.
  - **Judgment call:** the reference stays the plateau start, not its centre. This keeps every existing config's
    meaning (§5). For a flat top centred at the focus, use `laser_delay = −flat/2`.
- **Beam positions** (`beam_center_*`, the cylinder, lattice, stack, modulation): the electrons' **free-flight
  positions at t = 0**, i.e. where each electron would be at t = 0 without the field. With `laser_delay 0` and
  `beam_center_z 0`, the beam centre meets the peak at the focus.

## 3. Start of each electron

For electron j with configured position x_j and initial velocity v_j = c β_j (its own momentum, including any
spread), its phase in free flight is φ_j(t) = k[ct(1 − n̂·β_j) − n̂·x_j]. It starts at the pulse front, φ_j(t₀) =
φ_min:

  c t₀,j = (φ_min/k + n̂·x_j) / (1 − n̂·β_j),   x_j(t₀) = x_j + v_j t₀,j,   τ₀,j = t₀,j / γ_j.

- Head-on (n̂ = −ẑ): c t₀ = (φ_min/k − z_j)/(1 + β). With `laser_delay 0`, flat 0, z_j = 0:
  t₀ = −cutoff·σ·T_L/(1 + β), pull-back β/(1 + β)·cutoff·σ·λ_L (49.93 λ_L), and the peak is then at
  +cutoff·σ·λ_L/(1 + β), exactly cutoff·σ·λ_L from the electron.
- **Trajectory length unchanged:** `N_tau` and `d_tau` as now (`simulation.cpp`, `simulation_length` from
  φ_max − φ_min, `d_tau` scaled by the mean (p⁰ − n̂·p)/mc). Each electron covers the whole window from its own
  start; no lengthening.
- The pull-back is valid because the envelope at φ < φ_min is below e^{−cutoff²} (≈ 10⁻⁷ at cutoff 4).
- Plumbing: `generate_cylinder_beam` / `place_electron` already take `tau_0`. Replace the single `tau_0_traj` by
  the per-electron τ₀,j and x⁰ = c t₀,j; `simulation_parameters::tau_0_traj` goes, or becomes unused.

## 4. Positions that must stay at t = 0

The trajectory's start position is no longer the configured position (head-on, the start z is
z_j/(1 + β) + const: an axial structure is compressed by 1 + β ≈ 2). These users of `get_initial_position()` mean
the configured (t = 0) positions and must read a new stored value, e.g. `Electron::get_reference_position()`:

- `form_factor.cpp:45` (`b(q)`: q·r at a common time; the translation relation holds for positions at equal t).
- `electron_factory.cpp:192` (`axial_positions`, modulation and stack diagnostics).
- `electron_plotter.cpp:60` (`electron_beam_scatter.dat`, the configured beam).

`electron.dat` keeps the true trajectory (it starts at x(t₀), τ₀ < 0 possible); `plot_trajectory.py` subtracts the
first row, so it is unaffected.

## 5. Backward compatibility

Every tracked config has `laser_delay = cutoff·σ` (as the old comment required). With §2 that gives φ_min = 0 and,
for electrons with n̂·x_j = 0 (any flat disk at z = 0, head-on), t₀ = 0: **the same start as now.** Electrons with
n̂·x_j ≠ 0 (3D beams, stacks, a tilted laser) now start at their own front instead of all at t = 0. They follow
the same worldlines (the field before the front is negligible), sampled on a shifted τ grid.

So existing configs reproduce the old collision geometry, including its off-focus collision for LG. Opting in to
the focus is one edit: `laser_delay 0`. The seven reference configs and `config.cfg` switch to it (§9); report,
check and cross-check configs keep `cutoff·σ`.

## 6. Other changes

- `main.cpp` heatmap: t = delay/ω (plateau start at the origin), replacing the recomputed `cutoff·σ/ω` and its
  comment. Same value for old configs.
- `radiation_plotter.cpp:178` (`incident_field.dat`): uses `phi_mid` of the window, which follows `delay`
  automatically. Check only.
- `lg_decomposition.py` dressing: uses the envelope shape only. Unaffected.
- `run_log.txt`, Laser section: the `laser_delay` used and, for the beam centre, the collision time and point
  (where φ = delay, the plateau start, on its free-flight line; for a flat top also where φ = delay + flat), its
  distance from the focus, and for LG that distance over z_R.
  Beam section: the range of t₀,j.
- Reference configs (after the checkpoint, §8): `laser_delay 0` and a new comment in the seven reference configs
  and `config.cfg`, one identical edit (shared base). Re-measure the costs in `config/reference/README.md` only if
  they move (they should not: same trajectory length). Frozen report configs are not edited.
- CLAUDE.md: replace the "Known bugs, revisit together" bullet (keep the heatmap-window note) and add the collision-geometry
  warning at the top of this spec; `theory/reference/implementation_details.md` if it describes the start.

## 7. Checks

1. **Old geometry reproduced:** `checks/thomson_dipole` and a 2D head-on LG config with `laser_delay = cutoff·σ`:
   outputs agree with the pre-change binary (byte-identical if φ_min rounds to exactly 0; otherwise to rounding,
   ~1e-12 relative).
2. **3D beam, old geometry:** `random_3d` shrunk, `laser_delay = cutoff·σ`: agrees with the old binary to
   O(dτ²) (shifted τ grid). `b(q)` and `electron_beam_scatter.dat` byte-identical (they use t = 0 positions).
3. **Plane wave is translation invariant:** plane-wave config, `laser_delay` = cutoff·σ vs 0: |F| equal, phase
   differs by one global factor (relative phase between screen points unchanged, to rounding).
4. **Collision at the focus:** LG, `laser_delay 0`, cutoff 4 vs 6 (direct formula, the a0 = 1 `n2_residue`
   setup): N = 1 at ω₁ now independent of the cutoff (the 1.8× gap should close to the e^{−cutoff²} level).
   `run_log.txt` reports the collision point at z ≈ 0.
5. The standard set, since the start of every trajectory moved: `simplified` vs `direct` (O(dτ²)), rotation
   covariance (tilted laser: n̂·x_j ≠ 0 exercises the general t₀), `checks/nt_scan`, one-electron configs.
6. Delete `~/output_tests/lattice_stack/` and re-run the lattice-stack checks (CLAUDE.md: they reuse finished runs).
7. Python cross-check with `config_cross_check.cfg` (unchanged config, old geometry): agreement as before.

## 8. Checkpoint

Stop after §3–§6 and checks 1–4; report before §7.5–7.7 and before any config's `laser_delay` value is changed.

## 9. Decisions (user, 2026-11-24)

- **Reference configs switch to `laser_delay 0`** (collision at the focus), after the checkpoint.
- **Reports are left as they are.** They were preliminary, and their LG runs all use w0 = 25 λ_L (z_R = 1963 λ_L):
  the off-focus collision (≤ 50 λ_L for cutoff 4, σ = 25) costs 0.06% of the peak intensity. No note needed.
- **`laser_delay` refers to the plateau start**, as proposed in §2 (keeps the flat-top configs' meaning).
