# config/reference/

One reference config per beam situation, created 2026-11-23. **Don't edit them for a run:** copy one to
`config/work/` (ignored by git), change the copy, and run it:

```
cp config/reference/lattice_2d.cfg config/work/my_run.cfg
python3 py_scripts/compile_and_run.py config/work/my_run.cfg
```

All seven share one base and differ only in the beam section (and, for the spectrum, the screen and frequencies).
**A change to the shared base goes into all seven** (and into `config/config.cfg`, a copy of `random_2d.cfg`).

## Shared base

Run A of the 2026-10-29 progress report with `w0 = 10 lambda`: LG `(p, l) = (0, 2)`, `a0 = 0.1`, circular
polarization, pure Gaussian pulse (`laser_wing_sigma 25` cycles, cutoff 4, no flat top), `p_z = 13.483 mc` head-on,
harmonics N = 1, 2, `long_distance_direct_approx`, `trajectory_NT 100`, circular screen of radius `37.5 w0`
(64 × 64, equal-area rings) at `160000 w0`. Disk radius from the LG table (`beam_cylinder_radius -3`).
`laser_delay 0` since 2026-11-24 (spec 14): the pulse peak meets the beam centre at the LG focus. Before, it met it
49.9 λ downstream (2.5% below the peak intensity); the screen check below predates the change.

**Screen, checked on a `lattice_2d` run:** N = 1 fits (`w_eff = 12.4 w0`, 5e-5 of the flux within 10% of the edge);
N = 2 is clipped (`w_eff = 34 w0`, 37% near the edge). For N = 2 work, widen `circular_detector_R_max` (to ~100 w0)
and raise `circular_detector_N_R`, since the equal-area grid's first ring is at `R_max / sqrt(N_R - 1)`.

## The seven cases

| File | Beam | Trajectories |
|---|---|---|
| `random_2d.cfg` | 65536 random draws on the flat disk (`H = 0`) | 65536 |
| `lattice_2d.cfg` | Polar lattice, spacing 2.5 λ, flat disk | 284 |
| `random_3d.cfg` | 262144 random, `H = 10 λ₁`, `beam_sigma_z = 20 λ₁` (Gaussian ends) | 262144 |
| `random_3d_modulated.cfg` | As `random_3d`, modulated: `ζ = 1.11`, `Λ_mod = 20 λ₁` | 262144 |
| `lattice_stack_uniform.cfg` | Gaussian stack: `H = 10 λ₁`, `beam_stack_sigma_z = 20 λ₁`, cut 5 (`L = 210 λ₁`), 480 sheets (spacing rule: ≥ 430 for N_max = 2) | 136320 |
| `lattice_stack_modulated.cfg` | As the uniform stack, modulated as `random_3d_modulated`; 1120 sheets (aliasing rule: ≥ 1066 for N_max = 2) | 318080 |
| `spectrum_lattice_2d.cfg` | As `lattice_2d`; one screen point at `x = 7 w0` (the N = 1 ring: `w_eff/√3` for LG (0,2)), 10000 frequencies over 0.5–2.5 ω₁ | 284 |

`λ₁ = 2π/|q_axial|` (unit `lambda_1`, CLAUDE.md). The modulated cases warn that `H / Λ_mod = 0.5` is not a whole number
of periods: that only disables the cold-beam comparison in `run_log.txt`. A stack's cost is sheets × lattice sites
(`run_log.txt` prints both); it changes only through `beam_stack_sheets`, the spacing or the disk, not
`beam_particle_count`. Lowering the sheet count below the bounds above makes the solver warn (aliasing or coarse
spacing): see `theory/implemented/09_lattice_stack.md`.

Timing on a 24-thread laptop: `lattice_2d` takes 38 s; scale by trajectories.
