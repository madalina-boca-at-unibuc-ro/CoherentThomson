# Spec: beam-averaged effective a0 for the ponderomotive dressing

> **Status: DRAFT, NOT APPROVED FOR IMPLEMENTATION.**
> Do not implement anything from this spec, and do not edit any code, config or doc because of it, until the user
> explicitly says to start. When you are given this file: read it, read the code it names, and reply with your
> analysis only — open questions, problems with the plan, anything that disagrees with the code — then stop and
> wait. The user wants to discuss it first. Remove this notice only when the user approves.

Target code: `~/Dropbox/work/bin/c++/codes/CoherentThomson` (C++). Read its `CLAUDE.md` first and follow its
conventions (class + factory, all config parsing in factories/IoUtils, `clang-format`, Release builds, no test
suite: verification by the checks at the end).

## 1. Why

The radiation frequencies and the modulation length unit `lambda_1` are built from the ponderomotively dressed
electron momentum, `PhysUtils::dressed_momentum(P, k_L, xi)`, with `<a^2> = xi^2/2` and `xi = laser.get_a0()`
(the global laser `a0`, the same for every electron).

For a focused LG laser each electron sees a different local amplitude, `a(r) = a0 * M(r)`, where `M` is the code's
own transverse mode magnitude (normalization `Npn` included, so `M` is not 1 on the peak: e.g. `sqrt(2)` on axis for
LG_0^0). Electrons where the field is strong dominate the emission, so the line sits at the dressing those electrons
feel, not at the global `a0`.

This spec adds an **optional** single effective amplitude `xi_eff`, averaged over the electron beam's transverse
distribution and weighted by each electron's fundamental emission amplitude. It replaces `a0` **only in the
dressing**. The default stays the current behaviour.

## 2. Definition

Assumptions (deliberate, do not generalize):
- **Constant electron density** over the beam's transverse disk (radius `R` = the resolved `beam_cylinder_radius`,
  including the negative-radius LG-table option). Use the configured geometry, **never the generated electrons**:
  the result must not depend on `random_seed`, on `beam_particle_count`, or on lattice/stack settings.
- **Fundamental only.** Weight = `M` (the harmonic-1 amplitude is proportional to the local `a`). One `xi_eff` is
  used for **all** harmonics, so `omega_N = N * omega_1` stays exact and the frequency recurrence in `radiation.cpp`
  is untouched. No per-harmonic dressing.
- **Waist plane.** Evaluate `M` at `z_loc = 0` (waist, `w(z) = w0`). Longitudinal position, bunch height `H`,
  `beam_sigma_z` and modulation play no role.
- Magnitude only: the helical phase `exp(i l phi)` does not enter.

With local dressing `<a^2>(r) = xi(r)^2/2 = a0^2 M(r)^2/2` (the code's convention, valid for every polarization;
see the comment above `dressed_momentum` in `phys_utils.hpp`):

```
xi_eff^2 = a0^2 * ( ∫_disk M^3 dA ) / ( ∫_disk M dA )
```

and `dressed_momentum(P, k_L, xi_eff)` is used wherever `dressed_momentum(P, k_L, laser.get_a0())` is used now.

For the LG mode at the waist (from `LaguerreGaussLaser::complex_amplitude`, envelope and carrier removed):

```
M(rho) = Npn * (sqrt(2) rho / w0)^|l| * |1F1(-p, |l|+1, 2 rho^2 / w0^2)| * exp(-rho^2 / w0^2),
Npn    = sqrt(2) / |l|! * sqrt((p+|l|)! / p!)
```

For `PlaneWaveLaser`, `M = 1`, so `xi_eff = a0` exactly.

## 3. Config

One new optional key (absent = `global`):

```
dressing_a0   global          # global (default): xi = laser_a0, current behaviour
                              # beam_averaged: xi_eff of this spec
```

Any other value throws. Add it to the list of optional keys in `CLAUDE.md` ("Config file format").

## 4. Implementation

1. **Mode magnitude.** Add a public `virtual double transverse_mode_magnitude(double x_loc, double y_loc, double
   z_loc) const` to `LaserField` (plane wave: returns 1; LG: `M` above at general `z_loc`, with `w(z)` and the
   `w0/w` factor, exactly as in `complex_amplitude`). Share the expressions with `complex_amplitude` rather than
   writing an independent copy, so the two cannot drift apart.
2. **The average.** One function, e.g. `double effective_dressing_a0(const ConfigMap&, const LaserField&)` (place it
   where `init_simulation_parameters` and the beam factory can both reach it; `IoUtils` or `PhysUtils`, your call),
   returning `laser.get_a0()` for `global` and `xi_eff` for `beam_averaged`.
   - Integrate over the disk in the beam's own transverse plane: canonical point `(x, y, 0)` plus
     `beam_center_x/y`, rotated with the same `beam_axis_rotation` the generator uses, then projected on the laser
     frame: `x_loc = r . epsilon_1`, `y_loc = r . epsilon_2`, and `z_loc = 0` forced (waist plane). This handles an
     off-axis or tilted geometry; for the coaxial baseline it reduces to the 1D radial integral.
   - Deterministic quadrature: Gauss–Legendre in `r` on `[0, R]` times uniform trapezoid in azimuth (spectrally
     accurate for a periodic integrand), weight `r dr dphi`. Refine until doubling both counts changes `xi_eff` by
     less than `1e-10` relative; log the counts used.
   - `R = 0` (thin line beam): fall back to `M` on the axis point and say so in the log.
3. **One source of truth.** All three current uses of the dressing must take the value from that function, so the
   frequencies and the modulation stay matched:
   - `src/core/simulation/simulation.cpp:99` (frequencies),
   - `src/core/particle/electron_factory.cpp:82` (`axial_translation`, hence `lambda_1`),
   - `src/core/simulation/form_factor.cpp:36` (`b(q)` diagnostic).
   Compute it once per run and pass it in; don't recompute it at each call site.
4. **Run log.** In `run_log.txt` add: the `dressing_a0` mode, `a0`, `xi_eff`, `xi_eff/a0`, `<a^2>` both ways, the
   relative frequency change `omega_1(xi_eff)/omega_1(a0) - 1`, the corresponding `lambda_1` change, and the
   quadrature counts. Per `CLAUDE.md`, the log re-reads config keys the way the factories do.
5. **Python scripts.** `grep` `py_scripts/` for anything that rebuilds frequencies, `lambda_1` or the dressing from
   `laser_a0` (e.g. theory columns in `plot_observables.py`, `lg_decomposition.py`). Any such script must use the
   value the run actually used (read it from `run_log.txt` or write it into `incident_config.cfg`), not
   recompute it from `laser_a0`. List what you changed or found.
6. **Docs.** Update `CLAUDE.md` ("Frequencies" bullet: the dressing amplitude is `laser_a0` or `xi_eff`) and add a
   dated entry to `theory/dev_notes.md`.

Out of scope: per-harmonic dressing, using the generated electrons, averaging along `z`, any change to the
trajectories (the motion always uses the true local field).

## 5. Verification

1. **No-op by default.** Without the key (and with `dressing_a0 global`), the `.dat` outputs are byte-identical to
   the pre-change build (same config, same `random_seed`), for `config/config.cfg` and a lattice-stack config.
2. **Plane wave.** `beam_averaged` with `laser_type plane_wave` gives `xi_eff/a0 = 1` to rounding and
   byte-identical outputs to `global`.
3. **Reference values** (coaxial beam, uniform disk, waist plane), from an independent 1D integral:

   | LG (p, l) | R / w0 | xi_eff / a0 |
   |---|---|---|
   | (0, 0) | 2.37 | 0.817985 |
   | (0, 1) | 2.37 | 0.623694 |
   | (0, 2) | 2.37 | 0.550987 |
   | (1, 1) | 2.37 | 0.501581 |
   | (0, 0) | 1.86 | 0.829631 |

   The code must reproduce these to `1e-6` relative. Also independent of `random_seed`,
   `beam_particle_count`, `beam_distribution` and `beam_cylinder_height`.
4. **Mode magnitude consistency.** At a few points (on axis, on the ring, at `2 w0`; `z_loc = 0` and
   `z_loc = z_R`) `transverse_mode_magnitude` equals `|amplitude| / E0_c` from `complex_amplitude` at a phase on
   the flat top, to rounding. One-off check; record it in `dev_notes.md`.
5. **Off-axis sanity.** Shifting `beam_center_x` by `+d` and `-d` gives the same `xi_eff`, and the value moves in the
   expected direction for LG_0^0 (lower).
6. **Consistency of the three uses.** With `beam_averaged`, the log's `lambda_1`, the first entry of
   `frequencies_list` and the form factor's dressed momentum all correspond to the same `xi_eff`.
7. **Physics check (for the user, not required for acceptance).** `dense_frequency_spectrum true` at a screen point
   in the coherent cone, scanning about `±1e-2` in `omega/omega_1` with step `<= 5e-4`, once with `global` and once
   with `beam_averaged`: the coherent peak should sit closer to `omega/omega_1 = 1` with `beam_averaged`.

Report: what changed (files), the log lines of one `beam_averaged` run, and the results of checks 1–6.

## 6. Review comments (2026-11-10, Claude; spec postponed, still not approved)

Checked against the code: the reference table reproduces to 6 digits (independent scipy 1D integral); `M(rho)` matches
`LaguerreGaussLaser::complex_amplitude` (`src/core/laser/laser_field.cpp`); `beam_center` is added before the beam-axis
rotation, as §4.2 assumes; the three named call sites are right. No circular dependency: `xi_eff` needs the beam
radius, which never takes the `lambda_1` unit.

Open points, to settle before implementing:

1. **Temporal envelope missing, and as large as the spatial effect.** The argument of §1 applies in time too. With
   `a = a0 M(r) g(phi)` it factorizes: `xi_eff^2 = a0^2 * S_space * S_time`, `S_time = ∫g^3 / ∫g`. For the pure
   Gaussian pulse of the current configs (`laser_flat_duration 0`) `S_time = 1/sqrt(3) ≈ 0.577`, against
   `S_space = 0.818^2 ≈ 0.67` for LG_0^0 at `R = 2.37 w0` (flat top: `S_time ≈ 1`). `lg_decomposition.py`'s
   `dressing_factor` already integrates `g^2` over the pulse. Without it (or an explicit, reasoned exclusion),
   check 7 is ambiguous. Decide this first: it changes the definition of `xi_eff`.
2. **Two more dressing uses, not listed in §4.3:** `generate_cylinder_beam` calls `resolve_beam_lengths` itself
   (`electron_factory.cpp`, so the beam factory's signature must take `xi_eff`), and `run_log.cpp` calls
   `resolve_beam_lengths` again to log `lambda_1`; it would silently log the global-`a0` value. Five sites in all.
3. **Switching the mode moves the beam.** Lengths in `lambda_1` (height, `beam_sigma_z`, `beam_stack_sigma_z`,
   modulation period) change with `xi_eff`, so check 7's `global` vs `beam_averaged` also compares two different
   beams. Use absolute units there, or state it.
4. **Weight.** `M` (amplitude) is right for the coherent peak; a random beam's incoherent background goes as `M^2` and
   would give a different `xi_eff`. Say the spec targets the coherent line.
5. **Signature.** `effective_dressing_a0(const ConfigMap&, ...)` would parse the radius (negative-radius LG table) and
   centre a second time; take the parsed `CylinderBeamParams` instead.
6. **Waist plane** holds only for `H`, `beam_center_z` small against `z_R`; log `H/z_R` or warn.
7. **Checks.** Check 4's "phase on the flat top" doesn't exist for Gaussian-only pulses: use the envelope peak.
   Check 2's bit-exactness is safe (`sqrt(a0^2 * 1) == a0` in IEEE).
8. **Python.** `lg_decomposition.py` builds its per-radius dressing from `laser_a0` at the evaluated frequency, so it
   is independent of this choice and needs no change; it could predict check 7's peak shift. `plot_observables.py`
   does not rebuild the dressing.
