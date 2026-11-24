# theory/: index

Every document, with its status. Status lives here, not in folder names, so a status change never breaks a link.
Layout since 2026-11-23 (`archive/reorganization_plan.md`): `physics/` formulas and derivations, `reference/` living
descriptions of the code, `specs/NN_*.md` task specs and `implemented/NN_*.md` what was done for each (same number),
`work/` logs and open investigations, `archive/` superseded or one-off. Reusable checks are in `../checks/`, reports
and the experiments behind them in `../reports/`.

## Open issues

Detail in `work/dev_notes_2026-07-28_to_2026-11-23.md` unless noted.

- **Old runs need their observables re-run:** OAM numbers before 2026-10-29 (centered-difference `∂φ`, ~5% low at
  N = 2) and absolute spectral densities / photon numbers before 2026-11-21 (`2π` too small; ratios unaffected,
  `work/spectral_density_2pi_fix.md`). Re-run `plot_observables.py` / `plot_point_spectrum_observables.py` before
  quoting them.
- `b(q)` (`Simulation::compute_form_factor`, spec 06): harmonic 1 along the detector normal,
  `q = k_s − (K_s·P / K_L·P) k_L` with undressed mean `P`; a translated electron's field gets `e^{-i q·r}` and
  `q_N = N q_1` (not `N k_1 + k_L`).
- `LaguerreGaussLaser` matches the standard LG mode numerically but isn't physics-reviewed (Allen 1992; Siegman
  ch. 17); `PhysUtils::non_linear_Thomson_formula` is flagged in-code for review.
- The flat-top envelope leaks N = 1 into N = 2 (`laser_flat_duration 0` removes it); the LG on-axis field is √2 × the
  plane wave's at the same `a0`: `work/future_dev/pulse_envelope_and_lg_normalization.md`. Unacted plans:
  `work/future_dev/README.md`.
- **Python cross-check in the old geometry only** (spec 14): the Python reference still starts the electrons at t = 0
  with the pulse front at them, so `config_cross_check.cfg` keeps `laser_delay = cutoff·σ`. Port spec 14 there before
  comparing collisions at the focus.
- Ponderomotive dressing acts as a lens, so `N ≥ 2` isn't a single LG mode: `work/future_dev/ponderomotive_lens.md`.
- **Lattice stack (spec 09):** Phases 1–3a done; the beam rework (2026-11-09) and Gaussian stack (2026-11-10) await
  the user's review; Phase 3b next. State and "Resume here": `work/lattice_stack-work.md`.

## physics/

| File | Content | Status |
|---|---|---|
| `FT_faraday_tensor_v1.md` | FT of the Faraday tensor, direct and simplified forms: notation, boundary term, Jacobian | Superseded by v2 for the simplified short-range formula; still cited by the code |
| `FT_faraday_tensor_v2.md` | The same, current simplified form | Current |
| `angular_momentum_numerics.md` | Numerical AM density and flux | Current |
| `scaling_laws.md` | Dimensionless parameters for the flat-disk results | Draft |

## reference/

| File | Content | Status |
|---|---|---|
| `implementation_details.md` | Radiation modes, optimizations, precision, timings | Current; keep in sync with `radiation.cpp` |
| `lg_decomposition_reference.md` | `lg_decomposition.py` options and outputs | Current |

## specs/ and implemented/

Next free number: 16.

| NN | Spec | Implemented | Status |
|---|---|---|---|
| 01 | Long-distance direct and simplified formulas | Reading of the coding guide | Done |
| 02 | Gaussian-smoothed electron cylinder | Stub | Done |
| 03 | LG beam radius table (+ `.png`) | Stub | Done |
| 04 | Radial phase removal (post-processing stage 1) | Stub | Done |
| 05 | Propagation phase split `e^{ik|X|} · e^{ik(R − |X|)}` | Stub | Done |
| 06 | Single-wavevector form-factor diagnostic | Stub | Done (with the measured `q`, see stub) |
| 07 | LG decomposition of the screen field | Implementation notes | Done |
| 08 | Gaussian waist from the detector Poynting flux | Stub | Done |
| 09 | Lattice stacks and ponderomotive microbunching | What was implemented | Phases 1–3a done; rework and Gaussian stack await review; 3b next |
| 10 | Beam-averaged effective a0 | — | Draft |
| 11 | Scaled parameters of a run | — | Draft |
| 12 | Paper run campaign (`reports/2026_thomson_paper/`) | — | Draft; depends on 11 |
| 13 | Helicity-component plots `E₊`, `E₋` | — | Draft |
| 14 | Collision timing: honour `laser_delay`, start each electron at the pulse front (collision at the LG focus) | Implementation notes | Done (Python cross-check deferred: the reference keeps the old arrangement) |
| 15 | Dense spectrum on a screen line: dP_z/dω, \|E₊\|, \|E₋\| maps in (θ, ω) | Implementation notes | Done |

## work/

| File | Content | Status |
|---|---|---|
| `dev_notes.md` | Development log | Current (new entries) |
| `dev_notes_2026-07-28_to_2026-11-23.md` | Investigation history up to the reorganization | Frozen |
| `lattice_stack-work.md` | Spec 09's work log, "Resume here" | Current |
| `lg_mode_decomposition_plan.md` | Staged plan behind specs 04, 07, 08 | Done stages; kept for the plan |
| `spectral_density_2pi_fix.md` | The `2π` fix in the spectral densities (2026-11-21) | Done |
| `future_dev/README.md` | Findings and plans not yet acted on (index of the five below) | Current |
| `future_dev/pulse_envelope_and_lg_normalization.md` | Flat-top leak, LG normalization | Open |
| `future_dev/numerical_parameters_and_cost.md` | Numerical parameters, cost, minimal settings | Open |
| `future_dev/coherent_incoherent_decomposition.md` | Coherent and incoherent parts without every electron | Open |
| `future_dev/spectral_line_and_vortex_core.md` | Mode content across the line; purity near the core | Open |
| `future_dev/ponderomotive_lens.md` | The dressing as a dented source | Open (tests 2–4) |

## Other

| File | Content | Status |
|---|---|---|
| `cluster.md` | Slurm commands, limits, long-run recipe; head-node report | Current |
| `archive/reorganization_plan.md` | The 2026-11-23 moves, old path → new path | Done |
