# 16: Screen fields at arbitrary frequencies (implemented)

Spec: `theory/specs/16_omega_list.md`.

- **Status:** done 2026-11-24. The user asked for the full implementation and a commit on `main`, so the checkpoint
  (§7) was a self-check, not a stop.
- **Commit:** see `git log` (on `main`).

## Code

- `IoUtils::get_omega_list` (`io_utils.hpp`): optional `omega_list`, `N:value` entries then one unit
  (`first_harmonic_frequency` or `omega_laser`); throws on a missing or unknown unit, a malformed entry, N < 1, a value
  ≤ 0 or trailing characters, naming the entry.
- `Simulation::init_simulation_parameters`: with `omega_list` (and `dense_frequency_spectrum false`) the frequencies
  are the values × `fundamental_frequency` or × ω_L/c, in the given order; repeated values warn; with the dense mode it
  throws. New `simulation_parameters::harmonics` (N per frequency; 0 in the dense mode) and `frequency_mode`.
- `run_log.txt`, "Frequency spectrum": `frequency mode` (with, for `omega_list`, that `N_harmonics*` are ignored and
  whether the list is evenly spaced), and `N = n` on each frequency line; harmonic and `omega_list` runs list every
  frequency (the dense mode keeps its range summary above 20).
- `plot_observables.read_harmonic_indices(run_dir)` and `harmonic_index(...)`: N from those lines, else round(ω/ω₁)
  (dense runs, older logs). Used by `add_theoretical_oam_energy_ratio` (and its caption) and `lg_decomposition.py`
  (J = N(l + σ), its `harmonic` filter; the "not a harmonic" warning now fires only when N was guessed).
  `plot_point_spectrum_observables.py` keeps round(ω/ω₁): it runs on dense spectra, where no N is logged.
- `plot_spectrum_map.py`: per harmonic window the energy-weighted centroid and the local maxima of the integrated
  spectrum above 10% of the window maximum (a sample above both neighbours, flat tops once), and one ready
  `omega_list` line of all of them.
- No change to the radiation code: uneven lists already fall back (no recurrence; `_approx` → unstepped).
- CLAUDE.md (Frequencies, optional keys), `implementation_details.md` §4.1 (cost).

## Judgment calls

- Malformed `omega_list` throws (the spec), unlike the general warn-on-unknown-unit rule.
- Duplicates only warn (spec), from `init_simulation_parameters` (called once), not the parser (called several times).

## Checks (configs and logs: `~/Dropbox/work/AI-outputs/Claude/2026-11-24/spec16_checks/`; runs in
`~/output_tests/spec16/`)

| # | Check | Result |
|---|---|---|
| 1 | no `omega_list` (`lattice_2d`, 16×16), previous binary vs new | `radiation_field.dat` byte-identical |
| 2 | `omega_list 1:1 2:2` vs the harmonic run | identical (max rel diff 0) |
| 3 | uneven `1:1.0 1:1.03 2:1.97` on a 3-point line vs a dense line (stepped `_approx`, grid containing them) | 2.9e-13, 9.7e-7, 1.3e-7: the stepped phase vs the unstepped fallback |
| 4 | 8 frequencies, even vs uneven, `lattice_2d` 16×16 | `long_distance_direct_approx` 4.2 → 29.8 s (×7.2); `simplified` 16.8 → 46.0 s (×2.7) |
| 5 | a0 = 1 screen run (`lattice_2d`, 32×32, cutoff 6) with the 12-frequency line printed from `~/output/20261124_175505` | `plot_observables.py`: N = 1 Λ_zz/P_z = theory to 5–6 digits at all 7 N = 1 peaks; N = 2 rows use N = 2 (round would give 3), off theory because the screen clips N = 2 (13% beyond 0.9 R). `lg_decomposition.py`: J = −3 and −6, no "not a harmonic" warning. Older and dense runs: reader returns None, behaviour as before |
| 6 | 7 malformed lists and `omega_list` with a dense run | all refused, exit 1, entry named |
| 7 | peaks on the a0 = 1 line run | match the integrated spectrum (1.2715, 1.3305, 1.373, 1.4305, 1.4665, 1.493, …); the printed line ran as check 5 (262 s) |

## Possible follow-up

An uneven list could keep the stepped speed-up by calling the stepped function once per frequency (a single
frequency is trivially evenly spaced): N_freq trajectory passes instead of the unstepped path. Not done.
