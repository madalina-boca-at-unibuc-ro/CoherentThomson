# Screen fields at arbitrary frequencies: `omega_list`, harmonic index in the run log, spectrum peaks

> **Status: APPROVED 2026-11-24.** Written 2026-11-24. The user asked for the whole implementation and a commit on
> `main`, so the checkpoint (§7) is a self-check, not a stop.

## 1. Purpose

A normal run images the screen at ω = N ω₁ only (`frequencies_list`, the solid lines of the spectrum maps). At large
a0 that misses the emission: harmonic N fills the band N × [(1 + a0²/2)/(1 + g a0²/2), 1 + a0²/2] ω₁ (g = the laser's
peak |a|²/a0², 0.541 for the code's LG (0, 2)), and ω₁ is below it for a0 ≳ 0.3; at a0 = 1 the energy within ±1% of
ω₁ is 7e-11 of the line (`theory/implemented/15_spectrum_map.md`, a0 scan). Wanted: screen runs at frequencies the
user picks, guided by the line spectrum.

Workflow: line spectrum (spec 15) → `plot_spectrum_map.py` writes the peaks and a ready `omega_list` line into
`run_log.txt` → paste it into a screen config → normal run → `plot_field.py`, `plot_observables.py`,
`lg_decomposition.py` as usual.

## 2. Config: `omega_list`

Optional, read only when `dense_frequency_spectrum false`; without it the harmonics stay as now (byte-identical).

```
omega_list    1:1.47 1:1.27 2:2.51   first_harmonic_frequency
```

- Entries `N:value`: the harmonic N (integer ≥ 1) the frequency belongs to, and the frequency. N is needed because it
  can't be inferred: the scripts' N = round(ω/ω₁) fails for broad lines (at a0 = 1 the N = 2 band is 2.36–3.0,
  so ω = 2.9 would read as N = 3), and N sets the angular-momentum prediction J_z = N(l + σ).
- Unit (one, last token): `first_harmonic_frequency` (ω/ω₁, the unit of the maps and `run_log.txt`: ω₁ is the
  default-dressed fundamental, unchanged) or `omega_laser`. No unit, an unknown unit, a malformed entry, N < 1 or a
  value ≤ 0: **throw** (unlike the general warn-on-unknown-unit rule: a wrong frequency silently images the wrong
  thing).
- Kept in the given order (i_omega = position in the list). Duplicate entries: a warning.
- `omega_list` with `dense_frequency_spectrum true`: throw (two frequency modes at once).
- `N_harmonics` / `N_harmonics_min` stay required (as today) but are ignored when `omega_list` is set; the log
  says so.

Parsing: `IoUtils::get_omega_list(config)` → optional vector of `{N, value, unit}`, next to `get_omega_range`.
`Simulation::init_simulation_parameters` converts with the same scaling as the dense mode
(`first_harmonic_frequency`: × `fundamental_frequency`; `omega_laser`: × ω_L / c) and also returns the harmonic index
of every frequency (`simulation_parameters::harmonics`): N_harmonics_min + i in the harmonic mode, the given N with
`omega_list`, and 0 (= none) in the dense mode.

**Solver:** nothing else. Every radiation mode already accepts a frequency list that is not evenly spaced: the
phase recurrence is skipped (`are_evenly_spaced`), and `_approx` falls back to the unstepped long-distance function.
The cost per frequency rises; it is measured (§6.4) and noted in `theory/reference/implementation_details.md`.

## 3. Run log

"Frequency spectrum" section:
- `frequency mode    harmonics | omega_list | dense` (and, for `omega_list`, that `N_harmonics*` are ignored, and
  whether the list is evenly spaced, i.e. whether the fast recurrence applies);
- the per-frequency lines gain the harmonic: `[i] <omega/c>  (<x> omega_1)  N = <n>` for the harmonic and
  `omega_list` modes (the dense mode keeps its range summary).

## 4. Scripts

- A shared reader, `read_harmonic_indices(run_dir)` (in `plot_observables.py`, next to
  `_read_fundamental_frequency_au`): i_omega → N from the run log's lines; **None** for a dense run or a run log
  without them (older runs), in which case the callers keep today's N = round(ω/ω₁) with its warning.
- Users of the guess, switched to the reader: `plot_observables.add_theoretical_oam_energy_ratio` (L_N/ω_N),
  `lg_decomposition.py` (J = N(l + σ), its `harmonic` filter and its "not a harmonic" warning, which now fires
  only when N was guessed). Grep for other `round(omega` uses when implementing.
- Titles and file names that say "ω index"/"harmonic" stay per i_omega; where a title shows ω/ω₁ it adds N.

## 5. Spectrum peaks in `plot_spectrum_map.py`

Per harmonic window (the bands of spec 15's summary), added to the run-log section:
- the local maxima of the line-integrated spectrum dW/dω **above 10% of the window's maximum**, as ω/ω₁, sorted by
  height, each with its height relative to the window maximum;
- the energy-weighted centroid of the window;
- one ready line covering every harmonic in the scan:
  `omega_list    1:1.4665 1:1.2715 ... 2:2.511   first_harmonic_frequency`.

A local maximum is a sample larger than both neighbours (on the scan's own grid, no smoothing); plateaus count once.
Caveat in the docstring: a peak of the integrated spectrum is where the most power is in total, not where a given
island of the map is; pick others from the map.

## 6. Checks

1. **No `omega_list`:** a harmonic run (`lattice_2d.cfg` shrunk) byte-identical to the previous binary, apart from
   the run log's new lines.
2. **List = harmonics:** `omega_list 1:1 2:2 first_harmonic_frequency` against the harmonic run: equal to rounding
   (same frequencies; evenly spaced, so the recurrence path).
3. **Uneven list = dense line:** a short uneven list on a 3-point line against a dense-spectrum run on the same line
   whose grid contains those frequencies: equal to rounding (the non-recurrence path vs the recurrence path).
4. **Cost:** time per frequency, evenly spaced vs not, for `long_distance_direct_approx` and `simplified`, on one
   config; recorded.
5. **Scripts:** on an `omega_list` run at a0 = 1 with ω = 2.9 tagged N = 2, `plot_observables.py` and
   `lg_decomposition.py` use N = 2 (J = 2(l + σ)) and don't warn; on an old run they behave as before.
6. **Refusals:** each malformed `omega_list` of §2 throws with a message naming the entry.
7. **Peaks:** on the a0 = 1 line run (`~/output/20261124_175505`) the listed peaks match the integrated spectrum
   plot; the printed `omega_list` line pastes into a config that runs.

## 7. Checkpoint

Stop after §2–§3 and checks 1–4 and 6; report before the scripts.

## 8. Decisions (user, 2026-11-24)

- Entries as `N:value` pairs.
- The harmonic index lives in the run log; `radiation_field.dat` keeps its format.
- Peaks: local maxima above 10% of the window maximum.
