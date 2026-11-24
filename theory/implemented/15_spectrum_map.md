# 15: Spectrum along a screen line (implemented)

Spec: `theory/specs/15_spectrum_map.md`.

- **Status:** done 2026-11-24. The user asked for the full implementation and a merge, so the checkpoint (§7) was a
  self-check, not a stop.
- **Commit:** see `git log` (branch `spec15-spectrum-map`, merged into `main`).

## Code

- `Detector_2D::is_line()`: rectangular with exactly one of `Nx`, `Ny` = 1; circular with `N_phi` = 1 and `N_R` > 1.
- `main.cpp`: a dense spectrum runs on one point or a line; any other grid keeps its first point (message as
  before, now naming the grid). Before the run it prints the memory estimate and throws above
  `Simulation::kDenseSpectrumMemoryCap` (8 GB of buffers).
- `Simulation::resolve_num_threads` (the thread count `run_simulation` uses, now shared) and
  `Simulation::estimate_field_memory` (per-thread `PackedFaraday` buffers, the reconstructed `Faraday` field,
  `radiation_field.dat` at 1300 bytes per row, measured 1295).
- `run_log.txt`, "Frequency spectrum" (dense runs): the screen (one point or a line, its θ range) and the memory
  estimate.
- `py_scripts/radiation/plot_spectrum_map.py`: maps of dP_z/dω, |E₊|, |E₋| over (θ, ω/ω₁) (log; |E₊| and |E₋| share one
  scale), the line-integrated spectrum (annulus weights, edges halfway between points), and a summary per harmonic
  window written to `run_log.txt` ("Spectrum map (plot_spectrum_map.py)"). Reuses `compute_point_observables`,
  `get_observation_directions`, `helicity_components`, `_get_E_total`, `_write_run_log_section`.
- `plot_point_spectrum{,_all,_observables}.py`: refuse runs with more than 4 screen points
  (`plot_point_spectrum.require_few_points`).
- `config/reference/spectrum_lattice_2d.cfg`: a commented line example. CLAUDE.md: the script, and the line in the
  Frequencies bullet.

## Deviations and judgment calls

- **`|E₊|` and `|E₋|` share one colour scale** (not in the spec): separate scales hid that one helicity is ~10⁶ weaker.
- **Spec §6.5 expected |E_x| to depend on φ.** That holds only when both helicities are comparable (the a0 = 1 N = 2
  residue). At a0 = 0.1 one helicity dominates at N = 1 and 2, so |E_x| ≈ |E_dominant|/√2 is φ-independent too
  (equal on both lines to 1e-3 and 2%). The reason to map E± stands for the general case.
- `resolve_num_threads` returns at least 1 (was `min(threads, electrons)`, 0 for an empty beam).

## Checks (configs and logs: `~/Dropbox/work/AI-outputs/Claude/2026-11-24/spec15_checks/`; runs in
`~/output_tests/spec15/`)

Base: `config/reference/spectrum_lattice_2d.cfg` (a0 = 0.1, LG, 284 electrons, direct LD approx) shrunk.

| # | Check | Result |
|---|---|---|
| 1 | single point, spec-14 binary vs new (`N_omega` 200) | `radiation_field.dat` byte-identical |
| 2 | 3-point line (x = 0, 7, 14 w0) vs single-point runs at 7 and 14 w0 | rows byte-identical |
| 3 | dense (ω/ω₁ 1 to 2, 101 points) vs a harmonic run on the same line | identical at ω₁, 5.3e-12 at 2ω₁ |
| 4 | dP_z/dω of the line vs `compute_point_observables` on the single-point run | identical |
| 5 | the same 3 points on the φ = 90° line (y axis) | N = 1: dP_z/dω, \|E₊\|, \|E₋\| to ≤1e-3 (lattice discreteness, 2.5 λ). N = 2: dP_z/dω 3.3e-2, dominant \|E₋\| 1.6e-2; the weak \|E₊\| (~1e-4 of \|E₋\|, lattice-dominated) differs by O(1) |
| 6 | memory guard: 3 points × 400000 ω | refused, 9.2 GB (24 × 345.6 MB + 921.6 MB) > 8 GB, exit 1 |
| 6 | 64 points × 2000 ω: estimate vs measured | buffers 0.98 GB estimated, peak RSS 1.11 GB (+13%); file 166.4 MB estimated, 166.2 MB actual; 3 min 51 s |
| 7 | a0 = 0.1, 64 points to 87.5 µrad, ω/ω₁ 0.5–2.5 | integrated peaks at 1.0033 (FWHM 0.014) and 2.0058 (FWHM 0.021): within the line width of N ω₁; 82% (N = 1) and 95% (N = 2) of the window's energy within ±1% of N ω₁ |

## After the merge: a0 scan and the summary windows (2026-11-24)

`config/work/spectrum_line_a0_{0p1,0p2,0p5,1p0}.cfg` (copies of `reference/spectrum_line_lattice_2d.cfg`, cutoff 6,
ω range up to 2(1 + a0²/2) + 0.3); runs `~/output/20261124_173749`, `_174310`, `_174840`, `_175505`.

- **The windows [N − ½, N + ½) were wrong for broad lines** (at a0 = 1 the N = 1 line reaches 1.5 ω₁, so "N = 2"
  peaked at 1.5). Now each harmonic's band is N × [(1 + a0²/2)/(1 + g a0²/2), 1 + a0²/2] ω₁, with g the laser's peak
  |a|²/a0² (`peak_a2_fraction`: 1 for a plane wave, the C++ LG amplitude's maximum otherwise, 0.541 for (0, 2)), and
  windows split halfway between neighbouring bands. The map and spectrum draw the red edge as a dotted line, and the
  summary gives the energy fraction in the band.
- **Results** (N = 1; N = 2 alike, scaled):

  | a0 | integrated peak (ω/ω₁) | band | energy within ±1% of ω₁ | energy in the band |
  |---|---|---|---|---|
  | 0.1 | 1.0035 | 1.0023–1.005 | 0.82 | 0.18 (band ≪ pulse bandwidth 0.014) |
  | 0.2 | 1.015 | 1.009–1.02 | 0.25 | 0.57 |
  | 0.5 | 1.0785 | 1.054–1.125 | 1e-4 | 0.91 |
  | 1.0 | 1.4665 | 1.180–1.5 | 7e-11 | 0.97 |

  ω₁ from `frequencies_list` is below the red edge for a0 ≳ 0.3: there is essentially nothing there at a0 = 1, the
  reason the earlier fixed-ω₁ maps looked like thin-ring tails.
