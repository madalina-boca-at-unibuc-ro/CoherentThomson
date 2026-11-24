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
