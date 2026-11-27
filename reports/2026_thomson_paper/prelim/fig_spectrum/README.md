# Preliminary runs: Sec. IV.B spectrum figure

Dense spectra along one radial screen line, to draft the Sec. IV.B figure (local (θ, ω) maps and the
screen-integrated spectrum) before the wording. Decided 2026-11-26.

- LG (0, 2), circular polarization, head-on, polar-lattice flat sheet; `N_c` (`laser_wing_sigma`) per case, cutoff 6
  (cutoff 4 leaks ~4e-4 of N = 1 into N = 2 in the `direct` formulas).
- 10 µm (γ = 13.520, D = 27.003) and 800 nm (γ = 3.884, D = 7.637; λ_L/D² = 13.72 nm at both), `w0 = 10 λ_L` at both: the 800 nm runs are the
  scaled twins.
- a0 = 0.1 (N_c = 30) and a0 = 1 (N_c = 8) for the main text, a0 = 0.5 (N_c = 30) for the supplement; one run per
  harmonic (N = 1, N = 2). a0 = 1 with N_c = 30 was run first (job 32515): a0²N_c = 30 spreads the emission beyond its
  16 θ_* line, so it was replaced by N_c = 8 (decided 2026-11-26). Probe runs (line to 80 θ_*) put all but ~1e-3 of
  the flux of the new cases within 16 θ_*; their lines go to 24 θ_* with the lattice replica beyond 36 θ_*.
- Screen: `Z/Z_*` = 2.81 (a0 = 0.1) and 2.85 (a0 = 1), the same at both wavelengths; the line reaches 8 θ_* and 16 θ_*.
  Earlier lines to 250 µrad clipped the flux: at a0 = 1, 10% (N = 1) and 17% (N = 2) of it was in the last 10% of the
  radius. The a0 = 0.1 N = 2 window was dominated by the 2.5 λ lattice's grating replica at ~280 µrad.
- Lattice spacing from the aliasing rule `θ_rep/θ_* = π w0 / ((ω/ω_1)_max d)`: 1 λ, and 0.5 λ for a0 = 1, N = 2
  (replica at 19.6 θ_*, beyond the 16 θ_* line).
- Frequencies in `ω/ω_1` (dressed with a0²/2): steps 0.0005 (a0 = 0.1) and 0.001 (a0 = 1, about 17 points per
  fringe).

`make_configs.py` writes `configs/<id>.cfg` from `config/reference/spectrum_line_lattice_2d.cfg`. Don't edit the
configs by hand: change the script and regenerate. Outputs go to
`~/output/paper_2026/prelim/fig_spectrum/<id>/<timestamp>/`.
`require_clean_build true`: commit before building.

Post-processing: `python3 py_scripts/radiation/plot_spectrum_map.py <run_folder>` (add `--linear` and
`--harmonic N` for the islands).

Estimated cost (from `~/output/20261124_175505`, a laptop run: 0.50 µs per electron, screen point, cycle and
frequency; 360 cycles at cutoff 6):

| runs | electrons | points × ω | CPU h | on 96 CPUs |
|---|---|---|---|---|
| a0 = 0.1, N = 1 or 2 (×4) | ~1800 | 128 × 401 | 2 | ~1 min |
| a0 = 1, N = 1 (×2) | ~1800 | 192 × 521 | 9 | ~6 min |
| a0 = 1, N = 2 (×2) | ~7200 | 192 × 1041 | 72 | ~45 min |

About 2 h in all on 96 CPUs. Running on the cluster: `run_all.sh` (header) builds once on the node, runs all 8,
then runs `plot_spectrum_map.py` on each. Use 96 CPUs, not `--exclusive`: above ~130 threads the a0 = 1, N = 2 field
buffers (~0.06 GB per thread) exceed the 8 GB cap.

Second batch (a0 = 1 N_c = 8, a0 = 0.5 N_c = 30; `run_all.sh` default ids): about 300 CPU h in all, measured cost on
the cluster ~(6.4 + 0.5 N_omega) µs per electron, screen point and cycle; ~2.5 h on 128 CPUs (N = 2 buffers ~0.04 GB
per thread).

**800 nm momentum changed 2026-11-27:** the twins first ran at γ = 3.952 (`average_pz 3.823 mc`, D = 7.775,
undressed λ_L/D² = 13.24 nm, not the 10 µm 13.71 nm). Now γ = 3.884 (`3.753 mc`): both twins give 13.72 nm in the
zero-intensity limit. The 800 nm runs were redone (third batch); the earlier ones stay under their older timestamps.
