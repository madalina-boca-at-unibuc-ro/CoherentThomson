# Preliminary runs: Sec. IV.C field and observable maps

Full-screen maps at the spectral peak, to decide how Sec. IV.C (angular momentum, selection rule) presents them.
Decided 2026-11-26: a0 = 1, N_c = 8, LG (0, 2), circular polarization, N = 1 and 2, 10 µm and its 800 nm scaled twin.

- **Physics identical to Fig. 2** (b), (c): each config is the matching `../fig_spectrum/configs/<id>.cfg` with only
  the screen, the frequency and the lattice spacing changed (`make_configs.py`; `diff` the two to see it).
- **Frequency:** one `omega_list` entry per run, the largest maximum of that run's screen-integrated photon spectrum
  (`peak_frequencies.py`): N = 1: 1.2683 (10 µm), 1.2682 (800 nm); N = 2: 2.5401, 2.5385 (ω/ω₁). Each twin uses its
  own peak; they differ by 1e-4 and 6e-4 relative, ≤ 2% of the line width.
- **Screen:** circular, Z/Z_* = 2.85 (as the spectrum lines), R = 16 θ_* (45.6 w0). At the peak, 4e-6 of the flux is
  beyond 16 θ_*, 5e-4 beyond 12 θ_*. N_R = 1024 equal-area rings (first ring off axis at 0.5 θ_*), N_phi = 65
  (64 distinct azimuths; `plot_observables` takes ∂_φ by FFT, J = −3 and −6 need far fewer).
- **Lattice:** 1 λ (N = 1, 1736 electrons) and 0.5 λ (N = 2, 7088 electrons), replica at 1.55 × 16 θ_*.
- **Cost (estimate):** ~2–3 µs per electron, screen point and cycle (laptop smoke run: 1.9 µs), 96 cycles:
  ~9 CPU h per N = 1 run and ~38 per N = 2 run, ~95 CPU h in all: ~45 min on 128 CPUs, plus post-processing.
  Field buffers ~19 MB per thread.

`make_configs.py` writes `configs/<id>.cfg`; don't edit those by hand. `require_clean_build true` (inherited):
commit before building. Outputs: `~/output/paper_2026/prelim/fig_maps/<id>/<timestamp>/`. `run_all.sh` (header) builds
on the node, runs all four, then post-processes each run: `scaled_parameters.py --write-log`, `plot_observables.py`,
`plot_field_carrier_removed.py total|plus|minus`, `lg_decomposition.py`.

`10um_a0_0p1_LG02.cfg` is an earlier hand-made test config (a0 = 0.1, N = 1, 2 at ω_N,ref), kept for reference.
