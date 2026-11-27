# Runs for Sec. IV.D (LG content) and IV.E (ponderomotive lens)

Decided 2026-11-27: the same four cases serve both sections; which go to the paper and which to the supplement is
chosen from the results.

- **LG (0, 2) and LG (2, 2)** (both l = 2: one ring against three), **a0 = 0.1 and 0.25**, circular polarization with
  the helicity of Figs. 2 and 3 (lab σ = −1, m = −2; selection rule l = −2, −5 for N = 1, 2), harmonics N = 1, 2.
- **10 µm, γ = 13.520 (route A), N_c = 30**, cutoff 6, `w0 = 10 λ_L`, polar-lattice sheet (radius from the LG table:
  2.37 w0 and 3.22 w0). Lens strength a0²N_c = 0.3 and 1.9, the range where the lens shows but the LG decomposition
  still works (the earlier a0²N_c = 2 runs: fidelity without d_N 0.92 / 0.56, with it 0.9999 / 0.9986). N_c = 8 would
  give 0.08 and 0.5. No 800 nm twin planned (optional, supplement).
- **Estimates** (Eqs. of manuscript III.A; κ = 0.541, 0.481): dressing span 0.27% / 0.24% (a0 = 0.1) and 1.7% / 1.5%
  (a0 = 0.25) against a line width of 1.25% (N = 1), 0.88% (N = 2); fringe number ≤ 1.3.
  ω_N,ref lies below the line by ~0.25% (a0 = 0.1, a fifth of the width) and ~1.5% (a0 = 0.25, more than the width).
- **Screen** L of `runs_sec4_plan.md`: Z/Z_* = 2.81. A weak-field Fresnel estimate (source u^N) puts 99.99% of the
  flux within 2.8–3.1 θ_* for (0, 2) and 3.7 / 4.8 θ_* (N = 1 / 2) for (2, 2). a0 = 0.1: lines and maps to 8 θ_*;
  a0 = 0.25: lines to 16 θ_* (the lens widens the emission), map radius set from them.

## Stages

1. `run_all.sh` (default ids): 8 line spectra (`*_line`, dense spectrum on a radial line, as `../fig_spectrum`) and
   the two a0 = 0.1 maps (`*_a0_0p1_map`, N = 1, 2 at ω_N,ref in one run). ~250 CPU h, ~2 h on 128 CPUs.
2. From the a0 = 0.25 lines: the peak of each screen-integrated spectrum (as `../fig_maps/peak_frequencies.py`) and the
   flux extent → `PEAKS` and `MAP['0p25']['theta_max']` in `make_configs.py`, regenerate, run the four
   `*_a0_0p25_N{1,2}_map` configs (one harmonic per run, `omega_list`).

Map screen: circular, 512 equal-area rings × 65 azimuths (64 distinct; `plot_observables` takes ∂_φ by FFT).
Lattice: 1 λ, and 0.5 λ for the a0 = 0.25 N = 2 lines (grating replica beyond 1.5 × the screen radius).

`make_configs.py` writes `configs/<id>.cfg` from `config/reference/spectrum_line_lattice_2d.cfg` (the IV.B base;
a line config differs from the matching `../fig_spectrum` one only in a0, p and the window). Don't edit the configs by
hand. `require_clean_build true`: commit before building. Outputs: `~/output/paper_2026/prelim/fig_lg/<id>/<timestamp>/`.

**Note:** the earlier test `../fig_maps/10um_a0_0p1_LG02` (LG (0, 2), a0 = 0.1, N_c = 30, same screen) has the
**opposite** laser helicity (`laser_zeta_2_im -1`): it is the anti-aligned case, not a repeat of
`10um_LG02_a0_0p1_map`.
