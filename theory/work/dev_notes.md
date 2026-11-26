# Development notes

New entries, newest last. The record up to 2026-11-23 is frozen in `dev_notes_2026-07-28_to_2026-11-23.md` (paths
there predate the reorganization; `theory/archive/reorganization_plan.md` maps them).

## 2026-11-24: collision timing (spec 14)

- LG runs collided off the focus: the pulse peak met the beam at z = β/(1+β)·cutoff·σ·λ_L (49.9 λ_L for cutoff 4,
  74.9 for 6) while the focus is at the origin. With w0 = 10 λ_L that is 2.5–5.4% below the peak intensity, and N = 1
  at ω₁ (the red edge of the line at a0 = 1) changed 1.8× between cutoff 4 and 6. Found while chasing a fake N = 2
  signal; record: `~/Dropbox/work/AI-outputs/Claude/2026-11-24/screen_distance_and_harmonics_findings.md`.
- Spec 14 implemented to its checkpoint: `laser_delay` honoured, each electron starts at the pulse front, b(q) and
  the beam diagnostics use t = 0 positions. At the focus N = 1 no longer depends on the cutoff (ratio 1.000037).
  `theory/implemented/14_collision_timing.md`.
- Side result (n2_residue tests): the near-axis N = 2 of `simplified` at a0 = 1 was an O(dτ²) residue of the LR/BR
  cancellation (J_z = 0); use `long_distance_direct*` and cutoff 6 for weak harmonics. Record:
  `~/Dropbox/work/AI-outputs/Claude/2026-11-24/n2_residue/`.

## 2026-11-24: spectrum along a screen line (spec 15)

- `dense_frequency_spectrum` now runs on a line of screen points (rectangular `Nx` or `Ny` 1, circular `N_phi` 1),
  with a memory estimate and an 8 GB cap; `plot_spectrum_map.py` maps dP_z/dω, |E₊|, |E₋| over (θ, ω/ω₁) and gives
  the line-integrated spectrum. At a0 = 0.1 the integrated lines peak within their width of N ω₁ (offsets +0.0033,
  +0.0058). Next: the a0 = 1 scan the spec suggests, to see how far ω₁ is from the real spectrum.
  `theory/implemented/15_spectrum_map.md`.

## 2026-11-24: a0 scan and screen fields at arbitrary frequencies (spec 16)

- a0 scan on the line reference (0.1–1): harmonic N fills N × [(1 + a0²/2)/(1 + g a0²/2), 1 + a0²/2] ω₁ with g = 0.541
  for the code's LG (0, 2) (peak 0.736 a0), so ω₁ is below the line for a0 ≳ 0.3; at a0 = 1 the line breaks into
  islands in θ and fringes in ω (linear maps, `plot_spectrum_map.py --harmonic 1 --linear`).
- Spec 16: `omega_list` (`N:value`), N logged per frequency and read by the scripts, peaks and a ready `omega_list`
  line from `plot_spectrum_map.py`. `theory/implemented/16_omega_list.md`.

## 2026-11-24: N = 2 screen at a0 = 1 (end of session)

- `config/work/screen_N2_a0_1p0.cfg` (lattice 0.625 λ, 80 w0 screen at 160000 w0, 64 × 32, the five N = 2 peaks
  2.511–2.893 ω₁ of `~/output/20261124_184819`), run `~/output/20261124_191533` (65 min): Λ_zz/P_z = theory to 6
  digits at all five, pure E₋, J = −6. Main ring ~110–130 µrad and slowly decaying outer rings: 4–11% of the flux beyond
  0.9 R_max (30% at 2.893). Next: containment (R_max ~150 w0, N_R 128, below the lattice order at ~1.1 mrad; one
  frequency per run keeps the fast path), or dressing ω_N with the laser's real peak a² (a spec).
- Distance test (a0 = 1 line spectrum at 640000 w0 = 64 m, line to 160 w0 so 0–250 µrad as at 16 m; run
  `~/output/20261124_202914`): the (θ, ω) maps are the same as at 16 m (N = 1 indistinguishable, N = 2 marginally
  crisper), |E| ∝ 1/D and dP_z/dω ∝ 1/D² exactly, integrated peaks within 0.001 ω₁ and 3% in height. The island
  structure is the far-field pattern of the source, not a Fresnel effect: 16 m is already far field for it.

## 2026-11-26: agent instructions in one file

- `AGENTS.md` is now the single copy of the agent instructions; `CLAUDE.md` only imports it (`@AGENTS.md`), so edit
  `AGENTS.md`. Per-run findings (a0 = 1 N = 2 containment, the n2_residue record path) moved out of it to this log.
  `theory/README.md` + this file stand in for the Dropbox-wide `AI_NOTES.md`.

## 2026-11-26: trajectory_NT at a0 = 1 (fig_spectrum runs)

- Reduced copies of `reports/2026_thomson_paper/prelim/fig_spectrum/configs/10um_a0_1p0_N{1,2}.cfg` (lattice 2 λ,
  24 line points; `config/work/fig_spectrum_nt/`, runs `~/output_tests/fig_spectrum_nt/`) at NT 100, 200, 400: NT 100
  differs from NT 400 by ≤ 8e-6 (relative L2 of the dP_z/dω and |E₋| maps and the integrated spectrum), NT 200 by
  ≤ 4e-7; ratio ~20-25, i.e. close to 4th order (RK4), not the trapezoid's 2nd. NT 100 is converged at a0 = 1 for
  N ≤ 2 in `long_distance_direct_approx` (closes the a0 = 1, N ≤ 2 part of run-plan check C3). Script:
  `~/Dropbox/work/AI-outputs/Claude/2026-11-26/nt_convergence_a0_1.py`.

## 2026-11-26: fig_spectrum line runs (10 um and 800 nm twins, a0 = 0.1 and 1, N_c = 30, LG (0, 2))

- Cluster job 32515 (hpc-node13, 96 CPUs, commit 45a0e4f), 8 runs, ~1.6 h; copied to
  `~/output/paper_2026/prelim/fig_spectrum/`. Comparison script: `~/Dropbox/work/AI-outputs/Claude/2026-11-26/fig_spectrum_scaling.py`.
- a0 = 0.1: the scaled twins agree. Normalized integrated photon spectra differ by 9e-5 (N = 1) and 2e-4 (N = 2);
  radial flux profiles identical to 3 digits; edge flux 2e-5. Integrated peak +0.0035 (N = 1, FWHM 0.012) and +0.0055
  (N = 2) from N omega_1. Absolute photon numbers (same 1736 electrons) scale as (D_10um/D_800nm)^(2N) to 0.3-0.5%:
  12.02 and 144.8, against 12.06 and 145.5 (consistent with on-axis single-electron density ~ D^2 per harmonic order
  times a coherent solid angle ~ theta_*^2 ~ D^-4; to be checked).
- a0 = 1 (a0^2 N_c = 30): the emission fills the whole 16 theta_* line (flux per tenth of the radius roughly flat,
  rising at the edge for N = 2), so the integrated spectra are incomplete; normalized shapes differ between the twins by
  5% (N = 1) and 12% (N = 2). Estimate: the dressing phase across the sheet (~2 pi N sqrt(pi/2) N_c kappa a0^2/2 ~ 64 rad
  for N = 1) deflects light out to ~60 theta_*; at N_c = 8 it would be ~4x less.
- Correction: in those a0 = 1 runs the lattice replica sat at ~19.6 theta_*, just beyond the line, so part of the
  edge flux was probably its skirt. Probe runs (10 um, line to 80 theta_*, replica just beyond it; job 32528,
  `config/work/fig_spectrum_probe/`, `~/output_tests/fig_spectrum_probe/`, script `probe_extent.py` in the same
  AI-outputs folder): below 60 theta_*, a0 = 1 N_c = 8 and a0 = 0.5 N_c = 30 keep all but ~1e-3 of the flux within
  16 theta_* (N = 2 the widest: 0.17 and 0.36 in 8-16 theta_*). The probes' own edge pile-up (75-80 theta_*) is the
  replica: keep it beyond ~1.5x the line. Decided: a0 = 1 at N_c = 8 for the main text, a0 = 0.5 at N_c = 30 for the
  supplement; lines to 24 theta_*.
