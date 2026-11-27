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

## 2026-11-26: second fig_spectrum batch (a0 = 1, N_c = 8; a0 = 0.5, N_c = 30)

- Job 32529 (hpc-node04, 128 CPUs, a622484), 8 runs, 2.8 h; copied to `~/output/paper_2026/prelim/fig_spectrum/`.
  All contained (edge flux 1e-6 to 1e-3), aliasing margins 1.57-1.64, no `scaled_parameters.py` warnings.
- Band (Eq. harmonic_band): 93-98% of the energy inside; 5-95% spans match the edges (a0 = 1, N = 1: 1.188-1.490
  against 1.181-1.500).
- 10 um vs 800 nm: normalized integrated spectra differ by 0.5% (a0 = 1, N = 1), 1.5% (N = 2), 1.7% and 4.7% (a0 = 0.5,
  N = 1, 2), against 1e-4 at a0 = 0.1; photon ratios (D_10/D_800)^(2N) to 0.5-3%. The lines reach 24 theta_*, where
  the D^2 theta^2 term of omega_N(theta) is ~3e-3 at 800 nm (`angular_detuning`), a fraction of a fringe: the likely
  cause (not yet separated from epsilon = 2/(D k_L w0) = 4e-3 at 800 nm).
- Fringes: integrated spectra show 3 (a0 = 1, N = 1; N_f 2.7) and 2 (a0 = 0.5, N = 1; N_f 2.5) peaks above 10%, but
  one smooth peak at N = 2 (N_f 5.4, 5.1): the N = 2 fringes tilt with theta and most N = 2 flux is at 4-8 theta_*,
  so the angular integral washes them out. Local spectra at theta = theta_* show 3 (N = 1) and 4 (N = 2) peaks.

## 2026-11-26: teorie.ro as a stand-in server

- The cluster's NetBird SSO `sso.acc.unibuc.ro` refuses TCP 443 (ping works; `vpn.acc.unibuc.ro` answers): server side.
- Key login to teorie.ro set up from the Dell G15 (SELinux label fix on the CephFS `~/.ssh`); details in
  `theory/cluster.md` §3. Its checkout was clean at dde5ad6.
- Cluster reachable from the Dell G15 without NetBird: `ssh hpc` (ProxyJump teorie.ro, user `madalinaboca`);
  the cluster checkout was clean at a622484 (one commit behind).

## 2026-11-26: twin agreement of the fig_spectrum runs, quantified

- Normalized integrated spectra, 800 nm vs 10 um (max |diff| of peak; centroid shift; width change):
  a0 = 0.1 N_c = 30: 9e-5 / 2.5e-4 (N = 1 / 2), shifts ~1e-6; a0 = 1 N_c = 8: 7.3e-3, -4.1e-4, -0.2% (N = 1),
  1.7e-2, -1.1e-3, -1.0% (N = 2); a0 = 0.5 N_c = 30: 3.0e-2 (N = 1), 5.9e-2, -2.7% width (N = 2).
  800 nm is always red-shifted and narrower. (Supersedes the 0.5-5% of the earlier entry, which was a looser metric.)
- The mismatch grows with the angular cut of the integral (a0 = 1, N = 2: 2.7e-3 to 3 theta_*, 7.3e-3 to 6, 1.7e-2 to
  12): consistent with the (D theta)^2 term of omega_N(theta), ~1/D^2 at fixed theta/theta_*, i.e. 4x larger at
  800 nm. D theta_* ~ epsilon, so the "D^2 theta^2 vs epsilon" question of the earlier entry may be one parameter
  (not checked). Scaling in theta/theta_* holds up to O((D theta)^2), not exactly.
- Fig. 2 caption: note that the twins differ by <= 0.7% / 1.7% in (b), (c) (user: caption only, no figure change).
- The 12 fig_spectrum runs were copied to the Dell G15 (`~/output/paper_2026/prelim/fig_spectrum/`, 1.8 GB).
- Scripts (run from anywhere, read `~/output/paper_2026/prelim/fig_spectrum/`):
  `reports/2026_thomson_paper/figures/checks/spectrum_{twin_diff,twin_cut,fringes,photon_ratio}.py`.

## 2026-11-26: fringes in Fig. 2 (b), (c): local vs integrated (a0 = 1, N_c = 8, LG (0, 2), 10 um)

- Fewer fringes than at N_c = 30 is expected: Eq. fringe_number gives N_f = 2.7 (N = 1), 5.4 (N = 2) at N_c = 8,
  against ~10 and ~20 at N_c = 30 (N_f proportional to N_c; same a0^2 N_c ~ 8 as a0 = 0.5, N_c = 30).
- Local spectra (prominence > 0.05) have the fringes at every angle. N = 1: 2-3 peaks, spacing ~0.13-0.14 in
  omega/omega_1 (theta = 3: 1.224, 1.426; 6: 1.316, 1.448; 8: 1.220, 1.366, 1.494). N = 2: 3-4 peaks, spacing
  ~0.13-0.18 (theta = 4: 2.418, 2.594, 2.832; 6: 2.394, 2.576, 2.728, 2.860; 8: 2.540, 2.688, 2.838; 10: 2.648,
  2.770, 2.890).
- Flux in theta (10%-50%-90% of the screen-integrated dP_z): N = 1 at 0.8-2.5-6.7 theta_*, N = 2 at 2.7-5.5-8.7.
- Integrated: N = 1 keeps the fringes partly (peaks 1.268, 1.392, second prominence 0.17, plus a shoulder ~1.5,
  i.e. ~3 bumps = N_f); N = 2 is one smooth peak at 2.540.
- Reading: the fringes move with theta. Over the flux-carrying range a fringe moves by ~0.1-0.15 at N = 2 (about one
  spacing), so the angular integral fills the minima; at N = 1 the flux sits at smaller theta and the shift is about
  half a spacing, so they survive. Why they move (not checked): different theta collect mainly from different sheet
  radii with a different peak dressing, plus the angular term of omega_N(theta).
- For IV.B: local spectra keep N_f fringes (Eq. fringe_number); the integrated spectrum keeps them only where the
  fringe shift across the flux-carrying angles is below the spacing. Possible supplement panel: local spectra at
  theta = 2, 6, 10 theta_* for N = 2 (not made).

## 2026-11-26: absolute photon-number ratio of the twins (10 um / 800 nm)

- Same N_e within each twin pair (polar lattice, weight 1; 1736 at a0 = 0.1, 3526 / 14314 at a0 = 0.5 N = 1 / 2,
  7088 / 28056 at a0 = 1 N = 1 / 2: the lattice spacing follows the aliasing condition). D from omega_1/omega_L = D^2/(1 + a0^2/2):
  D = 27.0 (10 um), 7.775 (800 nm), D_10/D_800 = 3.473.
- Normalization: `dN` of `fig_spectrum.load` is flux per unit *screen area*; photons per solid angle need x Z^2,
  and the twins' screens differ (Z = 43524 w0 vs 3608 w0). Total = integral over omega and solid angle (omega_1
  cancels); peak = max of dN/d(omega/omega_1).
- N_ph(10 um)/N_ph(800 nm) against (D_800/D_10)^(2N) = 0.08290 (N = 1), 0.006872 (N = 2) [total; deviation of
  total; deviation of peak]:
  a0 = 0.1 N_c = 30: 0.08318, +0.35%, +0.35% (N = 1); 0.006907, +0.52%, +0.53% (N = 2).
  a0 = 1 N_c = 8: 0.08331, +0.49%, +0.08% (N = 1); 0.006965, +1.4%, +0.08% (N = 2).
  a0 = 0.5 N_c = 30: 0.08400, +1.3%, -0.01% (N = 1); 0.007101, +3.3%, -0.54% (N = 2).
- So 800 nm gives 12.0x (N = 1) and 144x (N = 2) more photons at fixed N_e. The peaks follow D^(-2N) to <= 0.5%
  (0.1% at a0 = 1); the totals carry the 800 nm redshift/narrowing (entry above) and deviate by up to 3%.
- Only at equal N_e: with w0 = 10 lambda the 800 nm sheet is 12.5x smaller in radius, so at equal areal density it
  holds 12.5^2 fewer electrons (coherent yield ~ N_e^2). The D^(-2N) explanation in the paper's AI_NOTES is still
  "to check".


## 2026-11-26: manuscript III.B against spec 11; `linewidth` fix

- III.B used $k_{N,\rm ref}$ for $Z_*,\theta_*$ and $\omega/\omega_{N,\rm ref}$ for spectra; the code (spec 11) uses $k_{1,\rm ref}$
  for every N and $\omega/\omega_{1,\rm ref}$. Green corrections proposed in the TeX (see the paper's AI_NOTES), including
  the angular-detuning equation, which is the module's `angular_detuning` ($(\theta/\theta_*)^2/[(k_1/k_L)(k_Lw_0)^2]$).
- `scaled_parameters.py`: `linewidth` = $0.375/(\sqrt N N_c)$ (Eq. `pulse_linewidth`) instead of $1/(NN_c)$; recorded in
  `theory/implemented/11_scaled_parameters.md`. Checked: 0.0468 (N = 1) and 0.0331 (N = 2) at $N_c = 8$.

## 2026-11-26: Sec. IV.C field-map runs prepared (not run)

- `reports/2026_thomson_paper/prelim/fig_maps/` (README there): a0 = 1, N_c = 8, LG (0, 2), N = 1, 2, 10 um and 800 nm;
  configs = the fig_spectrum ones with a circular screen (Z/Z_* = 2.85, R = 16 theta_*, 1024 x 65), one `omega_list`
  frequency at the spectrum's largest maximum (`peak_frequencies.py`: 1.2683 / 1.2682, 2.5401 / 2.5385), lattice 1 and
  0.5 lambda (replica margin 1.55). Flux beyond 16 theta_* at the peak: 4e-6.
- Smoke test (shrunk copy, laptop): config parses, frequency mode `omega_list` (N logged), 1.9 us per electron, point
  and cycle; `scaled_parameters --write-log`, `plot_observables`, `plot_field_carrier_removed`, `lg_decomposition` all
  run on it. Estimate ~95 CPU h in all (~45 min on 128 CPUs).

## 2026-11-27: Sec. IV.C field-map runs done (job 32536)

- Job 32536 (hpc-node04, 128 CPUs, 522d106 clean), 54 min with post-processing; runs copied to the Dell G15
  (`~/output/paper_2026/prelim/fig_maps/`, 780 MB). Solver 2.97-3.14 us per electron, point and cycle (wall 280 s per
  N = 1 run, 1145 s per N = 2 run).
- Selection rule, from `plot_observables` (Lambda_zz/P_z against (N(m + s) - s)/omega_N, m = -2, s = -1):
  10 um: exact to 7 digits at N = 1 (L = -2) and N = 2 (L = -5); 800 nm: 2.8e-6 (N = 1), 2.3e-6 (N = 2).
  Sigma_zz/Lambda_zz = 1/2 (N = 1), 1/5 (N = 2). Flux/density ratios = c to 1e-7 (10 um), 5e-6 to 1e-5 (800 nm).
- Helicity: E_- holds the power; E_+ fraction 2.6e-9 / 6.7e-9 (10 um, N = 1 / 2), 3.9e-7 / 9.8e-7 (800 nm; ~100x larger,
  as epsilon^2 and the angular spread grow at smaller D). Azimuthal purity of E_-: l = -2 (N = 1), -5 (N = 2), other l
  below 1e-10 (N = 1), 1e-15 (N = 2); E_+ (tiny) at l = J - 1 = -4, -7 as required.
- Containment: 0.3-0.6% of the E_- power beyond 0.9 R (R = 16 theta_*).
- Photon ratio from the screen-integrated P_z (10 um / 800 nm): 0.08296 (N = 1), 0.006882 (N = 2), against
  (D_800/D_10)^(2N) = 0.08290, 0.006872: +0.07%, +0.14% (cf. the spectral peaks, entry of 2026-11-26).
- Twin maps (`figures/checks/maps_twin_compare.py`: E_- with each run's carrier removed, global phase fitted, same
  scaled grid): N = 1: 1 - overlap 8e-7 / 1.3e-5 / 3.4e-5 within 4 / 8 / 16 theta_*, relative L2 0.13% / 0.5% / 0.8%;
  N = 2 (frequencies 6e-4 apart): 1.4e-4 / 2.1e-4 / 3.3e-4, L2 1.7% / 2.1% / 2.6%. Grows with theta, as the spectra.
- Maps: E_- is a vortex on axis (charge -2, -5) inside many concentric rings with a strong spiral (radial) phase: the
  ponderomotive lens at a0 = 1. N = 1: bright first ring at ~0.8 theta_*, rings fading outwards; N = 2: ~10 rings out to
  ~6 theta_*.
- `lg_decomposition.py` with its defaults (p <= 3, automatic w0 scan, z0 = 0) fails here: E_- captured power 57% (N = 1,
  best-eps waist at the scan edge) and 0.7% (N = 2, overlap maximum at the scan edge), series waist NaN.
- Why (2026-11-27, Legion, scratch script on the azimuthal projection f_l(theta) of E_-, origin carrier removed): the
  radial phase is **not** a curvature (no z0 fixes it; quadratic fit leaves 1.9 / 3.6 rad rms). It is tent-shaped and
  piecewise **linear in theta** (N = 1: ~ +2 rad per theta_* out to the peak region, then falling; N = 2: ~ +5 rad per
  theta_* up to 6 theta_*, then -4), with steps at the near-nodes, and the rings are **equally spaced in theta**
  (N = 2: ~0.55 theta_*), not in theta^2. So the lens is strongly aberrated (not parabolic), consistent with the
  paper's IV.E d_N(r) picture; the evenly spaced rings suggest two interfering contributions (hypothesis, untested:
  e.g. rising and falling a^2). No single-w0 LG set with small p represents it. The 10 um and 800 nm profiles coincide in theta/theta_*. Runs copied to the Legion as well.

## 2026-11-27: 800 nm twins moved to gamma = 3.884 (all 800 nm runs redone)

- The twins ran at gamma = 3.952 (`average_pz 3.823 mc`, D = 7.775): undressed lambda_L/D^2 = 13.24 nm, not the 10 um
  13.71 nm (3.823 had been chosen to hit 13.5 nm with the a0 = 0.2 dressing). Decided: gamma = 3.884 (`3.753 mc`,
  D = 7.637), 13.72 nm at both in the zero-intensity limit. Everything above that compares twins (spectra, peak
  frequencies, maps, photon ratios (D_800/D_10)^(2N)) used the old gamma; the 800 nm results are to be replaced.
- Commit d510cf2 (fig_spectrum configs regenerated; only `average_pz` and the scaled screen distance change), brought
  to the cluster as a git bundle (not yet pushed to GitHub). Job 32562 (hpc-node14, 128 CPUs): the six 800 nm spectrum
  runs (a0 = 0.1 N = 1, 2; a0 = 1 N_c = 8; a0 = 0.5 N_c = 30); the superseded a0 = 1 N_c = 30 pair is not redone. Then:
  peak frequencies, fig_maps configs, the two 800 nm map runs.
- Done (job 32562, spectra, d510cf2; then fig_maps 800 nm, e64ad70; both clean, hpc-node14, 128 CPUs; runs copied to
  the Legion). New 800 nm peaks = old ones to 4 digits (1.2682, 2.5385).
- Photon ratio 10 um / 800 nm against (D_800/D_10)^(2N) = 0.07999, 0.006398 (`figures/checks/spectrum_photon_ratio.py`):
  a0 = 0.1: +0.36%, +0.54% (line-integrated); a0 = 1 N_c = 8: +0.09% at the peak (N = 1, 2), +0.5%, +1.4% over the
  band; a0 = 0.5 N_c = 30: -0.01%, -0.6% at the peak, +1.4%, +3.5% over the band. Maps (screen-integrated P_z per unit
  omega/omega_1 at the peak): +0.073%, +0.085%. All of order epsilon (4.2e-3 at 800 nm) or of the angular detuning.
- Normalized twin spectra (`spectrum_twin_diff.py`, max |difference|): a0 = 0.1: 1e-4, 2.6e-4; a0 = 1 N_c = 8: 0.76%,
  1.8% (were 0.7%, 1.7%); a0 = 0.5: 3.1%, 6.2%.
- Twin maps (`maps_twin_compare.py`, relative L2 within 4 / 8 / 16 theta_*): N = 1: 0.14% / 0.54% / 0.86%; N = 2:
  1.6% / 2.0% / 2.6% (unchanged from the old gamma). 800 nm AM per photon: spin -0.999994, orbital -2.000006 (N = 1);
  -0.999988, -5.000012 (N = 2); E_+ power fraction 4.2e-7, 1.1e-6.
- Fig. 2 regenerated from the new runs (labels now N_gamma); manuscript IV.B numbers updated.

## 2026-11-27: paper runs copied to the Dell G15 laptop

- `~/output/paper_2026/prelim/` (fig_spectrum: 16 ids, fig_maps: 4 ids + `10um_a0_0p1_LG02`) now on the Dell G15 for the
  next days' work. Checked against `hpc:~/output/paper_2026/prelim/`: same files and sizes (sizes only, no checksums),
  every run `STATUS: finished` with `radiation_field.dat`. Two run logs (800 nm `a0_1p0_Nc8_N1/N2`, 20261127 runs)
  lacked the later "Spectrum map" section: re-copied from the cluster. `fig_maps/10um_a0_0p1_LG02/` and
  `fig_maps_preview.png` exist only on the laptop copy, not on the cluster.

## 2026-11-27: Sec. IV.D/E runs prepared (not run)

- `reports/2026_thomson_paper/prelim/fig_lg/` (README there): LG (0, 2), (2, 2); a0 = 0.1, 0.25; 10 um, N_c = 30,
  cutoff 6, Z/Z_* = 2.81. Stage 1: 8 line spectra + 2 a0 = 0.1 maps (~250 CPU h). Smoke runs of shrunk LG (2, 2) copies
  on the laptop: configs parse, table radius 3.22 w0, frequencies as intended.
- Weak-field Fresnel estimate (u^N): 99.99% of the flux within 2.8-3.1 theta_* (LG (0, 2)), 3.7 / 4.8 theta_*
  (LG (2, 2), N = 1 / 2).
- The old `fig_maps/10um_a0_0p1_LG02` test has `laser_zeta_2_im -1`, the opposite helicity of Figs. 2-3.
