# Future development: findings and plans not yet acted on

Collected 2026-10-29 while preparing the first stage report (runs `20261029_084449`, Gaussian laser, and the
unfinished (2,2) run). Nothing here is implemented unless marked **done**. Each note gives the measurements, how
they were obtained, and the proposed next step.

| Note | Topic | Status |
|---|---|---|
| [pulse_envelope_and_lg_normalization.md](pulse_envelope_and_lg_normalization.md) | The flat-top envelope leaks N = 1 radiation into N = 2 (the N = 2 crescent / 9% `l = 0`); the LG laser's on-axis field is √2 × the plane wave's at the same `a0` | Findings verified; fix = Gaussian envelope (config only); `a0` convention to decide |
| [numerical_parameters_and_cost.md](numerical_parameters_and_cost.md) | `trajectory_NT` convergence vs γ and formula; cost model; screen defined by angle; electron number and speckle law | `trajectory_NT` verified; the rest are estimates to test |
| [coherent_incoherent_decomposition.md](coherent_incoherent_decomposition.md) | Replacing random electron sampling by quadrature nodes, with coherent and incoherent parts computed separately; the three project steps (flat disk, Gaussian bunch, microbunching); LG decomposition of a partially coherent beam | Plan |
| [spectral_line_and_vortex_core.md](spectral_line_and_vortex_core.md) | Mode content across the emission line; purity near the vortex core (relevance for twisted-photon absorption experiments) | Deferred by choice |

**Done the same day:** `plot_observables.py`'s azimuthal derivative is now spectral (FFT), exact for
`|l| < (N_phi−1)/2`. The old centered difference read `L_z`/`Lambda_zz` low by `1 − sin(lΔφ)/(lΔφ)`: 0.17% at
`l = 1`, 4.1% at `l = 5` with `N_phi = 64` (details in `dev_notes.md`).

**Suggested order:**
1. Production settings that need no code change: Gaussian envelope, `long_distance_direct_approx` with
   `trajectory_NT 50`, `Z ≳ 10 z_R`, screen sized by angle.
2. Two cheap checks: the two-electron test that settles the form-factor `q`, and checking the coherent/incoherent
   decomposition against an existing random-sampled run.
3. Lattice/quadrature beam option, then the parameter scans (beam radius, screen half-angle, `N_R`, `N_phi`).
4. Incoherent part and partially coherent LG decomposition (steps 2 and 3 of the project).
5. Spectral line and vortex core.

Test inputs are in [tests/](tests/): `one_electron_on_axis.cfg` (the (0,0) run's config with one electron at
`r = 0`), `one_electron_on_axis_spectrum.cfg` (the same, single on-axis screen point, dense spectrum 0.5–2.5 ω₁)
and `nt_scan.py`. Each run takes seconds. Their output goes to `~/output_tests/`, so the plotting scripts' "latest run"
in `~/output` is not affected.
