# Coherent and incoherent parts without sampling every electron

The project has three steps:
1. **Flat (zero-thickness) disk:** the ideal, fully coherent case.
2. **Realistic Gaussian bunch length:** shows the loss of coherence.
3. **Gaussian bunch with ponderomotive microbunching** (a technique from the literature): how much coherence is
   gained back.

A lattice beam fixes step 1 (no artificial speckle) and serves as the baseline. For steps 2 and 3 the electron
number problem stays unless the coherent and incoherent parts are computed separately, as below.

## 1. Why random sampling with N_sim electrons is wrong for steps 2 and 3

For N independent electrons with phase-space density ρ (normalized to 1), and single-electron field e(x; X₀) for
an electron starting at X₀ = (r₀, z₀, p₀):

  ⟨|Σⱼ eⱼ|²⟩ = N ⟨|e|²⟩ + N(N−1) |⟨e⟩|²,  ⟨·⟩ = ∫ ρ(X₀) (·) dX₀.

- First term: incoherent (∝ N). Second: coherent (∝ N²).
- A real bunch has N ~ 10⁸–10¹⁰. A random sample of N_sim (10⁴–10⁵) gets the coherent/incoherent ratio wrong by
  N_real/N_sim, so a simulated "coherence loss" would describe the simulation, not the beam.
- The current zero-thickness disk is fully coherent longitudinally at λ_s ≈ 14 nm. Its speckle (0.05–1%) is a
  sampling artefact, far above the physical shot noise of a real disk. Reports should say which physical
  regime the disk represents.

## 2. Weighted macroparticles: right for the coherent part only

Each simulated electron ("node") can stand for a phase-space cell of w = N_real/N_sim electrons that move almost the
same way. The cell's coherent field is w·e(node), as long as the density is smooth over the cell: quadrature of ⟨e⟩.

But squaring the weighted sum treats each cell as one point charge of w electrons. Its incoherent intensity
then comes out as w²|e|² instead of the correct w|e|², overestimating the incoherent part by N_real/N_sim, the same
error as random sampling. Inside a real cell the electrons have random positions on the λ_s scale and add
incoherently.

So two separate sums over the same nodes:
- **coherent field:** F_coh = Σ_nodes w·e, a genuine field (the bunch's mean field);
- **incoherent intensity:** Σ_nodes w·|e|² = N⟨|e|²⟩, and the same for each quadratic observable, computed per
  node from that electron's own field;
- **total** = |F_coh|² + incoherent (N(N−1) ≈ N²).

For any quadratic quantity Q (energy flux, `Σ_zz`, `Λ_zz`, …): Q_total = Q[F_coh] + N⟨Q[e]⟩. "Total F" of a single
shot is never computed: these are shot-averaged results. A single-shot speckle realization can be generated
afterwards from the same per-electron fields with random phases if fluctuations are needed.

Same idea as FEL codes: regular "quiet start" loading suppresses artificial shot noise, and the physical noise is then
added in a controlled way (Fawley, PRST-AB 5, 070701 (2002)). Here only averages are needed, so the incoherent term is
computed directly.

## 3. The three steps

1. **Flat disk** (**implemented 2026-10-30**: `beam_distribution polar_lattice` with `beam_lattice_spacing`, rings
   at `r_j = j d` with `round(2πj)` equally spaced nodes, `(J + ½) d = R`; `IoUtils::polar_lattice_geometry`).
   Validated on the rest setup (`ponderomotive_lens.md`, Results):
   - the random beam is unchanged byte for byte;
   - spacing 2.44 λ (1321 electrons, 39 s on 256 CPUs) against 1.22 λ: the screen field changes by
     6e-3 (N = 1) and 9e-3 (N = 2) in amplitude, and the scale follows the electron count to 0.2%;
   - against the random 65536-electron run (1601 s): residual power 1.1e-3 (N = 1) and 6.8e-3 (N = 2), that run's
     speckle;
   - power at other `l`: 1e-7 (N = 1) and 1e-4 (N = 2), against 7e-4 and 6e-3 for the random run.

   The plan: F_coh from regularly placed nodes, equal-weight hexagonal or polar rings (the C++ sums electrons
   without weights). Spacing: grating orders at `θ = λ_s/d` must fall outside the screen, `d < λ_s/(2θ_max)` (≈ 11 λ
   for the (0,0) run), and the amplitude and dressing phase must be sampled, `d ≲ w0/5` (≈ 5 λ). That gives
   `N ≈ π R_b²/d² ≈ 300` instead of 28,000–65,000, with no speckle. Convergence: stability under halving d. Polar
   rings with equal azimuthal spacing keep the sampling exactly axisymmetric, so `l` purity becomes a clean test.
2. **Gaussian bunch length:**
   - *Coherent:* an electron starting at z₀ follows the same trajectory shifted in time, so `e(z₀) = e(0) e^{iqz₀}`
     (valid while the bunch length ≪ z_R ≈ 2000 λ; real bunches are 3–30 λ_L). The z-integral is then the analytic
     form factor `ρ̃(q)`: coherent part = step 1 × `|ρ̃(q)|²`. This is the `b(q)` diagnostic the code already prints
     (`Simulation::compute_form_factor`).
   - **Open:** the correct q. The dev notes leave it between `q = k_s − k_L` and the dressed
     `q = k_s − (K_s·P / K_L·P) k_L`. **Cheap test:** two electrons at different z₀ (same r₀), read the phase difference of
     their fields at each screen point. Seconds of compute.
   - A few extra z₀ samples test the time-shift factorization (focusing changes across the bunch).
   - *Incoherent:* N⟨|e|²⟩. In the far field an electron at transverse offset r₀ gives the same angular pattern
     apart from a phase, so `|e|²` depends only on |r₀| and the azimuthal average is a rotation. So one electron per
     radius: ~10–20 single-electron runs, then Python.
3. **Microbunching:** the bunched density enters through ρ. Its form factor has peaks `|bₙ|²` at `q = n k_b`, so
   coherent = `|b(q)|²` × the single-electron field. If the bunching comes from an energy modulation, momenta are
   correlated with z₀: sum over a few momentum slices, each with its own trajectory and form factor. The incoherent
   term is unchanged. All deterministic, a few hundred to a few thousand trajectories in total.

Assumptions: electrons independent apart from their one-particle distribution (no space charge, no correlations
beyond ρ); shot-averaged results; bunch length ≪ z_R.

**Check of the decomposition itself (no code change):** the formula holds for any N. From a single-electron
computation of ⟨|e|²⟩, it must reproduce the speckle fraction of the existing random runs (0.05% for
`20261029_084449`, 0.9% for `20261027_181823`) and the empirical `(kθR)²/(2πN)` law
(numerical_parameters_and_cost.md §4).

## 4. LG decomposition of a partially coherent beam

- **Coherent part, unchanged:** projection is linear, so `a_pl^coh = ⟨u_pl|F_coh⟩ = N⟨a_pl(e)⟩`. Largest overlap,
  least squares, series waist and waist-position scan all apply as now.
- **Incoherent part, a mode matrix:** project each node's field separately (vector a_j) and accumulate
  `ρ_inc = Σⱼ wⱼ aⱼ aⱼ†`. Diagonal = mode powers (additive); off-diagonal = coherences between modes. The beam:
  `ρ_total = a_coh a_coh† + ρ_inc`, the LG projection of the cross-spectral density
  `W(x₁, x₂) = F_coh(x₁)F_coh*(x₂) + N⟨e(x₁)e*(x₂)⟩`.
- **What carries over:**
  - power in a mode: `⟨u|W|u⟩/Tr W` (the diagonal), the analogue of the largest-overlap diagnostic;
  - OAM per photon: `Σ_l l ρ_ll / Tr ρ` + spin, additive, so it splits into coherent and incoherent contributions;
  - degree of coherence: `Tr(ρ²)/(Tr ρ)²`, or the coherent power fraction: the number step 2 is about;
  - least squares applies only to a field (the coherent part). For the total, diagonalize `ρ_total`: its
    eigenvectors are the beam's natural coherent modes (Wolf's coherent-mode representation), its eigenvalues their
    powers.
- **Expected physics:** a single electron's pattern is ~1/γ wide (~74 mrad) against a ~60 μrad screen. On the screen
  each electron's field is almost a plane wave whose tilt depends on where it started, so the incoherent part is a nearly
  uniform background spread over the ~(kθ_max R_b)² modes the screen resolves. Consequences:
  - the coherent fraction depends on the angular acceptance: the coherent part is fully collected once contained,
    the incoherent power keeps growing ∝ θ_max². Quote coherence for a stated aperture (physically right: an
    experiment collects through one);
  - for an axisymmetric beam `ρ_inc` is block-diagonal in l: an electron at azimuth φ₀ is a rotated copy of one at
    φ₀ = 0, so its coefficients pick up phases `e^{−i(l+s)φ₀}`, and the average over φ₀ removes the l ≠ l'
    coherences. One electron per radius suffices. Within one l, different p stay coherent.
- **Implementation (Python only):** `lg_decomposition.py` already computes the projections. Accumulate `w a a†` over
  the nodes, report the diagonal, the l-resolved OAM, the coherent fraction for the screen's acceptance and optionally
  the eigenmodes of `ρ_total`. Truncate the basis at the modes the screen resolves and report the incoherent
  power outside it, as `Σ|a|²` does now.

## 5. Code changes needed

- C++: a beam option with lattice/quadrature nodes next to `generate_cylinder_beam` (equal weights, or add
  per-electron weights to the coherent sum).
- Python: the incoherent average from single-electron runs, the form factor and combining with N_real, and the
  partially coherent LG decomposition.
