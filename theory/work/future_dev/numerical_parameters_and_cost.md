# Numerical parameters, cost and minimal settings

Context: the (2,2) runs became expensive (an `R_max = 40 w0` case ran ~7 h on 256 cores). Goal: estimate the
minimal value of each parameter, then test around it.

## Cost model

From `run_log.txt` of `20261029_084449` (machine `ad1`, GCC 11.5, `long_distance_simplified_approx`, 2 harmonics,
`trajectory_NT 100`): **3.2 μs CPU per electron × screen point × laser cycle**, i.e.

  wall time ≈ 3.2 μs × N_e × N_screen × N_cycles × (trajectory_NT/100) / N_cores.

This reproduces that run's 2670 s (65536 electrons, 4096 points, 100 cycles, 32 threads). Every parameter enters
linearly.

## 1. `trajectory_NT` (verified)

`dτ = T / (trajectory_NT × (p⁰ − n·p)/mc)` (`simulation.cpp`): `trajectory_NT` points per laser cycle **as seen by
the electron**. In a plane wave `k·p` is conserved, so the laser phase is linear in τ even with dressing. Each laser cycle
seen by the electron becomes one emitted N = 1 period, so harmonic N gets `trajectory_NT/N` points per period,
independent of γ.

Scan (`tests/nt_scan.py`): 16 electrons on a 1 w0 disk, Gaussian pulse σ = 25 cycles (cutoff 4), screen scaled to each γ's
emitted angular width, relative L2 difference of the total F (6 components) against NT = 1600 with the same formula. At
`p_z = 13.48 mc`:

| NT | N = 1 simplified | N = 1 direct | N = 2 simplified | N = 2 direct |
|---|---|---|---|---|
| 25 | 1.5×10⁻⁴ | 1.5×10⁻⁴ | 9.4×10⁻³ | 3.0×10⁻⁴ |
| 50 | 5×10⁻⁶ | 5×10⁻⁶ | 2.5×10⁻³ | 1.1×10⁻⁵ |
| 100 | 3×10⁻⁷ | 3×10⁻⁷ | 6.3×10⁻⁴ | 6×10⁻⁷ |
| 200 | 2×10⁻⁷ | 4×10⁻⁸ | 1.6×10⁻⁴ | 3×10⁻⁷ |
| 400 | 2×10⁻⁷ | 5×10⁻⁹ | 4×10⁻⁵ | 3×10⁻⁷ |

("simplified" = `long_distance_simplified_approx`, "direct" = `long_distance_direct_approx`.)

- **γ:** no real dependence for `p_z` = 2, 13.48, 50 mc. At rest the simplified N = 2 error is ~3× smaller. At
  50 mc the direct N = 2 has a floor of ~7×10⁻⁶ (precision, `k_s Z ≈ 10¹¹`).
- **Formula:** `direct` converges spectrally (~30× per doubling), `simplified` as O(dτ²) (4× per doubling). In
  `simplified` the exact endpoint term cancels a large trapezoid sum, so the trapezoid's O(dτ²) endpoint error
  remains. In `direct` the integrand vanishes at both ends. The two agree with each other to the expected accuracy at
  every NT and γ.
- **Recommendation:** `long_distance_direct_approx` with `trajectory_NT 50` (≤ 10⁻⁵ for N ≤ 2). Rule of thumb
  NT ≳ 25 N. Checked only for N ≤ 2 and a0 = 0.1; recheck for higher harmonics or larger a0 (more harmonic content in
  the trajectory). With `simplified`, keep NT = 100 (6×10⁻⁴ at N = 2).

## 2. Screen: define it by angle (estimates, to test)

- **Distance:** the far-field angular scaling (the image at 2Z on a 2R screen equals the image at Z on an R
  screen) needs `Z ≳ 10 z_R` of the emitted beam to hold to ~1%. For `w0 = 25 λ` and `λ_s = λ/725.5`:
  `z_R = π w0²/λ_s ≈ 1.4×10⁶ λ`, so the current `Z = 4×10⁶ λ` is only 2.8 z_R. The long-distance formulas don't
  cost more at larger Z: use `Z ≈ 1.5×10⁷ λ`.
- **Half-angle θ_max:** about 3× the beam's angular radius `θ_N ≈ λ_s/(π w0,eff)`, times `√(2p+|l|+1)` for an LG mode.
  For (0,0) N = 1: `θ_N ≈ 17.5 μrad` (25 λ source), ~24 μrad with the dressing. The current 62.5 μrad is ≈ 2.6 θ_N,
  with 0.3% of the power beyond 0.9 R_max: adequate.
- **N_R:** set by the smooth envelope structure, **not** by the wavefront fringes. The exported field carries
  `e^{ikR}`, whose fringe count grows ∝ Z at fixed angle, but the intensities and observables use only |·|² and ∂φ,
  and the LG basis carries the same curvature, so the projection integrand stays smooth. Information about the waist
  position survives at any Z (the residual phase `k θ² z0/2` doesn't depend on Z). Criterion: `Δ(θ²)/θ_N² ≲ 0.15`
  (the decomposition flags > 0.25), i.e. `N_R ≈ (θ_max/θ_N)²/0.15 × (p+1)`: 64–128.
  - Note: the old rule "R_max × 2 needs N_R × 4" holds at fixed Z, where the fringes must stay resolved for the
    Re/Im plots. In angle terms it is automatic.
- **N_phi:** set by the largest `|l|` (Nyquist) and the azimuthal-spectrum check: 32–64. With the FFT ∂φ (done) the
  OAM no longer needs a fine φ grid.
- The equal-area grid is coarse near the axis (first ring at `θ_max/√(N_R−1)`). That matters for future
  vortex-core work (spectral_line_and_vortex_core.md), not now.
- **Symmetry reduction (later):** for an on-axis LG laser with circular polarization and an axisymmetric beam,
  `E_±(θ, φ) = E_±(θ, 0) e^{ilφ}`, so one radial line (N_phi = 1) would be 64× cheaper. Only once the pure-`l`
  property is established (now true with a Gaussian envelope, pulse_envelope_and_lg_normalization.md), and it removes the purity check.

## 3. Beam radius

Electrons outside the laser emit nothing, so `R_b` only needs to cover the region that radiates. The source of
harmonic N goes as `a(r)^N`, so higher harmonics need less. The LG table radius (1.86 w0 for (0,0), 3.22 w0 for
(2,2)) is "just larger than the mode". To test once runs are cheap: N = 1 results at 1.86 vs 2.5 w0. The truncation
edge sits at a source amplitude of ~3% (`e^{−R_b²/w0²}` for N = 1), which may leave a weak diffraction tail.

## 4. Electron number: speckle law (two runs agree)

Incoherent (speckle) power fraction on the screen, from random positions:

  **f_inc ≈ (k θ_max R_b)² / (2π N_e)**,  `k = 2π/λ_s`.

The numerator is the number of independent speckle modes the screen resolves (space-bandwidth product).

| Run | N_e | R_b | θ_max | predicted | measured (power outside the expected `l`, N = 1) |
|---|---|---|---|---|---|
| `20261027_181823` | 8192 | 75 λ (3 w0) | 62.5 μrad | 0.89% | 0.93% |
| `20261029_084449` | 65536 | 46.5 λ (1.86 w0) | 62.5 μrad | 0.043% | 0.05% |

N_e needed for a target f_inc grows as `θ_max² R_b²`. For (0,0) at 10⁻³: ~28,000 electrons. For (2,2) (larger R_b and
θ_max): a few 10⁵, hence the 7 h. Random sampling is only a stopgap: see coherent_incoherent_decomposition.md for
why its speckle level isn't physical anyway, and the deterministic replacement.

## 5. Proposed production settings (no code change)

| | (0,0) | (2,2) |
|---|---|---|
| Envelope | `flat_duration 0`, `wing_sigma 25`, `cutoff 4` | same |
| Formula, NT | `long_distance_direct_approx`, 50 | same |
| Z | 1.5×10⁷ λ | same |
| Screen | θ_max ≈ 3 θ_N, N_R 64–128, N_phi 64 | larger θ_max for the N = 2 ring (`l = −5`) |
| N_e (random, f_inc ~10⁻³) | ~28,000: ~40 min on 32 cores | ~3×10⁵: ~8 h on 32 cores, ~1 h on 256 |
| N_e (lattice, once implemented) | ~300–1000: minutes | ~1000–3000: minutes |
