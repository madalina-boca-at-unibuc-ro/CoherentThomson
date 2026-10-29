# Pulse envelope and LG normalization

Two properties of the laser model, found 2026-10-29 while investigating the asymmetric N = 2 energy flux of run
`20261029_084449` (LG `p = l = 0`, `w0 = 25 λ`, `a0 = 0.1`, circular `ζ = (1, i)`, flat top 50 cycles + Gaussian
wings `σ = 5` cycles, `γ ≈ 13.5` head-on, forward screen at `4×10⁶ λ`).

## 1. The flat-top envelope leaks N = 1 radiation into N = 2

### Symptom

At N = 2, `E_−` should be pure `l = −1` (`J = N(m+σ) = −2`, spin −1). The (0,0) run has 90.75% at `l = −1` and
**8.9% at `l = 0`**. The `l = 0` part interferes with the `l = −1` ring, which explains the off-axis crescent in
`|F^{01}|` and `dP_z/dω` and the 9% low `Λ_zz/P_z` at N = 2 (0.913 × theory; the `l` content predicts 0.907).

### Not electron sampling

- The `l = 0` share does not change with electron density: 9.0% in the 8192-electron run `20261025_002302`, 8.9% in
  the 65536-electron run, at ~21× the surface density. Over the same change, the speckle (all `l`) fell from ~1% to 0.05%.
- **One electron exactly on axis** (`tests/one_electron_on_axis.cfg`), with no sampling and exact axial symmetry, gives N = 2
  `E_−`: 84.9% `l = −1`, **15.1% `l = 0`**.

### Not the radiation code

On-axis spectrum of that electron (`tests/one_electron_on_axis_spectrum.cfg`). On axis the true N = 2 (`l = −1`)
vanishes, so whatever is at 2ω₁ there is the forbidden part:

| Variant | Peak in 1.9–2.1 ω₁ (÷ N = 1 peak) |
|---|---|
| Current run settings (`long_distance_simplified_approx`, NT 100) | 4.95×10⁻⁹ |
| `long_distance_simplified`, `long_distance_direct_approx`, `simplified`, `direct` | 4.95×10⁻⁹, identical |
| `trajectory_NT 400` | 4.95×10⁻⁹ |
| `laser_type plane_wave` | 1.46×10⁻⁹ |
| `laser_wing_sigma 10` cycles | 2.33×10⁻¹⁰ |

The true N = 2 is ~4×10⁻⁹ of N = 1 in total (screen-integrated `P_z`: 9.4×10⁻⁹ vs 2.32), so this is enough to give
the 9–15%.

### Mechanism: the N = 1 line's far wing

- The flat top is joined to the Gaussian wings with a continuous first derivative only, so the N = 1 line's spectral
  wings fall off as a power law. Away from 2ω₁ the floor is a0-independent: plane wave, `[1.4, 1.6] ω₁`:
  2.7×10⁻⁸ (a0 0.02), 2.6×10⁻⁸ (0.05), 2.4×10⁻⁸ (0.1), 1.5×10⁻⁸ (0.2). That matches the estimate from the jumps in
  `g''` (≈ 2/σ² at each joint).
- Exactly at Δ = 1 harmonic spacing, the linear part cancels for this symmetric pulse with a whole number of cycles
  (50). The envelope's Fourier transform at the carrier is 8×10⁻¹² numerically. What remains at 2ω₁ is the dressing
  phase in the wings, `ψ ∝ a0² ∫ g²`, which breaks the cancellation. Plane wave, `[1.9, 2.1] ω₁`: 2.1×10⁻¹⁰ (0.02),
  2.7×10⁻¹⁰ (0.05), 1.5×10⁻⁹ (0.1), 1.9×10⁻⁸ (0.2). The excess grows roughly as a0⁴ relative to N = 1, while the true
  N = 2 grows as a0², so the contamination share grows as a0².
- This radiation carries the N = 1 angular momentum (J = −1, so `l = 0` in `E_−`).

### Why the axial symmetry doesn't forbid it

For circular polarization, rotating the system about z equals shifting the laser's carrier-envelope phase (CEP). The setup
is invariant under rotation **combined with** a CEP shift, not under rotation alone. Pure `l = J − s` per harmonic
additionally requires a CEP shift to act like a time shift, which an envelope with non-smooth joins violates.
Test: multiplying both ζ by `i` (CEP + π/2) rotates the N = 2 bump of the one-electron field from 90.4° to 180.4°,
with the same `|E|²` asymmetry (`l = 1` / `l = 0` harmonic of the strongest ring: 0.274). N = 1 stays symmetric.

### Fix: pure Gaussian envelope (config only)

`laser_flat_duration 0` makes the two wings meet at the peak, giving one exact Gaussian `exp(−φ²/σ²)`. Its
spectral wing at one harmonic spacing is `∝ exp(−σ²/4)` (σ in rad), negligible for σ ≳ 2 cycles.

| Pulse | Simulated length | N = 1 FWHM | Floor at 2ω₁ | N = 2 `E_−` at `l = 0` |
|---|---|---|---|---|
| flat 50 + wings 5 (current) | 100 cycles | 0.0149 ω₁ | 5×10⁻⁹ | 15% (1 e⁻), 9% (beam) |
| Gaussian σ = 25, cutoff 5σ | 250 cycles | 0.0153 ω₁ | ~10⁻²¹ | 4×10⁻¹⁵ |
| Gaussian σ = 25, cutoff 4σ | 200 cycles | — | 10⁻¹⁷ | — |
| Gaussian σ = 15, cutoff 5σ | 150 cycles | 0.0225 ω₁ (coarse grid) | ~10⁻²¹ | — |

Proposal: `laser_flat_duration 0`, `laser_wing_sigma 25 cycles_adim`, `laser_wing_sigma_cutoff 4`. It keeps the current
line width at 2× the trajectory length. With the `direct` formula at NT = 50 (other note) the total cost stays the
same as now.

Differences to expect from the flat top: the dressing now changes throughout the pulse. At a0 = 0.1 the on-axis line is
slightly broader and shifted (centroid 0.9988 vs 0.9951), and the dressing lens is averaged over the pulse, so waists
and waist positions change. At a0 = 0.02 this is negligible.

### Side observation (not understood)

In the LG beam, an electron exactly on axis leaves the pulse with a net transverse momentum of 2.6×10⁻⁵ of its peak
(flat top; 6.4×10⁻⁶ with 10-cycle wings). In the plane wave it is 7.9×10⁻¹², equal to the envelope's Fourier
transform at the carrier. The kick doesn't cause the floor (the plane wave has the floor without the kick), but where it
comes from in the focused beam (first-order paraxial `E_z`/`B_z`?) is open.

## 2. LG laser: on-axis field = √2 × plane wave at the same `a0`

- `laser_field_heatmap_z0.dat` at the waist: LG (0,0) on-axis `|E_⊥|` = 4.5528×10⁻⁴, plane wave 3.2216×10⁻⁴ (same
  units; ratio **1.4132**). This comes from the custom `Npn` normalization. So with circular polarization the
  on-axis dressing is `⟨a²⟩ = a0²` for the LG beam, while the frequency formula uses `a0²/2`, which is correct for the
  plane wave.
- One electron on axis, N = 1 line centre (fine grid, flat top):
  - LG, a0 = 0.1: **0.99515** (75% of peak intensity at the computed 1.000);
  - LG, a0 = 0.02: 0.99977;
  - plane wave, a0 = 0.1: 1.00008.

  The shift scales as a0², and its size, a0²/2, is what the √2 predicts.
- **For the beam it matters much less:** the 1024-electron beam's line at the on-axis screen point peaks at **1.00210**
  (95% at the computed frequency). Most of the coherent line comes from off-axis electrons, which feel a lower a(r)
  and are dressed less. The nominal frequency is resonant where `a(r)² = a0²`, i.e. `r ≈ 0.59 w0`.
- **To decide:** keep the config `a0` as the plane-wave-equivalent value and state the on-axis peak in reports and in
  `run_log.txt`, or redefine `a0` as the LG on-axis peak. Redefining changes the meaning of every earlier LG run. For
  `p, l ≠ 0` the mode's peak isn't on axis, so the definition needs care. Until then: "`a0 = 0.1`" for an LG (0,0)
  laser means a peak `a ≈ 0.141` on axis.
- CLAUDE.md's rule "`<a²> = a0²/2`, don't change it to `a0²`" was established for the plane-wave normalization and
  stays right there.
