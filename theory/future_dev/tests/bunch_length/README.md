# Bunch-length test: does the lens survive a Gaussian bunch, with only the amplitude reduced?

Prepared 2026-10-30. Step 2 of `../../coherent_incoherent_decomposition.md`, as a numerical check: run B of the
progress report (LG (0,0), a0 = 0.2, γ ≈ 13.5 head-on, R_max 25 w0, 65536 random electrons;
`config/expansion_a0_scan_server/config_expansion_p0l0_a0_02_R25.cfg`, run `~/output/20261029_143739`) with a
Gaussian bunch length `beam_sigma_z` (zero-height cylinder, so z is a pure Gaussian).

## Expectation

The bunch is far shorter than z_R, so an electron starting at z₀ follows the disk electron's trajectory shifted in
time: `e(z₀) = e(0) e^{iqz₀}`. Then for the same transverse sample

  F_σ = ρ̃(q) F_0 + R,   ρ̃(q) = exp(−q²σ²/2),   ⟨|R|²⟩ = N⟨|e|²⟩ (1 − ρ̃²),

with `q_N = N k_1 + k_L` (head-on, on axis; the screen spans ~1e-4 rad, so q is constant over it to ~1e-8). The
coherent field keeps its shape (waist, z0, `l`, radial phase), its amplitude falls by ρ̃, and the residual R is
incoherent noise whose power relative to the coherent part grows as (1 − ρ̃²)/ρ̃². The dressed q
(`q = k_s − (K_s·P/K_L·P) k_L`) differs from this one by ~6e-5 relative, which this test can't resolve.

## Runs

| Config | σ_z [λ] | σ_z [nm] | q₁σ | q₂σ | ρ̃ (N = 1) | ρ̃ (N = 2) |
|---|---|---|---|---|---|---|
| `p0l0_a0_02_R25_sigma_ref.cfg` | 1e-9 | ~0 | 4.5e-6 | 9e-6 | 1 | 1 |
| `p0l0_a0_02_R25_sigma_qs025.cfg` | 5.5581e-5 | 0.556 | 0.25 | 0.4997 | 0.9692 | 0.8827 |
| `p0l0_a0_02_R25_sigma_qs05.cfg` | 1.1116e-4 | 1.112 | 0.5 | 0.9993 | 0.8825 | 0.6070 |
| `p0l0_a0_02_R25_sigma_qs1.cfg` | 2.2232e-4 | 2.223 | 1.0 | 1.999 | 0.6065 | 0.1357 |

(λ_L = 10 µm; k₁ from run B's `frequencies_list`.) Only `output_folder` (`~/output_tests/bunch_length`) and
`beam_sigma_z` differ from run B's config.

**Why a near-zero reference instead of run B itself:** `generate_electron` draws the Gaussian z between the
transverse draws only when σ_z > 0. So σ_z = 0 (run B) has a different transverse sample from every σ_z > 0 run,
while all σ_z > 0 runs share one: the same standard normals, scaled by σ_z (checked with a standalone copy of the
draw order: transverse positions bit-identical, z/σ equal to 4e-16). The reference is the flat disk on that sample,
so the comparison is point by point, not up to speckle. Run B vs the reference is a free check that two
realizations agree within speckle. Check `electron_beam_scatter.dat` x, y are identical across the four runs.

## Running

Inside a compute-node job (`theory/srun.txt`; ~2.5 h for the four runs, so `interactive` (4 h) or `cpu-wide`):

```
bash theory/future_dev/tests/bunch_length/run_batch.sh
```

## Analysis

Per harmonic N, over the screen (all six F components, or E_−):
- **amplitude law:** `c = ⟨F_ref, F_σ⟩ / ‖F_ref‖²`; expect `|c| = ρ̃(q_N)`, arg c ≈ 0;
- **shape:** `1 − |⟨F_ref, F_σ⟩|² / (‖F_ref‖² ‖F_σ‖²)`; expect ≈ (1 − ρ̃²)/ρ̃² × the reference's incoherent
  fraction η (speckle, 1e-3–1e-2). The factor is 0.06/0.28/1.7 (N = 1) and 0.28/1.7/53 (N = 2) for the three runs,
  so the mismatch should grow by that ratio from run to run, and N = 2 at q₂σ = 2 is mostly noise;
- **lens:** `lg_decomposition.py` on each run (bounded `--z0` list, as in the server README): same dominant mode,
  w0, z0 and share as the reference; `radial_phase` maps unchanged.

## Results (2026-10-30)

Runs on `hpc-node02`, commit `088b7ad` (clean), Release, ~37 min each; config → run in
`~/output_tests/bunch_length/logs/index.txt`: ref `20261030_143857`, qs025 `_151536`, qs05 `_155222`, qs1 `_162911`.
Transverse positions are identical across the four runs.

**Amplitude law** (`compare_bunch_length.py` here; F total, six components; η = resid / (1 − ρ̃²)):

| run | N | qσ | ρ̃ | \|c\| | \|c\|/ρ̃ | arg c | mismatch | η |
|---|---|---|---|---|---|---|---|---|
| qs025 | 1 | 0.25 | 0.9692 | 0.9686 | 0.9994 | 2.6e-3 | 1.9e-4 | 2.9e-3 |
| qs025 | 2 | 0.50 | 0.8827 | 0.8798 | 0.9968 | 8.0e-3 | 8.5e-3 | 3.0e-2 |
| qs05 | 1 | 0.50 | 0.8825 | 0.8803 | 0.9975 | 5.2e-3 | 8.7e-4 | 3.1e-3 |
| qs05 | 2 | 1.00 | 0.6070 | 0.6002 | 0.9889 | 1.4e-2 | 5.5e-2 | 3.3e-2 |
| qs1 | 1 | 1.00 | 0.6065 | 0.5994 | 0.9882 | 9.4e-3 | 6.3e-3 | 3.6e-3 |
| qs1 | 2 | 2.00 | 0.1357 | 0.1270 | 0.9355 | −4.0e-2 | 0.67 | 3.4e-2 |

- The coherent amplitude falls by `exp(−q²σ²/2)` to ≤ 1.2% at qσ = 1. The C++ `|b|` at N = 1 (0.9693, 0.8828,
  0.6074) also matches ρ̃.
- η is constant across the runs (N = 1 ≈ 3e-3, N = 2 ≈ 3e-2), so the residual is the incoherent remainder
  ∝ (1 − ρ̃²). The mismatch follows (1 − ρ̃²)/ρ̃² · η.
- The remaining deviation depends on qσ alone (the same at N = 1 and N = 2 for equal qσ): `1 − |c|/ρ̃ ≈ 0.012 (qσ)²`
  and `arg c ≈ 0.010–0.015 qσ`. That is the expansion of `c = Σ_j w_j e^{iqσ g_j}` (w_j = electron j's share of
  F_ref, g_j its standard-normal z draw) with weighted moments Σ w g ≈ 0.01, Σ w g² ≈ 1.024. It's a property of the
  one z sample, well within the bound `|c − ρ̃| ≲ sqrt(η (1 − ρ̃²))`, not physics.

**Lens** (`lg_decomposition.py`, bounded z0 list; E_−; series waist from c₁/c₀; |a| relative to ref):

| run | N | dominant | share | Σ\|a\|² | `l` purity | series w0 [λ] | series z0 [λ] | \|a\|/\|a_ref\| | beyond 0.9 R |
|---|---|---|---|---|---|---|---|---|---|
| ref | 1 | `LG_0^0` | 0.919 | 0.995 | 0.9975 | 8.547 | 2.807e5 | 1 | 0.10% |
| qs025 | 1 | `LG_0^0` | 0.920 | 0.995 | 0.9973 | 8.555 | 2.811e5 | 0.969 | 0.11% |
| qs05 | 1 | `LG_0^0` | 0.921 | 0.995 | 0.9966 | 8.562 | 2.810e5 | 0.882 | 0.14% |
| qs1 | 1 | `LG_0^0` | 0.917 | 0.987 | 0.991 | 8.571 | 2.777e5 | 0.601 | 0.36% |
| ref | 2 | `LG_0^{−1}` | 0.808 | 0.950 | 0.976 | 4.045 | 2.048e5 | 1 | 0.75% |
| qs025 | 2 | `LG_0^{−1}` | 0.802 | 0.944 | 0.967 | 4.015 | 2.062e5 | 0.880 | 1.0% |
| qs05 | 2 | `LG_0^{−1}` | 0.764 | 0.894 | 0.919 | 4.027 | 2.055e5 | 0.600 | 2.4% |
| qs1 | 2 | `LG_0^{−1}` | 0.263 | 0.306 | 0.331 | 4.035 | 2.061e5 | 0.126 | 22% |

- **The lens survives the bunch length unchanged.** Same dominant mode, same basis waist (z0 = 3e5 / 2e5 λ), and the
  series waist moves by < 1% in w0 and ~1% in z0 even at q₂σ = 2, where the coherent part holds a third of the power.
  `|a|` scales by ρ̃.
- Share, purity and Σ|a|² fall only by the added incoherent power: share_σ ≈ share_ref · ρ̃² / (ρ̃² + (1 − ρ̃²) η)
  gives 0.765 (qs05, N = 2; measured 0.764) and ≈ 0.29 (qs1, N = 2; 0.263). The "beyond 0.9 R" rise is the
  speckle filling the screen, not a wider beam.
- `radial_phase` maps: at q₂σ = 2 the same spiral is visible on axis at ~1/7 the reference amplitude, with speckle of
  comparable size over the whole screen.
- The reference vs run B (`~/output/20261029_143739`, another transverse sample): N = 1 share 0.919 vs 0.915, series
  z0 2.807e5 vs 2.815e5; N = 2 share 0.808 vs 0.800, z0 2.048e5 vs 2.029e5. They agree within speckle.

**Conclusion:** a Gaussian bunch length only multiplies the coherent field by `exp(−q_N²σ²/2)` and adds the
incoherent remainder; the ponderomotive lens (mode, waist, waist position) is unaffected for σ ≪ z_R.
