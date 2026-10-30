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

## Results

Not run yet.
