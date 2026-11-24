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
  cancellation (J_z = 0); use `long_distance_direct*` and cutoff 6 for weak harmonics.

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
