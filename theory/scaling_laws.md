# Scaling laws: dimensionless parameters for the flat-disk results

> **Status: DRAFT, NOT APPROVED FOR IMPLEMENTATION.**
> Do not implement anything from this note (the collapse-test runs of §5, the detector units of §6, or any other
> change to code, configs or docs because of it) until the user explicitly confirms. When you are given this file:
> read it, check it against the code and the notes it cites, and reply with your analysis only (open questions,
> errors, anything that disagrees with the code), then stop and wait. Remove this notice only when the user approves.

Written 2026-11-13. A note, not a spec: nothing here changes the solver. Its purpose is the publication plan
(flat 2D lattice → noise in flat and 3D random beams → microbunching), where the first part should be stated in
variables in which the results don't depend on the laser wavelength or the beam energy. Only the microbunching part
needs real numbers.

**Claim.** For a cold, flat electron disk hit head-on by an LG pulse, with `γ ≫ 1`, every dimensionless result
(field maps in scaled angles, LG content, `l` purity, fitted waist in `w0`, waist position in emitted Rayleigh
ranges, photon number per electron) depends only on the parameters of the first table below. `λ_L` drops out
exactly. `γ` drops out up to corrections of order `ε²`, with `ε = 2/(D k_L w0)` (about 1e-3 for our runs). The
existing runs already show it (§4).

## 1. Why

**Classical electrodynamics has no length scale.** Lorentz force plus Maxwell are invariant under `x → s x`,
`t → s t`, `F → F/s`. `a0 = eE/(mcω)` is unchanged, and so are all phases. With every length given in `λ_L` or `w0`
(as the configs already do), `λ_L` is not a parameter. The classical photon number, `∫ dW/(ħω)`, is also scale-free;
only photon energies carry `ħω`.

**A boost along the axis removes `γ`.** Go to the frame where the beam is initially at rest, with Doppler factor
`D = γ(1+β) = u⁺` (`≈ 2γ`). In that frame:

- `a0`, the laser phase `φ = k_L·x`, and hence the envelope `g(φ)` and the number of cycles, are invariant;
- transverse lengths are invariant: `w0`, the disk radius `R`, the lattice spacing `d`, the electron positions;
- a flat disk (all electrons at the same `z` and velocity) stays flat;
- the laser wavenumber is `k' = D k_L`, so the rest-frame focusing parameter is `D k_L w0`;
- the laser's longitudinal fields, `~ θ_d E⊥` with `θ_d = 2/(k_L w0)`, become `~ (θ_d/D) E'⊥`: `E_z` is invariant
  while `E⊥` grows by `D`;
- emission at rest-frame angle `θ'` from the screen normal reaches the lab at `tan(θ/2) = tan(θ'/2)/D`, and
  `k_N = D k'_N`.

So an electron's dynamics is fixed by invariant functions of `φ`: the local amplitude `a0 M(ρ)`, the mode phase
`l ϕ`, and `g(φ)`. Two things are not functions of `φ` alone:

1. **How the mode changes along the worldline** (Gouy phase, curvature, `w(z)`), through the distance the electron
   travels through the focus while the pulse passes. That is the parameter `ξ` (§2). It depends on `k_L w0`, the
   envelope and `β/(1+β)`, and `β/(1+β) → 1/2` for `γ ≫ 1`.
2. **The longitudinal fields and the rest-frame emission pattern across the emitted beam**, both of order
   `ε = 2/(D k_L w0)`. When `ε ≪ 1`, the single-electron (dipole-like) pattern is flat over the coherent beam, whose
   rest-frame divergence is `~ε`. When `ε → O(1)`, the polarization structure of the single-electron pattern enters
   the coherent beam: spin–orbit coupling in the emission becomes physics, not a correction.

The coherent source of harmonic `N` at the disk, `S_N(ρ)`, is therefore universal in `ρ = r/w0`. Its far field is a
Fourier transform at `k_N`. Under the boost `θ k_N = θ' k'_N`, so the far field is universal in `θ/θ_N`, with
`θ_N = 2/(k_N w0)`. A paraxial beam is fixed by its far field, so the field at a finite distance is universal in
`x/w0` and `Z/Z_N`, with `Z_N = k_N w0²/2`.

## 2. The parameters

### Physics

| Parameter | What it sets | Planned runs (10 µm, `w0 = 8.4 λ`) | Report runs A–D |
|---|---|---|---|
| `a0` | Nonlinearity: harmonic `N` scales as `a0^N` in amplitude, plus the dressing | 0.1 | 0.1, 0.2 |
| Polarization `ζ`, mode `(p, l)` | Spin and orbital AM of the drive | circular; LG₀⁰, LG₀² | circular; LG₀⁰, LG₂² |
| Envelope `g(φ)`, in cycles | Line width `~1/n_eff`. `n_eff = ∫g² dφ / 2π`, `σ√(π/2)` for `g = e^{−φ²/σ²}` | `σ = 25`: `n_eff = 31.3` | same |
| `k_L w0` | Only through `ξ` (and `ε`) | 52.8 | 157 |
| `R/w0` | Truncation of the mode by the beam | 2.37 (fixed beam, all modes) | LG table (1.86 for LG₀⁰) |
| `d/w0` | Lattice sampling; grating replicas at about `π w0/d` in units of `θ_N` | 0.30 | — (random) |

Derived combinations that carry the physics:

| Combination | Meaning | Value |
|---|---|---|
| `Λ_N ~ N a0² M² n_eff` (× O(1)) | **Ponderomotive lens.** The dent's emission delay is a phase, so it is Lorentz invariant. The exact form is the dressing factor `d_N(r)` of `future_dev/ponderomotive_lens.md`, written in `φ`, `A(r, φ)` and `ω/ω_1`: all invariant. That note's "Not covered: a moving beam" is answered by this: `d_N(r)` from the rest frame applies unchanged at any `γ ≫ 1`. | `\|arg d_1(w0)/d_1(0)\|` = 0.85 rad, `d_2`: 1.69 rad (LG₀⁰, `a0 = 0.1`) |
| `ξ = [β/(1+β)] n_eff λ_L / z_R`, `z_R = k_L w0²/2` | Path through the focus during the interaction, in Rayleigh ranges. Gouy change about `(2p+\|l\|+1) ξ` | 0.071 → 0.21 rad for LG₀² | 0.008 |

### Observation

| Quantity | Unit | Notes |
|---|---|---|
| Angle | `θ_N = 2/(k_N w0) = w0/Z_N` | `k_N` = the dressed harmonic the run evaluates |
| Distance | `Z_N = k_N w0²/2` | Far field for `Z ≫ Z_N` |
| Frequency | `ω_1` | Already the `radiation_field.dat` convention |
| Transverse size, waists | `w0` | Fitted emitted waists in `w0` |
| Energy | photon number (invariant); energies `× ħω_N` | Ratios between harmonics: photon numbers |

`Z_1` for the configs in use (`k_1/k_L = D²/(1 + a0²/2)`, `a0 = 0.1`):

| `γ` | `D` | `k_1/k_L` | `Z_1`, `w0 = 8.4 λ` | `Z_1`, `w0 = 25 λ` | `ε`, 8.4 λ | `ε`, 25 λ |
|---|---|---|---|---|---|---|
| 13.52 (`p_z = 13.483 mc`) | 27.00 | 725.5 | 1.92e4 `w0` | 5.70e4 `w0` | 1.4e-3 | 4.7e-4 |
| 4 | 7.873 | 61.7 | 1628 `w0` | 4844 `w0` | 4.8e-3 | 1.6e-3 |
| 1 (at rest) | 1 | 1 | 26.3 `w0` | 78.2 `w0` | 3.8e-2 | 1.3e-2 |

The current screen (`circular_detector_distance 160000 w0`, `R_max 15 w0`) is at 8.4 `Z_1` and reaches
1.8 `θ_1` for `w0 = 8.4 λ`; for runs A–D (`w0 = 25 λ`) it was 2.8 `Z_1` and 5.3 `θ_1`.

## 3. What the scaling does not cover

- **`γ` near 1.** At rest `β/(1+β) = 0` and `ξ = 0`: an electron at rest doesn't cross the focus. Rest-frame runs are
  equivalent to moving-beam runs only where `ξ` is negligible (`w0 = 25 λ`: yes; `w0 = 8.4 λ`: 0.21 rad for LG₀²,
  probably not).
- **Large `ε`.** Tight focusing in the rest frame (`D k_L w0` of order 10). Not our case.
- **The paraxial laser model.** Its error, `O(θ_d²)` with `θ_d = 2/(k_L w0)` (1.4e-3 at 8.4 λ), is the same in every
  frame and does not break the scaling, but sets how well two "equivalent" runs can agree.
- **Beam momentum spread.** It adds `2 n_eff σ_γ/γ` (dephasing over the line) and the electron angles against
  `θ_N`. Not part of the cold flat-disk results.
- **3D and random beams** add `N_eff = n_⊥ A_eff` (transverse coherent/incoherent ratio; `n_⊥` is the areal density,
  invariant) and `H`, `σ_z` in units of `λ_1 = 2π/|q_axial|` (the code's unit). Those carry the second part of the
  paper without new real numbers.
- **Microbunching** brings in `b_n`, the energy spread and the charge: the real-number part.

## 4. Evidence already in hand

The at-rest lens runs (`future_dev/ponderomotive_lens.md`, Gaussian pulse; `γ = 1`) and the report's runs A and B
(`γ = 13.52`) have the same `a0`, pulse, LG₀⁰ laser with `w0 = 25 λ`, circular polarization and disk radius. They
differ in `γ` by a factor 13.5 and in screen distance (10.2 vs 2.8 `Z_1`); the rest runs use a lattice, A and B a
random beam. `ξ` is below 0.01 in both, so the scaling says the fitted beams should agree in `w0` units and `Z_N`
units.

| `a0` | `N` | waist, rest [`λ_L`] | waist, lab [`λ_L`] | `z0/Z_N`, rest | `z0/Z_N`, lab |
|---|---|---|---|---|---|
| 0.1 | 1 | 21.2 | 21.4 | 0.37 | 0.35 |
| 0.1 | 2 | 10.9 | 10.6 | 0.18 | 0.18 |
| 0.2 | 1 | 8.5 | 8.7 | 0.19 | 0.21 |
| 0.2 | 2 | 4.7 | 4.0 | 0.090 | 0.071 |

Rest values from the lens note (series waist); lab values from report Sec. 4.3, where `z0` is given to one digit.
`Z_N` uses the dressed `k_N` of each run (`k_1/k_L = 1/(1 + a0²/2)` at rest; 725.5 and 714.9 in the lab for
`a0 = 0.1, 0.2`). Agreement is within the lab fits' precision, except `a0 = 0.2`, `N = 2` (15% in the waist), where
the lab fit holds only 80% of the power in one mode. This is not yet a controlled test (different beams, screens and
`z0` precision), which is what §5 is for.

## 5. Proposed collapse test

Two runs that differ only in `γ`, with every scaled parameter equal:

- laser 10 µm, LG₀² (then LG₀⁰), `w0 = 8.4 λ`, `a0 = 0.1`, circular, Gaussian `σ = 25` cycles, cutoff 4;
- flat polar lattice, `R = 2.37 w0` (fixed), spacing 2.5 λ;
- `γ = 13.52` (`average_pz 13.483 mc`): `circular_detector_distance 160000 w0`, `R_max 15 w0`;
- `γ = 4` (`average_pz 3.873 mc`): `circular_detector_distance 13602 w0` (= 160000 × 61.68/725.5), `R_max 15 w0`
  (same `Z/Z_1` and the same reach in `θ_1`); same screen grid.

At `N = 2` the screen reaches half as far in `θ_2`, in both runs alike.

Compare: `|F|` and phase maps on the scaled screen; `l` purity and LG content; fitted waist in `w0` and `z0/Z_N`;
`Λ_zz/P_z`; photon number per electron at each harmonic. Expected agreement: set by `ε²` (about 2e-5 at `γ = 4`) and by the change
of `β/(1+β)` (1.3% of `ξ`, about 3e-3 rad of Gouy phase for LG₀²). Anything larger is a bug or a missing parameter.

Optional third run, `γ = 1` with `circular_detector_distance 220.5 w0`: equal in everything but `ξ` (0 against
0.07). The difference measures how much the Gouy/curvature change during the interaction matters at `w0 = 8.4 λ`.

## 6. Possible solver change (not approved; only if wanted)

Screen geometry is given in `λ`/`w0`, so a config at a new `γ` needs the distance recomputed by hand (§5). Detector
units `Z_1` and `theta_1`, resolved from the laser, beam and detector normal like `lambda_1`, would make a scaled
config portable between energies. Not needed for the collapse test.
