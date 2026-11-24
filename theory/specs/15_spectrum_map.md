# Spectrum along a screen line: dP_z/dω and |E±| maps in (θ, ω)

> **Status: APPROVED 2026-11-24.** Written 2026-11-24. The user asked for the whole implementation and a merge, so the
> checkpoint (§7) is a self-check, not a stop.

## 1. Purpose

Where does the computed ω_N sit in the real spectrum, how wide is the line, and does the spectrum change across the
screen? At a0 = 1 the emission at ω₁ (dressed with the peak a0) is the red edge of the N = 1 line: only electrons at
the LG maximum, at the pulse peak, reach it, and the far field there is the J₂ pattern of that thin ring
(`~/Dropbox/work/AI-outputs/Claude/2026-11-24/screen_distance_and_harmonics_findings.md`, §3). Today the dense
spectrum is computed at one screen point only, so the line can't be followed across the screen.

Wanted: the dense spectrum on a **line of screen points** (one azimuth, a range of θ = ρ/D), shown as maps over
(θ, ω/ω₁) of
- **dP_z/dω**, the spectral energy flux (the quantity that says where the line is and how wide), and
- **|E₊|** and **|E₋|**, the two helicity components (the structure per helicity).

Expected physics, for orientation (not a test): at these angles γθ ~ 10⁻³, so the angle dependence of ω_N is
negligible; a θ–ω correlation comes from the source. At a given ω only electrons whose local a(ρ_beam, t) makes it
resonant emit strongly (one ring in the beam, or two either side of the LG maximum), and a ring of radius ρ₀ puts
its far-field maximum near θ ≈ 3/(kρ₀). The map should show the ring moving and splitting between ω₁ and the
undressed limit (1 + a0²/2) ω₁. This is the spectral image of each electron's chirp, not the chirp in time.

## 2. Why these quantities, and the axisymmetry assumption

One azimuth stands for the whole screen only for quantities that don't depend on φ. For an axisymmetric source
(LG laser, axisymmetric beam, head-on) that holds for |E₊|, |E₋| and dP_z/dω. It does **not** hold for Cartesian
components when both helicities are present: the a0 = 1 N = 2 field had |E_x| ∝ |cos φ|. So the maps use E± and
dP_z/dω, not |F_μν|. The assumption is checked (§6.5).

Definitions, the existing ones:
- E₊ = (E_x − iE_y)/√2, E₋ = (E_x + iE_y)/√2 (`lg_decomposition.helicity_components`: spin ±1 about lab +z for
  e^{−iωt}); total field LR + SR + BR; E_i = c F^{i0}.
- dP_z/dω = 4π ε₀ c² Re(E_x B_y* − E_y B_x*) (`plot_point_spectrum_observables.compute_point_observables`, with the
  2π Parseval factor of 2026-11-21).

## 3. Solver: dense spectrum on a line

`main.cpp` now refuses a dense spectrum on more than one point (error message, then `restrict_to_first_point`). The
radiation loop handles any number of points and frequencies; the frequency recurrence only needs evenly spaced
frequencies, which the dense mode has.

- **Allowed:** one point (as now), or a **line**: rectangular with `rectangular_detector_Ny 1` (or `Nx 1`), circular
  with `circular_detector_N_phi 1`. Recommended: rectangular, `y_min = y_max = 0`, `x_min = 0`, so θ is evenly
  spaced; the circular line has the equal-area radial grid (`R_max/√(N_R − 1)` first ring), sparse near the axis
  where N = 1 peaks.
- **Not allowed:** a 2D grid or a spherical detector with `dense_frequency_spectrum true`: keep today's behaviour
  (message, first point only).
- **Memory check before the run:** per thread `points × N_omega × 288 B` (`PackedFaraday`: 3 bivectors of 6
  complex), plus the reconstructed field (`× 768 B`) and `radiation_field.dat` (~1.3 kB per row as text). Print the
  three in `run_log.txt` ("Frequency spectrum" section) and on stdout; refuse above a cap (e.g. 8 GB total buffers),
  with the numbers in the message. Example: 64 points × 2000 ω: 37 MB per thread (0.9 GB on 24 threads), 98 MB
  field, ~170 MB file.
- `run_log.txt`: the line's geometry (type, θ range and spacing, azimuth) next to `screen points`.

## 4. Plot script: `py_scripts/radiation/plot_spectrum_map.py [run_folder] [--range R] [--omega-min/max]`

Reads a dense-spectrum line run (refuses a single point or a non-line run, pointing to `plot_point_spectrum*.py`).
θ of each point from the local screen coordinates: θ = atan(ρ/D), ρ = √(x² + y²) (shown in µrad). PNGs in
`<run>/png_folder/radiation/` under the usual `MODULE_NAME` rules (`py_scripts/PLOTTING_NOTES.md`).

1. **Map figure:** three panels over (θ, ω/ω₁): dP_z/dω, |E₊|, |E₋|; log colour, floor at max/R (`--range`, default
   1e6); shared axes. Overlays: horizontal lines at the computed harmonics N ω₁ (solid) and at the head-on
   undressed limits N ω₁ (1 + a0²/2) (dashed; exact on axis for head-on, θ dependence negligible here, say so in the
   caption).
2. **Line-integrated spectrum:** dW/dω ≈ Σ_i dP_z/dω(θ_i) 2π ρ_i Δρ_i over the line, assuming axisymmetry (log and
   linear), with the same overlays. This is the beam's spectrum through the disk of the screen.
3. **Summary per harmonic window** (ω/ω₁ in [N − ½, N + ½), as `plot_point_spectrum_observables.py` does), printed and
   written as its own `run_log.txt` section (replaced on a re-run):
   - ω/ω₁ of the integrated spectrum's maximum, and its offset from N;
   - the window's FWHM and the band holding 90% of its energy;
   - the fraction of the window's energy within ±1% of N ω₁;
   - θ of the dP_z/dω maximum at N ω₁ and at the integrated spectrum's peak.

Reuse, don't copy: `compute_point_observables` for dP_z/dω (it already works per row and screen point),
`helicity_components`, `read_radiation_data`, the run-folder resolution of `plot_point_spectrum.py`. If spec 13 has
moved the helpers to `field_components.py` by then, import from there.

The existing `plot_point_spectrum*.py` draw one curve per screen point, unreadable for a line: they refuse a run with
more than 4 screen points, with a hint to use `plot_spectrum_map.py`. Spectra at a fixed θ come from single-point
runs (§8).

## 5. Config

A new reference? No: the reference set is one per beam situation. Add a commented example block to
`config/reference/spectrum_lattice_2d.cfg` (the line variant of its screen), and use a work copy for runs. Suggested
first scan, the a0 = 1 case of the findings: `rectangular`, `Nx 64`, `Ny 1`, x from 0 to 60 w0 at 250000 w0
(0–240 µrad), ω/ω₁ from 0.8 to 3.3 with `N_omega 2000` (step 1.25e-3 ω₁; refine later around a line),
`long_distance_direct_approx`, cutoff 6, `laser_delay 0`, lattice spacing 1.25 λ (~1100 electrons). Time it first on
16 points.

## 6. Checks

1. **Single point unchanged:** `spectrum_lattice_2d.cfg` (shrunk `N_omega`) gives a byte-identical
   `radiation_field.dat` before and after.
2. **Line = points:** two points of a line run against two single-point runs at the same positions: equal (byte-
   identical if the summation order per point is unchanged, else to rounding).
3. **Dense vs harmonic:** at ω = N ω₁ exactly (choose `omega_min/max` and `N_omega` so it is a grid point), the line
   run equals a normal run (`dense_frequency_spectrum false`) on the same line, to the frequency recurrence's
   rounding.
4. **dP_z/dω:** at one point, the map's value equals `plot_point_spectrum_observables.py`'s on the same run.
5. **Axisymmetry:** the same scan on the φ = 90° line (`x = 0`, y from 0): dP_z/dω and |E±| equal to the φ = 0 line to
   the lattice's discreteness; |E_x| not (for N = 2), which is why the maps don't use it.
6. **Memory guard:** a run over the cap is refused with the right numbers; the printed estimate matches the
   process's peak memory within ~20% on a medium run.
7. Physics sanity on the a0 = 0.1 config (narrow lines): the integrated spectrum peaks within the line width of
   N ω₁ for N = 1, 2.

## 7. Checkpoint

Stop after §3 and checks 1–3 and 6; report before the plot script.

## 8. Decisions (user, 2026-11-24)

- **A line only** for now; no 2D grid with a dense spectrum.
- **Fixed-θ spectra come from single-point runs** (`plot_point_spectrum*.py` as now); the map script gives the map
  and the line-integrated spectrum only.
