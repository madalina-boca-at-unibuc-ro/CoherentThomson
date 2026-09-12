"""
Computes and plots the screen-level angular-momentum DENSITY (not flux) split into intrinsic
kinetic, spin, and two independent orbital representations (residual and canonical), per
theory/angular_momentum_density_screen_from_Faraday.md.

This is a sibling theory doc to theory/angular_momentum_flux_screen_from_Faraday.md (implemented in
plot_angular_momentum_flux_screen.py) -- despite the superficially similar `x1*Q2 - x2*Q1` structure
and int/spin/orbital split, this doc's own Section 6 is explicit that these are "integrals of the
angular-momentum density on the screen, not angular-momentum fluxes through the screen": `L3_int`/
`L3_spin` here are built from `G1`/`G2` (Poynting-momentum-density components, using only
`F10`/`F20`/`F30`/`F12`/`F13`) rather than the flux script's `T13`/`T23` (built from
`F10`*`F30`/`F23`*`F12` products) -- genuinely different physical quantities from different formulas,
not expected to numerically agree with the flux script's `M33_*` even on the same run.

Reuses plot_angular_momentum_flux.py's/plot_angular_momentum_flux_screen.py's shared machinery
(rotation, screen coordinates, raw-omega reconstruction, area-quadrature weights) rather than
duplicating it.

Scope and caveats (read before trusting a result) -- inherited from plot_angular_momentum_flux_screen.py
unless noted otherwise:
- Same rectangular/circular flat-screen restriction, same X_c=Y_c=0 centroid assumption, same raw-omega
  reconstruction via run_log.txt's logged fundamental_frequency, same un-verified FT-normalization
  caveat (this script applies the density doc's own boxed formulas literally without independently
  re-deriving whether they match radiation_field.dat's actual export normalization).
- **The two orbital forms are NOT expected to agree pointwise** (the doc's own Section 5): `orb_res`
  (`L3_int - L3_spin`) is an exact pointwise decomposition of the intrinsic density; `orb_can` (the
  canonical azimuthal-derivative density) differs from it by a spatial-divergence term and only
  agrees with `orb_res` after integration over a *complete* screen, and only if the boundary term at
  the screen edge vanishes (not checked here) -- do not use pointwise agreement between `orb_res` and
  `orb_can` as a correctness test; only the doc's own listed checks (reality, `L3_int ==
  L3_spin + L3_orb_res`, formula agreement) are implemented as sanity checks.
- **The canonical orbital density needs the azimuthal derivative of each of the six Faraday
  components**, evaluated on the run's own screen grid, per frequency (unlike the algebraic
  `int`/`spin`/`orb_res` quantities, this makes `orb_can` grid-resolution-dependent, not exact):
  - `circular`: computed directly as a periodic centered difference in `phi` at fixed `r`
    (`_azimuthal_derivative_circular`) -- exact for a periodic function, since
    `Core::Detector::CircularDetector`'s own `phi = linspace(0, 2*pi, N_phi)` grid already samples at
    fixed `r` along `phi`, and `d/dphi` at fixed `rho` *is* the doc's azimuthal derivative directly
    (no `x1`/`x2` combination needed) -- the doc's own Cartesian-grid trick
    (`dphi_F = -x2*dF/dx1 + x1*dF/dx2`) is only needed when the grid isn't already a fixed-radius
    polar one.
  - `rectangular`: built from that Cartesian-grid identity via an FFT-based spectral derivative in
    `x` and `y` (`_azimuthal_derivative_rectangular`/`_fft_derivative`) -- exact for a band-limited
    signal that is periodic over the configured window, and far more accurate than a simple
    finite-difference scheme (`np.gradient`, used here previously) once the field's phase varies by
    an appreciable fraction of a cycle per grid step, which it does whenever the detector isn't deep
    in the far field (a real, not-especially-rare regime for this project -- see CLAUDE.md's Fresnel-
    number bullet). Found by direct comparison on a real run: `np.gradient`'s result and this
    spectral result correlated only ~0.75 with each other and disagreed in magnitude too, on a
    64x64 grid where the geometric/Fresnel phase changed by up to ~1.2 rad per grid step near the
    screen edges -- `np.gradient` was simply not resolving that phase curvature accurately.
    **Caveat inherent to the FFT approach**: it implicitly assumes periodicity over the window, i.e.
    that the field's amplitude has decayed close to zero at the screen edges -- a window too narrow
    for the actual field extent will show edge (Gibbs-type) artifacts. If seen, widen
    `rectangular_detector_x/y_min/max` so the field is genuinely small at the boundary, rather than
    only adding more points at the same width. The doc's own convergence-under-grid-refinement check
    is still not attempted automatically here.
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file
from plot_field import (
    read_config_value,
    get_rectangular_cell_edges,
    get_circular_cell_edges,
    get_detector_geometry_label,
)
from plot_angular_momentum_flux import (
    get_canonical_to_local_rotation,
    get_screen_coordinates_local_au,
    extract_rotated_faraday_fields,
)
from plot_angular_momentum_flux_screen import (
    get_fundamental_frequency_au,
    get_screen_area_weights_au,
)
from utils.w0_axes_utils import get_laser_lg_w0_in_axes_units, add_w0_secondary_axes

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

# Core::PhysUtils::AtomicUnits (phys_utils.hpp): c and epsilon_0 in the solver's own atomic units --
# c is 1/alpha (2018 CODATA), and 4*pi*epsilon_0=1 defines these units so epsilon_0=1/(4*pi). Mirrors
# phys_utils.hpp's own values independently (no shared constants module between C++ and Python) --
# keep in sync if those change. Matches plot_angular_momentum_flux_screen.py's own constants exactly.
C_LIGHT = 137.035999084
EPSILON_0 = 1.0 / (4.0 * np.pi)

QUANTITIES = (
    ('L3_int', 'intrinsic', r'$d\mathcal{L}^{\mathrm{int}}_{3,+}/d\omega$ (a.u.)'),
    ('L3_spin', 'spin', r'$d\mathcal{L}^{\mathrm{spin}}_{3,+}/d\omega$ (a.u.)'),
    ('L3_orb_res', 'orbital_residual', r'$d\mathcal{L}^{\mathrm{orb,res}}_{3,+}/d\omega$ (a.u.)'),
    ('L3_orb_can', 'orbital_canonical', r'$d\mathcal{L}^{\mathrm{orb,can}}_{3,+}/d\omega$ (a.u.)'),
)


def _fft_derivative(values, axis, coord):
    """
    Exact spectral derivative of `values` along `axis`, for a uniformly-spaced coordinate array
    `coord` (`coord[1]-coord[0]` sets the sample spacing). Assumes `values` is periodic over the
    full span of `coord` -- see `_azimuthal_derivative_rectangular`'s own docstring for the caveat
    this implies. A size-1 axis has no defined derivative (a single-column or single-row rectangular
    detector has no spatial variation to differentiate) -- treated as zero.
    """
    N = values.shape[axis]
    if N <= 1:
        return np.zeros_like(values)
    d_coord = coord[1] - coord[0]
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=d_coord)
    shape = [1] * values.ndim
    shape[axis] = N
    k = k.reshape(shape)
    return np.fft.ifft(1j * k * np.fft.fft(values, axis=axis), axis=axis)


def _azimuthal_derivative_rectangular(values, x1_grid, x2_grid):
    """
    dphi_F = -x2*dF/dx1 + x1*dF/dx2 (theory doc's Cartesian-grid identity), via `_fft_derivative` in
    `x1` and `x2` -- see the module docstring's rectangular bullet for why this replaced an earlier
    `np.gradient`-based finite-difference version. `values`/`x1_grid`/`x2_grid` all shape (Nx, Ny),
    meshgrid('ij') convention.
    """
    dF_dx1 = _fft_derivative(values, axis=0, coord=x1_grid[:, 0])
    dF_dx2 = _fft_derivative(values, axis=1, coord=x2_grid[0, :])
    return -x2_grid * dF_dx1 + x1_grid * dF_dx2


def _azimuthal_derivative_circular(values):
    """
    The azimuthal derivative at fixed rho reduces exactly to d/dphi for a grid already sampled at
    fixed r along phi (CircularDetector's own construction) -- computed as a periodic centered
    difference over the N_phi-1 *unique* columns (dropping phi=2*pi, which duplicates phi=0), exact
    for a periodic function and free of the seam-boundary artifact a naive one-sided derivative would
    have there. `values` shape (N_R, N_phi). N_phi<=1 (no azimuthal points to differentiate between)
    is treated as zero.
    """
    N_phi = values.shape[1]
    if N_phi <= 1:
        return np.zeros_like(values)
    unique = values[:, :-1]
    d_phi = (2.0 * np.pi) / unique.shape[1]
    dF_dphi = (np.roll(unique, -1, axis=1) - np.roll(unique, 1, axis=1)) / (2.0 * d_phi)
    return np.concatenate([dF_dphi, dF_dphi[:, :1]], axis=1)


def compute_angular_momentum_density_screen(radiation_filepath):
    """
    Returns (result, detector_type, config_path). `result` has one row per radiation_field.dat row
    (i_omega, omega, i_screen), the local-frame (x_local, y_local) position, and the four
    angular-momentum densities (L3_int, L3_spin, L3_orb_res, L3_orb_can) from
    theory/angular_momentum_density_screen_from_Faraday.md, evaluated on the physical total field
    (LR+SR+BR, see extract_rotated_faraday_fields -- BR is identically zero for
    radiation_formula="direct", so this is correct regardless of which formula produced the file).
    Processed one frequency at a time (unlike the flux script's fully-vectorized computation) since
    L3_orb_can needs each frequency's own field reshaped onto the 2D screen grid for its azimuthal
    derivative.
    """
    config_path = os.path.join(os.path.dirname(radiation_filepath), 'config.cfg')
    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"'{config_path}' not found -- this run predates per-run config snapshots "
            "(Core::IoUtils::copy_config_to_run_directory); re-run the solver to regenerate it.")

    run_dir = os.path.dirname(radiation_filepath)
    omega_1 = get_fundamental_frequency_au(run_dir)

    detector_type = read_config_value('detector_type', config_path)[0]
    if detector_type not in ('rectangular', 'circular'):
        raise ValueError(
            f"detector_type '{detector_type}' is not a flat screen with a single well-defined local "
            "(x1, x2) plane -- theory/angular_momentum_density_screen_from_Faraday.md assumes one "
            "flat transverse screen, which only 'rectangular' and 'circular' detectors are.")
    x_local, y_local = get_screen_coordinates_local_au(detector_type, config_path)

    if detector_type == 'rectangular':
        grid_shape = (int(read_config_value('rectangular_detector_Nx', config_path)[0]),
                      int(read_config_value('rectangular_detector_Ny', config_path)[0]))
    else:
        grid_shape = (int(read_config_value('circular_detector_N_R', config_path)[0]),
                      int(read_config_value('circular_detector_N_phi', config_path)[0]))
    x1_grid = x_local.reshape(grid_shape)
    x2_grid = y_local.reshape(grid_shape)

    data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    R_to_local = get_canonical_to_local_rotation(config_path)

    n_skipped = 0
    rows = []
    for i_omega, subset in data.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_ratio = subset['omega'].iloc[0]
        omega_raw = omega_ratio * omega_1
        i_screen = subset['i_screen'].to_numpy()
        x1 = x_local[i_screen]
        x2 = y_local[i_screen]

        f = extract_rotated_faraday_fields(subset, R_to_local)
        Ex = f['Ex_l'] + f['Ex_s'] + f['Ex_b']
        Ey = f['Ey_l'] + f['Ey_s'] + f['Ey_b']
        Ez = f['Ez_l'] + f['Ez_s'] + f['Ez_b']
        Bx = f['Bx_l'] + f['Bx_s'] + f['Bx_b']
        By = f['By_l'] + f['By_s'] + f['By_b']
        Bz = f['Bz_l'] + f['Bz_s'] + f['Bz_b']

        # theory doc's component conventions: E_i = c*F^{i0}, B1=-F^{23}, B2=F^{13}, B3=-F^{12}.
        F10, F20, F30 = Ex / C_LIGHT, Ey / C_LIGHT, Ez / C_LIGHT
        F23, F13, F12 = -Bx, By, -Bz

        G1 = -4.0 * np.pi * EPSILON_0 * C_LIGHT * np.real(F20 * np.conj(F12) + F30 * np.conj(F13))
        G2 = 4.0 * np.pi * EPSILON_0 * C_LIGHT * np.real(F10 * np.conj(F12) - F30 * np.conj(F23))
        L3_int = x1 * G2 - x2 * G1

        if omega_raw > 0:
            spin_complex = (np.conj(F10) * F20 - np.conj(F20) * F10
                            + np.conj(F13) * F23 - np.conj(F23) * F13)
            L3_spin = (2.0 * np.pi * EPSILON_0 * C_LIGHT**2 / omega_raw) * np.imag(spin_complex)

            derivative = _azimuthal_derivative_rectangular if detector_type == 'rectangular' \
                else _azimuthal_derivative_circular
            dphi = {}
            for name, values in (('F10', F10), ('F20', F20), ('F30', F30),
                                 ('F12', F12), ('F13', F13), ('F23', F23)):
                grid = values.reshape(grid_shape)
                dphi_grid = (derivative(grid, x1_grid, x2_grid) if detector_type == 'rectangular'
                            else derivative(grid))
                dphi[name] = dphi_grid.ravel()

            orb_can_complex = (np.conj(F10) * dphi['F10'] + np.conj(F20) * dphi['F20']
                               + np.conj(F30) * dphi['F30'] + np.conj(F12) * dphi['F12']
                               + np.conj(F13) * dphi['F13'] + np.conj(F23) * dphi['F23'])
            L3_orb_can = (2.0 * np.pi * EPSILON_0 * C_LIGHT**2 / omega_raw) * np.imag(orb_can_complex)
        else:
            n_skipped += 1
            L3_spin = np.full_like(L3_int, np.nan)
            L3_orb_can = np.full_like(L3_int, np.nan)

        L3_orb_res = L3_int - L3_spin

        rows.append(pd.DataFrame({
            'i_omega': i_omega, 'omega': omega_ratio, 'i_screen': i_screen,
            'x_local': x1, 'y_local': x2,
            'L3_int': L3_int, 'L3_spin': L3_spin,
            'L3_orb_res': L3_orb_res, 'L3_orb_can': L3_orb_can,
        }))

    if n_skipped:
        print(f"Warning: {n_skipped} frequency/frequencies have omega<=0; spin/orbital densities left "
              "as NaN there (theory/angular_momentum_density_screen_from_Faraday.md: undefined at "
              "omega=0).")

    result = pd.concat(rows, ignore_index=True)
    return result, detector_type, config_path


def integrate_angular_momentum_density_screen(result, detector_type, config_path):
    """
    Screen integration of Section 6 of theory/angular_momentum_density_screen_from_Faraday.md:
    L^q_3(omega) = sum_p Delta A_p * L3^q_p, q in {int, spin, orb_res, orb_can} -- reuses
    plot_angular_momentum_flux_screen.py's get_screen_area_weights_au (the same physical screen area,
    regardless of whether it's weighting a density or a flux). Returns a DataFrame with one row per
    unique i_omega.
    """
    area_weights = get_screen_area_weights_au(detector_type, config_path)

    rows = []
    for i_omega, subset in result.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        if len(subset) != len(area_weights):
            raise ValueError(
                f"i_omega={i_omega}: {len(subset)} screen points in radiation_field.dat but "
                f"{len(area_weights)} expected from the detector's own configured grid -- geometry "
                "mismatch (e.g. the config was edited after this run).")
        rows.append({
            'i_omega': i_omega,
            'omega': subset['omega'].iloc[0],
            'L3_int': np.sum(area_weights * subset['L3_int'].to_numpy()),
            'L3_spin': np.sum(area_weights * subset['L3_spin'].to_numpy()),
            'L3_orb_res': np.sum(area_weights * subset['L3_orb_res'].to_numpy()),
            'L3_orb_can': np.sum(area_weights * subset['L3_orb_can'].to_numpy()),
        })
    return pd.DataFrame(rows).sort_values('i_omega').reset_index(drop=True)


_RUN_LOG_SECTION_HEADER = "Angular-momentum density screen integration (plot_angular_momentum_density_screen.py)"

# Printed (stdout) and logged (run_log.txt) whenever this script runs on incident_field.dat: the
# incident-beam "phasor" (Radiation::export_incident_field_fourier) has NO Fourier-transform
# normalization applied -- it's LaguerreGaussLaser::complex_amplitude's own raw spatial envelope,
# evaluated at one instant, relying on a global phase common to every component cancelling out of
# every bilinear formula here. That makes RATIOS between these quantities meaningful (orb_can/spin,
# orb_res/spin, this-script's-int/flux-script's-int, ...) but leaves the ABSOLUTE magnitudes on an
# arbitrary scale -- found the hard way when a real run's L3_spin ~ 1e13 a.u. looked nonsensical next
# to an actual scattered-radiation run's own L3_spin ~ 1e-6 a.u.: the two are not on comparable
# footing, by construction, and neither should be read as "the beam's real angular momentum in units
# of hbar." See CLAUDE.md's "analytic incident-field cross-check" note for the ratios that ARE meant
# to be trusted from this mode.
_UNCALIBRATED_MAGNITUDE_CAVEAT = (
    "NOTE: incident_field.dat's absolute field magnitude has no Fourier-transform normalization "
    "applied (see this script's _UNCALIBRATED_MAGNITUDE_CAVEAT) -- only RATIOS between the columns "
    "below (e.g. L3_orb_can/L3_spin, L3_orb_res/L3_spin) are physically meaningful; the absolute "
    "values are on an arbitrary scale and are NOT comparable to a real scattered-radiation run's own "
    "L3_* values."
)


def append_integration_to_run_log(run_dir, detector_type, integrated, label_suffix=""):
    """
    Appends (or, on a re-run, replaces) an 'L3' summary table in the run's own run_log.txt --
    Section 6 of theory/angular_momentum_density_screen_from_Faraday.md integrated over the whole
    screen, one row per frequency. Only this script's own previously-appended section (identified by
    _RUN_LOG_SECTION_HEADER + label_suffix, distinct from plot_angular_momentum_flux_screen.py's own
    header) is ever replaced. `label_suffix` (e.g. " -- incident beam", set when run on
    incident_field.dat -- see plot_angular_momentum_density_screen) keeps that table in its own
    section rather than overwriting the actual scattered-radiation run's own logged section, and also
    triggers _UNCALIBRATED_MAGNITUDE_CAVEAT being written into that section.
    """
    header = _RUN_LOG_SECTION_HEADER + label_suffix
    log_path = os.path.join(run_dir, "run_log.txt")
    if not os.path.exists(log_path):
        raise FileNotFoundError(
            f"'{log_path}' not found -- this run predates run_log.txt (Logging::write_run_log); "
            "re-run the solver to regenerate it.")

    with open(log_path) as f:
        content = f.read()

    marker = f"\n{header}\n"
    marker_index = content.find(marker)
    if marker_index != -1:
        content = content[:marker_index]

    lines = [content.rstrip("\n"), "", header, "-" * len(header)]
    if label_suffix:
        lines += [f"  {_UNCALIBRATED_MAGNITUDE_CAVEAT}", ""]
    lines += [
             f"  Rough sum_p Delta_A_p * L3_p screen integration (theory/angular_momentum_density_screen_from_Faraday.md",
             f"  Section 6), {detector_type} detector, trapezoidal quadrature in "
             f"{'(x, y)' if detector_type == 'rectangular' else '(r^2, phi)'} (same quadrature as "
             f"plot_angular_momentum_flux_screen.py's own J3 table; see this script's own module",
             f"  docstring for the un-verified FT-normalization caveat and the orb_res-vs-orb_can caveat).",
             f"  {'i_omega':>8} {'omega/omega_1':>14} {'L3_int [a.u.]':>16} {'L3_spin [a.u.]':>16} "
             f"{'L3_orb_res [a.u.]':>18} {'L3_orb_can [a.u.]':>18}"]
    for _, row in integrated.iterrows():
        lines.append(f"  {int(row['i_omega']):>8d} {row['omega']:>14.6f} {row['L3_int']:>16.6e} "
                     f"{row['L3_spin']:>16.6e} {row['L3_orb_res']:>18.6e} {row['L3_orb_can']:>18.6e}")

    with open(log_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Appended angular-momentum density screen integration to {log_path}")


def _plot_one_quantity(result, detector_type, config_path, radiation_filepath, png_dir, column, name, axis_label,
                       label_suffix="", title_suffix=""):
    """
    Renders `column` (one of L3_int/L3_spin/L3_orb_res/L3_orb_can): a line plot vs. omega/omega_1 for
    a 1x1 detector (the dense_frequency_spectrum workflow), or one heatmap PNG per frequency on the
    detector's own native grid otherwise -- same rendering convention as
    plot_angular_momentum_flux_screen.py's own _plot_one_quantity. `label_suffix`/`title_suffix` (set
    when run on incident_field.dat -- see plot_angular_momentum_density_screen) keep those PNGs from
    overwriting the actual scattered-radiation run's own plots and mark the figure titles accordingly.
    """
    if result['i_screen'].nunique() == 1:
        subset = result.sort_values('omega')
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(subset['omega'], subset[column], marker='.')
        ax.set_xlabel("$\\omega / \\omega_1$ (units of the fundamental)")
        ax.set_ylabel(axis_label)
        ax.grid(True)
        fig.suptitle(f"Spectral angular-momentum density ({name}){title_suffix}", fontsize=13, fontweight='bold')
        plt.tight_layout()

        output_img = os.path.join(png_dir, f"angular_momentum_density_{name}{label_suffix}_spectrum.png")
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)
        return

    if detector_type == 'rectangular':
        Nx = int(read_config_value('rectangular_detector_Nx', config_path)[0])
        Ny = int(read_config_value('rectangular_detector_Ny', config_path)[0])
        grid_shape = (Nx, Ny)
        x_edges, y_edges = get_rectangular_cell_edges(config_path)
        axes_unit = read_config_value('rectangular_detector_x_min', config_path)[1]
    else:
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        grid_shape = (N_R, N_phi)
        x_edges, y_edges = get_circular_cell_edges(config_path)
        axes_unit = read_config_value('circular_detector_R_min', config_path)[1]

    x_label = "$x$" + (f" [{axes_unit}]" if axes_unit else "")
    y_label = "$y$" + (f" [{axes_unit}]" if axes_unit else "")
    detector_geometry_label = get_detector_geometry_label(detector_type, config_path)
    w0 = get_laser_lg_w0_in_axes_units(radiation_filepath, axes_unit)

    for i_omega, subset in result.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_value = subset['omega'].iloc[0]

        values = subset[column].to_numpy().reshape(grid_shape)
        finite = values[np.isfinite(values)]
        vmax = np.max(np.abs(finite)) if finite.size and np.any(finite) else 1.0

        fig, ax = plt.subplots(figsize=(7, 6))
        sc = ax.pcolormesh(x_edges, y_edges, values, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        fig.colorbar(sc, ax=ax, label=axis_label, pad=0.15)
        add_w0_secondary_axes(ax, w0)

        fig.suptitle(f"Angular-momentum density ({name}){title_suffix}, $\\omega$ index {i_omega} "
                     f"($\\omega/\\omega_1$={omega_value:.4g}), {detector_type} detector, "
                     f"{detector_geometry_label}", fontsize=11, fontweight='bold')
        plt.tight_layout()

        output_img = os.path.join(png_dir, f"angular_momentum_density_{name}{label_suffix}_omega{i_omega}.png")
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)


def plot_angular_momentum_density_screen(radiation_filepath):
    """
    Computes and renders the four screen-level angular-momentum densities (intrinsic, spin,
    orbital-residual, orbital-canonical) from theory/angular_momentum_density_screen_from_Faraday.md.
    Also integrates each over the screen (Section 6) and appends the resulting L3(omega) table to
    the run's own run_log.txt (see append_integration_to_run_log).

    `radiation_filepath` may point at either a run's radiation_field.dat (the actual coherently-summed
    scattered field) or its incident_field.dat sibling (Core::Radiation::export_incident_field_fourier
    -- the incident LG/plane-wave beam's own analytic field, evaluated at z=0, same screen geometry) --
    both share the identical column format, so every formula/plot here is agnostic to which produced
    it. Detected by filename so incident-beam runs get their own PNG names and run_log section instead
    of overwriting the real run's: see the CLAUDE.md "analytic incident-field cross-check" note for
    why this is useful (a clean, exactly-known reference field with no coherent-sum noise or far-field
    aliasing, to validate these formulas/derivatives against independently of the scattered field's
    own numerical behavior).
    """
    is_incident = os.path.basename(radiation_filepath) == "incident_field.dat"
    label_suffix = "_incident" if is_incident else ""
    title_suffix = " (incident beam)" if is_incident else ""
    run_log_header_suffix = " -- incident beam" if is_incident else ""

    result, detector_type, config_path = compute_angular_momentum_density_screen(radiation_filepath)

    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder", MODULE_NAME)
    os.makedirs(png_dir, exist_ok=True)

    for column, name, axis_label in QUANTITIES:
        _plot_one_quantity(result, detector_type, config_path, radiation_filepath, png_dir, column, name, axis_label,
                           label_suffix=label_suffix, title_suffix=title_suffix)

    integrated = integrate_angular_momentum_density_screen(result, detector_type, config_path)
    if is_incident:
        print(_UNCALIBRATED_MAGNITUDE_CAVEAT)
    print(integrated.to_string(index=False))
    append_integration_to_run_log(os.path.dirname(radiation_filepath), detector_type, integrated,
                                  label_suffix=run_log_header_suffix)


if __name__ == "__main__":
    args = sys.argv[1:]
    use_incident = "--incident" in args
    if use_incident:
        args.remove("--incident")
    field_filename = "incident_field.dat" if use_incident else "radiation_field.dat"

    if len(args) >= 1:
        input_file = os.path.join(args[0], field_filename)
        if not os.path.exists(input_file):
            print(f"Error: '{input_file}' not found")
            sys.exit(1)
    else:
        try:
            input_file = find_latest_output_file(field_filename)
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            sys.exit(1)
        print(f"No folder given; using latest run: {input_file}")

    try:
        plot_angular_momentum_density_screen(input_file)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)
