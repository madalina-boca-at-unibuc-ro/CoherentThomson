"""
Computes and plots the screen-level angular-momentum flux DENSITY split into intrinsic, spin, and
orbital-intrinsic parts, per theory/angular_momentum_flux_screen_from_Faraday.md -- distinct from
plot_angular_momentum_flux.py's own (differently normalized, un-split) flux density from
theory/angular_momentum_flux_density.md. Both are Python-only post-processing of an existing run's
radiation_field.dat; no new C++ output is needed.

Scope and caveats (read before trusting a result):
- Same rectangular/circular flat-screen restriction as plot_angular_momentum_flux.py, and for the same
  reason: the theory doc assumes one flat screen parallel to the local x1x2 plane with a single
  well-defined normal, which a spherical detector's points don't generally share.
- The doc's formulas assume the screen's energy-flux centroid sits at the local origin (X_c=Y_c=0),
  which makes the "intrinsic" and "total" J3 flux densities coincide -- no centroid is computed here.
- The screen-level distribution (the doc's Sections 1-4) and its Section 5 area integration into a
  total J3(omega) spectrum are both implemented. The integration is a rough trapezoidal quadrature,
  adapted to whichever flat-screen geometry the run used: (x, y) for 'rectangular', (r^2, phi) for
  'circular' (matching CircularDetector's own equal-area radial spacing) -- see
  get_screen_area_weights_au. Both are exact-in-total-area regardless of grid resolution (a
  trapezoidal rule integrates a constant function exactly), which is a useful sanity check but not a
  guarantee that a coarse grid resolves the flux density's actual spatial variation well. Running this
  script APPENDS a 'J3' summary table to the run's own run_log.txt (replacing any such section from a
  previous run of this script on the same run_log.txt, so re-running doesn't pile up duplicates).
- The spin term's explicit 1/omega factor needs the *raw* angular frequency, not radiation_field.dat's
  own omega column (which is normalized to the fundamental, omega/omega_1 -- see CLAUDE.md's
  radiation_field.dat omega-column bullet). This script reconstructs the raw value by reading the
  run's own run_log.txt for its logged 'fundamental_frequency' (already in raw a.u.,
  Logging::write_run_log/run_log.cpp) and multiplying by the omega/omega_1 ratio -- both are exactly
  proportional (the shared k=omega/c convention cancels between them), so this reconstruction is
  exact, not approximate. Rows with reconstructed omega<=0 are left as NaN for the spin/orbital terms,
  per the theory doc's explicit "do not evaluate at omega=0" instruction.
- The theory doc's boxed formulas assume a specific Fourier normalization (its own "Fourier and tensor
  conventions" section: F(omega) = 1/(2*pi) * int dt e^{i*omega*t} F(t), with the doc's 4*pi
  prefactors already folding in the two-sided-to-one-sided factor of 2). This script applies those
  formulas literally to the Faraday tensor radiation_field.dat already exports; whether that export's
  own normalization (Simulation::run_simulation's general_factor) matches the doc's assumed convention
  exactly has not been independently re-derived here -- flag this before comparing absolute
  magnitudes against an external reference, the same kind of open question as the FT sign-convention
  check elsewhere in this repo.
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
    convert_unit_to_number,
    get_canonical_to_local_rotation,
    get_screen_coordinates_local_au,
    extract_rotated_faraday_fields,
)
from utils.w0_axes_utils import get_laser_lg_w0_in_axes_units, add_w0_secondary_axes

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

# Core::PhysUtils::AtomicUnits (phys_utils.hpp): c and epsilon_0 in the solver's own atomic units --
# c is 1/alpha (2018 CODATA), and 4*pi*epsilon_0=1 defines these units so epsilon_0=1/(4*pi). Mirrors
# phys_utils.hpp's own values independently (no shared constants module between C++ and Python) --
# keep in sync if those change. Matches plot_angular_momentum_flux.py's own constants exactly.
C_LIGHT = 137.035999084
EPSILON_0 = 1.0 / (4.0 * np.pi)

QUANTITIES = (
    ('M33_int', 'intrinsic', r'$d\mathcal{M}^{\mathrm{int}}_{33,+}/d\omega$ (a.u.)'),
    ('M33_spin', 'spin', r'$d\mathcal{M}^{\mathrm{spin}}_{33,+}/d\omega$ (a.u.)'),
    ('M33_orb_int', 'orbital', r'$d\mathcal{M}^{\mathrm{orb,int}}_{33,+}/d\omega$ (a.u.)'),
)


def get_fundamental_frequency_au(run_dir):
    """
    Returns the run's raw fundamental angular frequency omega_1, in atomic units -- NOT the omega/c
    convention Simulation::simulation_parameters::fundamental_frequency stores internally (run_log.cpp
    already multiplies by c before logging it, see Logging::write_run_log). Needed to convert
    radiation_field.dat's own omega column (normalized to omega/omega_1) back into a raw angular
    frequency for the spin flux density's explicit 1/omega factor.
    """
    log_path = os.path.join(run_dir, "run_log.txt")
    if not os.path.exists(log_path):
        raise FileNotFoundError(
            f"'{log_path}' not found -- this run predates run_log.txt (Logging::write_run_log); "
            "re-run the solver to regenerate it.")
    with open(log_path) as f:
        for line in f:
            parts = line.split()
            if parts and parts[0] == 'fundamental_frequency':
                return float(parts[1])
    raise ValueError(f"No 'fundamental_frequency' entry found in '{log_path}'")


def compute_angular_momentum_flux_screen(radiation_filepath):
    """
    Returns (result, detector_type, config_path). `result` has one row per radiation_field.dat row
    (i_omega, omega, i_screen), the local-frame (x_local, y_local) position, and the three flux
    densities (M33_int, M33_spin, M33_orb_int) from
    theory/angular_momentum_flux_screen_from_Faraday.md, evaluated on the physical total field
    (LR+SR+BR, see extract_rotated_faraday_fields -- BR is identically zero for
    radiation_formula="direct", so this is correct regardless of which formula produced the file).
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
            "(x1, x2) plane -- theory/angular_momentum_flux_screen_from_Faraday.md assumes one flat "
            "transverse screen, which only 'rectangular' and 'circular' detectors are.")
    x_local, y_local = get_screen_coordinates_local_au(detector_type, config_path)

    data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    x1 = x_local[data['i_screen'].to_numpy()]
    x2 = y_local[data['i_screen'].to_numpy()]
    omega_raw = data['omega'].to_numpy() * omega_1

    R_to_local = get_canonical_to_local_rotation(config_path)
    f = extract_rotated_faraday_fields(data, R_to_local)

    Ex = f['Ex_l'] + f['Ex_s'] + f['Ex_b']
    Ey = f['Ey_l'] + f['Ey_s'] + f['Ey_b']
    Ez = f['Ez_l'] + f['Ez_s'] + f['Ez_b']
    Bx = f['Bx_l'] + f['Bx_s'] + f['Bx_b']
    By = f['By_l'] + f['By_s'] + f['By_b']
    Bz = f['Bz_l'] + f['Bz_s'] + f['Bz_b']

    # theory doc's component conventions: E_i = c*F^{i0}, B1=-F^{23}, B2=F^{13}, B3=-F^{12}.
    F10, F20, F30 = Ex / C_LIGHT, Ey / C_LIGHT, Ez / C_LIGHT
    F23, F13, F12 = -Bx, By, -Bz

    T13 = -4.0 * np.pi * EPSILON_0 * C_LIGHT**2 * np.real(F10 * np.conj(F30) + F23 * np.conj(F12))
    T23 = 4.0 * np.pi * EPSILON_0 * C_LIGHT**2 * np.real(-F20 * np.conj(F30) + F13 * np.conj(F12))
    M33_int = x1 * T23 - x2 * T13

    valid_omega = omega_raw > 0
    spin_complex = -np.conj(F23) * F10 + np.conj(F13) * F20 + np.conj(F12) * F30
    M33_spin = np.full_like(M33_int, np.nan)
    M33_spin[valid_omega] = (4.0 * np.pi * EPSILON_0 * C_LIGHT**3 / omega_raw[valid_omega]
                             * np.imag(spin_complex[valid_omega]))
    if not np.all(valid_omega):
        print(f"Warning: {np.sum(~valid_omega)} row(s) have omega<=0; spin/orbital flux left as NaN "
              "there (theory/angular_momentum_flux_screen_from_Faraday.md: undefined at omega=0).")

    M33_orb_int = M33_int - M33_spin

    result = data[['i_omega', 'omega', 'i_screen']].copy()
    result['x_local'] = x1
    result['y_local'] = x2
    result['M33_int'] = M33_int
    result['M33_spin'] = M33_spin
    result['M33_orb_int'] = M33_orb_int
    return result, detector_type, config_path


def _linear_trapezoidal_weights(N, a, b):
    """
    Standard trapezoidal quadrature weights for N points evenly spaced over [a, b] (half weight at
    each of the two endpoints, full weight at interior points): sum(weights) == b - a exactly,
    regardless of N, since a trapezoidal rule integrates a constant function exactly. N == 1 has no
    well-defined spacing, so the single point is given the whole interval's weight (b - a) as a
    fallback -- an edge case not expected to occur for the detector types this script supports (a
    true single-point flat detector has no meaningful screen area to integrate).
    """
    if N == 1:
        return np.array([b - a])
    d = (b - a) / (N - 1)
    w = np.full(N, d)
    w[0] *= 0.5
    w[-1] *= 0.5
    return w


def get_screen_area_weights_au(detector_type, config_path):
    """
    Returns an array of per-screen-point quadrature area weights Delta A_p (atomic length^2 units),
    in the same row-major i_screen order radiation_field.dat uses -- the discretization Section 5 of
    theory/angular_momentum_flux_screen_from_Faraday.md needs for
    J3^q(omega) ~= sum_p Delta A_p * M33^q_p.

    'circular': trapezoidal quadrature in (r^2, phi), matching Core::Detector::CircularDetector's own
    equal-area radial spacing (dA = r dr dphi = (1/2) d(r^2) dphi, so a constant-r^2-step trapezoidal
    rule is exact in area for any N_R) and its phi = linspace(0, 2*pi, N_phi) convention (phi=0 and
    phi=2*pi are the same physical direction, both included as the first/last grid points -- the
    ordinary trapezoidal rule's own half-weighting at both endpoints handles this correctly without
    double-counting the shared wedge, since it reduces to the exact periodic-trapezoidal rule for a
    function with equal values at both ends).
    'rectangular': a plain 2D trapezoidal rule in (x, y).

    Row order matches get_screen_coordinates_local_au/get_screen_coordinates's own construction:
    'rectangular' is row-major over (i outer over x, j inner over y); 'circular' is row-major over
    (i outer over r, j inner over phi).
    """
    if detector_type == 'rectangular':
        Nx = int(read_config_value('rectangular_detector_Nx', config_path)[0])
        Ny = int(read_config_value('rectangular_detector_Ny', config_path)[0])
        x_min_s, x_unit = read_config_value('rectangular_detector_x_min', config_path)
        x_max_s, _ = read_config_value('rectangular_detector_x_max', config_path)
        y_min_s, y_unit = read_config_value('rectangular_detector_y_min', config_path)
        y_max_s, _ = read_config_value('rectangular_detector_y_max', config_path)
        x_min = float(x_min_s) * convert_unit_to_number(x_unit, config_path)
        x_max = float(x_max_s) * convert_unit_to_number(x_unit, config_path)
        y_min = float(y_min_s) * convert_unit_to_number(y_unit, config_path)
        y_max = float(y_max_s) * convert_unit_to_number(y_unit, config_path)

        w_x = _linear_trapezoidal_weights(Nx, x_min, x_max)
        w_y = _linear_trapezoidal_weights(Ny, y_min, y_max)
        return (w_x[:, None] * w_y[None, :]).ravel()

    elif detector_type == 'circular':
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        R_min_s, r_unit = read_config_value('circular_detector_R_min', config_path)
        R_max_s, _ = read_config_value('circular_detector_R_max', config_path)
        R_min = float(R_min_s) * convert_unit_to_number(r_unit, config_path)
        R_max = float(R_max_s) * convert_unit_to_number(r_unit, config_path)

        w_r_sq = _linear_trapezoidal_weights(N_R, R_min**2, R_max**2)
        w_phi = _linear_trapezoidal_weights(N_phi, 0.0, 2.0 * np.pi)
        return (0.5 * w_r_sq[:, None] * w_phi[None, :]).ravel()

    else:
        raise ValueError(
            f"detector_type '{detector_type}' is not a flat screen with a single well-defined "
            "local (x1, x2) plane -- only 'rectangular' and 'circular' detectors are supported.")


def integrate_angular_momentum_flux_screen(result, detector_type, config_path):
    """
    Rough numerical area integration of Section 5 of
    theory/angular_momentum_flux_screen_from_Faraday.md: J3^q(omega) = sum_p Delta A_p * M33^q_p,
    q in {int, spin, orb_int}, using get_screen_area_weights_au's quadrature weights. Returns a
    DataFrame with one row per unique i_omega (columns i_omega, omega, J3_int, J3_spin, J3_orb_int),
    sorted by i_omega. A NaN M33_spin/M33_orb_int row (omega<=0, see
    compute_angular_momentum_flux_screen) propagates to NaN in the corresponding J3.
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
            'J3_int': np.sum(area_weights * subset['M33_int'].to_numpy()),
            'J3_spin': np.sum(area_weights * subset['M33_spin'].to_numpy()),
            'J3_orb_int': np.sum(area_weights * subset['M33_orb_int'].to_numpy()),
        })
    return pd.DataFrame(rows).sort_values('i_omega').reset_index(drop=True)


_RUN_LOG_SECTION_HEADER = "Angular-momentum flux screen integration (plot_angular_momentum_flux_screen.py)"

# See plot_angular_momentum_density_screen.py's own _UNCALIBRATED_MAGNITUDE_CAVEAT for the full
# reasoning -- identical issue here: incident_field.dat's field has no Fourier-transform normalization
# applied, so M33_* absolute magnitudes from --incident are on an arbitrary scale, not comparable to a
# real scattered-radiation run's own M33_* values. Only ratios (e.g. M33_int/M33_spin, or against the
# density script's own L3_* values) are meaningful.
_UNCALIBRATED_MAGNITUDE_CAVEAT = (
    "NOTE: incident_field.dat's absolute field magnitude has no Fourier-transform normalization "
    "applied (see this script's _UNCALIBRATED_MAGNITUDE_CAVEAT) -- only RATIOS between the columns "
    "below (e.g. M33_orb_int/M33_spin) are physically meaningful; the absolute values are on an "
    "arbitrary scale and are NOT comparable to a real scattered-radiation run's own M33_* values."
)


def append_integration_to_run_log(run_dir, detector_type, integrated, label_suffix=""):
    """
    Appends (or, on a re-run, replaces) a 'J3' summary table in the run's own run_log.txt --
    Section 5 of theory/angular_momentum_flux_screen_from_Faraday.md integrated over the whole
    screen, one row per frequency. Only this script's own previously-appended section (identified by
    _RUN_LOG_SECTION_HEADER + label_suffix) is ever replaced; everything Logging::write_run_log itself
    wrote is left untouched, so re-running this script doesn't pile up duplicate sections.
    `label_suffix` (e.g. " -- incident beam", set when run on incident_field.dat -- see
    plot_angular_momentum_flux_screen) keeps that table in its own section rather than overwriting the
    actual scattered-radiation run's own logged section.
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
             f"  Rough sum_p Delta_A_p * M33_p screen integration (theory/angular_momentum_flux_screen_from_Faraday.md",
             f"  Section 5), {detector_type} detector, trapezoidal quadrature in "
             f"{'(x, y)' if detector_type == 'rectangular' else '(r^2, phi)'} (exact in total area regardless of grid",
             f"  resolution; see this script's own module docstring for the un-verified FT-normalization caveat).",
             f"  {'i_omega':>8} {'omega/omega_1':>14} {'J3_int [a.u.]':>16} {'J3_spin [a.u.]':>16} "
             f"{'J3_orb_int [a.u.]':>18}"]
    for _, row in integrated.iterrows():
        lines.append(f"  {int(row['i_omega']):>8d} {row['omega']:>14.6f} {row['J3_int']:>16.6e} "
                     f"{row['J3_spin']:>16.6e} {row['J3_orb_int']:>18.6e}")

    with open(log_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Appended angular-momentum flux screen integration to {log_path}")


def _plot_one_quantity(result, detector_type, config_path, radiation_filepath, png_dir, column, name, axis_label,
                       label_suffix="", title_suffix=""):
    """
    Renders `column` (one of M33_int/M33_spin/M33_orb_int) the same way
    plot_angular_momentum_flux.py's own plot_angular_momentum_flux renders its single flux_total: a
    line plot vs. omega/omega_1 for a 1x1 detector (the dense_frequency_spectrum workflow), or one
    heatmap PNG per frequency on the detector's own native grid otherwise. `label_suffix`/
    `title_suffix` (set when run on incident_field.dat -- see plot_angular_momentum_flux_screen) keep
    those PNGs from overwriting the actual scattered-radiation run's own plots.
    """
    if result['i_screen'].nunique() == 1:
        subset = result.sort_values('omega')
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(subset['omega'], subset[column], marker='.')
        ax.set_xlabel("$\\omega / \\omega_1$ (units of the fundamental)")
        ax.set_ylabel(axis_label)
        ax.grid(True)
        fig.suptitle(f"Spectral angular-momentum flux density ({name}){title_suffix}", fontsize=13, fontweight='bold')
        plt.tight_layout()

        output_img = os.path.join(png_dir, f"angular_momentum_screen_{name}{label_suffix}_spectrum.png")
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
        # pad=0.15 (rather than matplotlib's default ~0.05) leaves room for the secondary right
        # y/w0 axis added below -- otherwise the colorbar crowds it into invisibility.
        fig.colorbar(sc, ax=ax, label=axis_label, pad=0.15)
        add_w0_secondary_axes(ax, w0)

        fig.suptitle(f"Angular-momentum flux density ({name}){title_suffix}, $\\omega$ index {i_omega} "
                     f"($\\omega/\\omega_1$={omega_value:.4g}), {detector_type} detector, "
                     f"{detector_geometry_label}", fontsize=11, fontweight='bold')
        plt.tight_layout()

        output_img = os.path.join(png_dir, f"angular_momentum_screen_{name}{label_suffix}_omega{i_omega}.png")
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)


def plot_angular_momentum_flux_screen(radiation_filepath):
    """
    Computes and renders the three screen-level flux densities (intrinsic, spin, orbital-intrinsic)
    from theory/angular_momentum_flux_screen_from_Faraday.md -- three separate plots (or, for a
    multi-point detector, three separate sets of per-frequency plots), one per quantity. Also
    integrates each over the screen (Section 5) and appends the resulting J3(omega) table to the
    run's own run_log.txt (see append_integration_to_run_log).

    `radiation_filepath` may point at either a run's radiation_field.dat or its incident_field.dat
    sibling (Core::Radiation::export_incident_field_fourier) -- see
    plot_angular_momentum_density_screen's own doc comment for why/how; detected by filename so an
    incident-beam run gets its own PNG names and run_log section instead of overwriting the real run's.
    """
    is_incident = os.path.basename(radiation_filepath) == "incident_field.dat"
    label_suffix = "_incident" if is_incident else ""
    title_suffix = " (incident beam)" if is_incident else ""
    run_log_header_suffix = " -- incident beam" if is_incident else ""

    result, detector_type, config_path = compute_angular_momentum_flux_screen(radiation_filepath)

    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder", MODULE_NAME)
    os.makedirs(png_dir, exist_ok=True)

    for column, name, axis_label in QUANTITIES:
        _plot_one_quantity(result, detector_type, config_path, radiation_filepath, png_dir, column, name, axis_label,
                           label_suffix=label_suffix, title_suffix=title_suffix)

    integrated = integrate_angular_momentum_flux_screen(result, detector_type, config_path)
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
        plot_angular_momentum_flux_screen(input_file)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)
