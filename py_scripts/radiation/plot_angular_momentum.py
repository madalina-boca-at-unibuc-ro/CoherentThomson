"""
Computes and plots the z-component of the spin (SAM), orbital (OAM), and total (TAM) angular
momentum spectral densities from a run's coherently-summed Faraday tensor (radiation_field.dat) or
the incident beam's own analytic field (incident_field.dat, via --incident) -- see
theory/numerical_calculation_of_angular_momentum.md ("Applications for the angular momentum
density") for the derivation. Supersedes plot_spin_angular_momentum.py, which only computed S_z and
(per that doc's own correction) was missing the epsilon_0 = 1/(4*pi) prefactor -- see
Core::PhysUtils::AtomicUnits::epsilon_0 (phys_utils.hpp) for the same constant on the C++ side.

$$\\frac{d\\mathcal{S}_z}{d\\omega} = \\frac{4\\epsilon_0}{\\omega}\\mathrm{Im}[\\tilde E_x^*\\tilde E_y]$$
$$\\frac{d\\mathcal{L}_z}{d\\omega} = \\frac{2\\epsilon_0}{\\omega}\\sum_i\\mathrm{Im}\\left[\\tilde E_i^*\\,\\hat L_z\\tilde E_i\\right],
\\quad \\hat L_z = x\\partial_y - y\\partial_x \\;(\\text{rectangular}) = \\partial_\\phi \\;(\\text{circular})$$
$$\\frac{d\\mathcal{J}_z}{d\\omega} = \\frac{d\\mathcal{L}_z}{d\\omega} + \\frac{d\\mathcal{S}_z}{d\\omega}$$

**Always uses the total (physical) Faraday tensor, LR+SR+BR summed** -- there is no
long-range/short-range/boundary decomposition of an angular-momentum observable (that split only
has meaning for the derivation's intermediate F^{mu nu} terms, see
theory/FT_Faraday_tensor-direct_and_simplified_forms.md); only the sum is the actual radiated field
E appearing in the S_z/L_z formulas above. BR is treated as zero if the file predates the boundary
term; for incident_field.dat, SR/BR are already identically zero so the sum reduces to LR alone.

Ex/Ey/Ez are recovered directly from the stored F^{10}/F^{20}/F^{30} columns (E_i=c*F^{i0}, the same
F^{mu nu}<->E/B convention faraday_frame_utils.py's extract_rotated_faraday_fields documents).

**Restricted to rectangular and circular detectors** (theory doc's own "Target Screens"
instruction): a spherical detector's (theta, phi) grid has no flat local (x, y) plane for the OAM
operator's x*d/dy - y*d/dx (rectangular) or d/dphi (circular) construction to act on, so this script
rejects detector_type=spherical outright -- unlike the old, now-superseded
plot_spin_angular_momentum.py, which did support spherical for S_z alone (no spatial derivative
needed there).

**OAM spatial derivatives**: rectangular uses numpy.gradient in x/y (central differences, taking the
grid spacing into account, per the theory doc's own numerical-implementation instructions); circular
uses an exact periodic centered difference in phi (numpy.gradient has no periodic-boundary mode, so
this is implemented by hand -- see _periodic_phi_derivative), which is what x*d/dy - y*d/dx reduces
to exactly at fixed r on a polar grid. Both are finite-difference approximations, not exact for a
rapidly-varying near-field phase (the same caveat this project documents elsewhere for
np.gradient-based derivatives, e.g. CLAUDE.md's Fresnel-number bullet) -- a coarse grid close to the
beam (small Fresnel number) may under-resolve the true spatial variation. Note the L_z density itself
is invariant under an isotropic rescaling of x/y (x*d/dy-y*d/dx is the rotation generator, scale-free
under x->c*x, y->c*y), so it does not matter that get_screen_coordinates returns x/y in the config's
own length unit (e.g. 'lambda') rather than atomic units -- only the screen-integrated total below
needs an explicit unit conversion, since an area DOES scale with the unit choice.

**Frequency convention**: neither .dat file exports raw atomic-unit omega, only the dimensionless
omega/omega_1 ratio (see CLAUDE.md's radiation_field.dat 'omega' column bullet) -- recovering true
fundamental_frequency in Python would mean re-parsing it out of run_log.txt's human-readable dump.
So the 1/omega prefactor above uses that same ratio, making the plotted densities
d.../d(omega/omega_1) rather than d.../domega in raw atomic units -- consistent with every other
frequency-axis quantity this project already plots in omega/omega_1 units, and differing from the
true d.../domega only by the constant factor omega_1 (fundamental_frequency), which does not affect
the spatial pattern at any given omega index or the relative comparison across harmonics.

**Screen integration**: appends an S_z/L_z/J_z(omega) summary table to the run's own run_log.txt
(theory doc's "Surface Integration" step), replacing any such table from a previous run of this
script on the same run_log.txt so re-running doesn't pile up duplicates. Trapezoidal quadrature in
(x, y) for rectangular, (r^2, phi) for circular (matching CircularDetector's own equal-area radial
spacing, so a constant-r^2-step trapezoidal rule is exact-in-area for any N_R) -- area weights are
converted to atomic-unit length^2 (get_screen_area_weights_au), unlike the density values themselves.

Reuses plot_field.py's screen-geometry/cell-edge reconstruction directly (same CLI shape:
[path_to_run_folder] [--incident]), and the same is_incident split into
png_folder/radiation/{emitted,incident}/ subfolders, so a real run's plots and the incident beam's
never collide and sit alongside the Faraday-tensor component plots.
"""
import sys
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file
from plot_field import (
    read_config_value,
    get_screen_coordinates,
    get_rectangular_cell_edges,
    get_circular_cell_edges,
    get_detector_geometry_label,
)
from faraday_frame_utils import C_LIGHT, convert_unit_to_number
from utils.w0_axes_utils import get_laser_lg_w0_in_axes_units, add_w0_secondary_axes

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

# Core::PhysUtils::AtomicUnits::epsilon_0 (phys_utils.hpp): 4*pi*epsilon_0=1 defines these atomic
# units, so epsilon_0=1/(4*pi). Mirrors the C++ constant independently (no shared constants module
# between C++ and Python, see CLAUDE.md) -- keep in sync if that one changes.
EPSILON_0 = 1.0 / (4.0 * np.pi)

# Printed (stdout) and logged (run_log.txt) whenever this script runs on incident_field.dat: the
# incident-beam "phasor" (Radiation::export_incident_field_fourier) has no Fourier-transform
# normalization applied -- it's LaguerreGaussLaser::complex_amplitude's own raw spatial envelope,
# evaluated at one instant, relying on a global phase common to every component cancelling out of
# every bilinear formula here. That makes ratios between these quantities meaningful (Lz/Sz, ...) but
# leaves the absolute magnitudes on an arbitrary scale -- see CLAUDE.md's "analytic incident-field
# cross-check" note for the ratios that ARE meant to be trusted from this mode.
_UNCALIBRATED_MAGNITUDE_CAVEAT = (
    "NOTE: incident_field.dat's absolute field magnitude has no Fourier-transform normalization "
    "applied -- only RATIOS between the S_z/L_z/J_z columns below are physically meaningful; the "
    "absolute values are on an arbitrary scale and are NOT comparable to a real scattered-radiation "
    "run's own values."
)

QUANTITIES = (
    ('dSz', 'spin', r'$d\mathcal{S}_z/d(\omega/\omega_1)$'),
    ('dLz', 'orbital', r'$d\mathcal{L}_z/d(\omega/\omega_1)$'),
    ('dJz', 'total', r'$d\mathcal{J}_z/d(\omega/\omega_1)$'),
)


def _complex_column(data, prefix, mu, nu):
    return (data[f'{prefix}_F{mu}{nu}_re'] + 1j * data[f'{prefix}_F{mu}{nu}_im']).to_numpy()


def _get_E_total(data):
    """
    Returns (Ex, Ey, Ez), each a complex numpy array (one entry per row of `data`), from the total
    (LR+SR+BR) Faraday tensor -- see this module's own doc comment for why only the total, never a
    single range, is used here. E_i = c*F^{i0}.
    """
    def component(mu):
        value = _complex_column(data, 'LR', mu, 0) + _complex_column(data, 'SR', mu, 0)
        if {f'BR_F{mu}0_re', f'BR_F{mu}0_im'}.issubset(data.columns):
            value = value + _complex_column(data, 'BR', mu, 0)
        return C_LIGHT * value

    return component(1), component(2), component(3)


def _periodic_phi_derivative(grid):
    """
    Exact periodic centered difference along axis=1 (phi), for `grid` shape (N_R, N_phi) sampled on
    phi = linspace(0, 2*pi, N_phi) -- CircularDetector's own convention, where phi=0 and phi=2*pi are
    duplicated as the first/last columns (the same physical direction). Computed over the N_phi-1
    *unique* columns (dropping the duplicate phi=2*pi column) so the derivative is exact for a
    periodic function and free of the one-sided-derivative artifact a naive numpy.gradient call would
    have right at that seam; the duplicate last column gets the same value as the first, consistent
    with phi=0 and phi=2*pi being the same point. N_phi<=1 (nothing to differentiate) is zero.
    """
    N_phi = grid.shape[1]
    if N_phi <= 1:
        return np.zeros_like(grid)
    unique = grid[:, :-1]
    d_phi = (2.0 * np.pi) / unique.shape[1]
    d_unique = (np.roll(unique, -1, axis=1) - np.roll(unique, 1, axis=1)) / (2.0 * d_phi)
    return np.concatenate([d_unique, d_unique[:, :1]], axis=1)


def _rectangular_oam_terms(Ex, Ey, Ez, x_grid, y_grid, x_vals, y_vals):
    """
    sum_i Im[E_i^* (x*dE_i/dy - y*dE_i/dx)] over a (Nx, Ny) grid (Eq. primary-definition-of-Lz-
    rectangular), via numpy.gradient in x (axis 0) and y (axis 1). A size-1 axis has no defined
    derivative (a single-row/column rectangular detector has no spatial variation to differentiate
    along that axis) -- treated as zero, matching numpy.gradient's own minimum-2-points requirement.
    """
    def term(E):
        dE_dx = np.gradient(E, x_vals, axis=0) if E.shape[0] > 1 else np.zeros_like(E)
        dE_dy = np.gradient(E, y_vals, axis=1) if E.shape[1] > 1 else np.zeros_like(E)
        return np.imag(np.conj(E) * (x_grid * dE_dy - y_grid * dE_dx))

    return term(Ex) + term(Ey) + term(Ez)


def _circular_oam_terms(Ex, Ey, Ez):
    """sum_i Im[E_i^* dE_i/dphi] over a (N_R, N_phi) grid (Eq. primary-definition-of-Lz-circular)."""
    def term(E):
        return np.imag(np.conj(E) * _periodic_phi_derivative(E))

    return term(Ex) + term(Ey) + term(Ez)


def compute_angular_momentum_density(radiation_filepath):
    """
    Returns (result, detector_type, config_path). `result` has one row per radiation_field.dat row
    with omega/omega_1 != 0 (i_omega, omega, i_screen, dSz, dLz, dJz), computed from the physical
    total field (LR+SR+BR, see _get_E_total).
    """
    try:
        data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    except Exception as e:
        print(f"Error reading file '{radiation_filepath}': {e}")
        sys.exit(1)

    required = {f'{p}_F{m}0_{part}' for p in ('LR', 'SR') for m in (1, 2, 3) for part in ('re', 'im')}
    if not required.issubset(data.columns):
        print(f"Error: File must contain columns {sorted(required)}")
        sys.exit(1)

    config_path = os.path.join(os.path.dirname(radiation_filepath), 'config.cfg')
    if not os.path.exists(config_path):
        print(f"Error: '{config_path}' not found -- this run predates per-run config snapshots "
              "(Core::IoUtils::copy_config_to_run_directory); re-run the solver to regenerate it.")
        sys.exit(1)

    detector_type, _ = read_config_value('detector_type', config_path)
    if detector_type not in ('rectangular', 'circular'):
        raise ValueError(
            f"detector_type '{detector_type}' is not a flat rectangular/circular screen -- "
            "theory/numerical_calculation_of_angular_momentum.md's angular-momentum density "
            "calculation is only implemented for those two (its own 'Target Screens' instruction).")

    if detector_type == 'rectangular':
        Nx = int(read_config_value('rectangular_detector_Nx', config_path)[0])
        Ny = int(read_config_value('rectangular_detector_Ny', config_path)[0])
        grid_shape = (Nx, Ny)
        x_centers, y_centers, _, _ = get_screen_coordinates(radiation_filepath, config_path)
        x_grid = x_centers.reshape(grid_shape)
        y_grid = y_centers.reshape(grid_shape)
        x_vals, y_vals = x_grid[:, 0], y_grid[0, :]
    else:
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        grid_shape = (N_R, N_phi)

    rows = []
    n_skipped = 0
    for i_omega, subset in data.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_ratio = subset['omega'].iloc[0]
        if omega_ratio == 0.0:
            n_skipped += 1
            continue

        Ex, Ey, Ez = _get_E_total(subset)
        dSz = (4.0 * EPSILON_0 / omega_ratio) * np.imag(np.conj(Ex) * Ey)

        Ex_grid, Ey_grid, Ez_grid = Ex.reshape(grid_shape), Ey.reshape(grid_shape), Ez.reshape(grid_shape)
        if detector_type == 'rectangular':
            oam_terms = _rectangular_oam_terms(Ex_grid, Ey_grid, Ez_grid, x_grid, y_grid, x_vals, y_vals)
        else:
            oam_terms = _circular_oam_terms(Ex_grid, Ey_grid, Ez_grid)
        dLz = (2.0 * EPSILON_0 / omega_ratio) * oam_terms.ravel()

        rows.append(pd.DataFrame({
            'i_omega': i_omega, 'omega': omega_ratio, 'i_screen': subset['i_screen'].to_numpy(),
            'dSz': dSz, 'dLz': dLz, 'dJz': dSz + dLz,
        }))

    if n_skipped:
        print(f"Warning: skipped {n_skipped} frequency/frequencies with omega/omega_1==0 "
              "(division by zero in the 1/omega prefactor).")
    if not rows:
        print("Error: no frequency with omega/omega_1 != 0 -- nothing to compute.")
        sys.exit(1)

    result = pd.concat(rows, ignore_index=True)
    return result, detector_type, config_path


def _linear_trapezoidal_weights(N, lo, hi):
    """Composite trapezoidal quadrature weights for N uniformly-spaced points spanning [lo, hi]."""
    if N <= 1:
        return np.full(N, hi - lo)
    step = (hi - lo) / (N - 1)
    weights = np.full(N, step)
    weights[0] *= 0.5
    weights[-1] *= 0.5
    return weights


def get_screen_area_weights_au(detector_type, config_path):
    """
    Returns an array of per-screen-point quadrature area weights Delta A_p (atomic length^2 units),
    in the same row-major i_screen order radiation_field.dat uses. 'rectangular': plain 2D
    trapezoidal rule in (x, y). 'circular': trapezoidal in (r^2, phi), matching CircularDetector's
    own equal-area radial spacing (dA = r dr dphi = (1/2) d(r^2) dphi) and its
    phi = linspace(0, 2*pi, N_phi) convention (the shared phi=0/2*pi wedge is handled correctly by
    the trapezoidal rule's own half-weighting at both endpoints, without double-counting).
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

    N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
    N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
    R_min_s, r_unit = read_config_value('circular_detector_R_min', config_path)
    R_max_s, _ = read_config_value('circular_detector_R_max', config_path)
    R_min = float(R_min_s) * convert_unit_to_number(r_unit, config_path)
    R_max = float(R_max_s) * convert_unit_to_number(r_unit, config_path)

    w_r_sq = _linear_trapezoidal_weights(N_R, R_min**2, R_max**2)
    w_phi = _linear_trapezoidal_weights(N_phi, 0.0, 2.0 * np.pi)
    return (0.5 * w_r_sq[:, None] * w_phi[None, :]).ravel()


def integrate_angular_momentum_density(result, detector_type, config_path):
    """
    Screen integration (theory doc's 'Surface Integration' step): Q_z(omega) = sum_p Delta A_p *
    (dQ_z/domega)_p, Q in {S, L, J}. Returns a DataFrame with one row per unique i_omega.
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
            'Sz_total': np.sum(area_weights * subset['dSz'].to_numpy()),
            'Lz_total': np.sum(area_weights * subset['dLz'].to_numpy()),
            'Jz_total': np.sum(area_weights * subset['dJz'].to_numpy()),
        })
    return pd.DataFrame(rows).sort_values('i_omega').reset_index(drop=True)


_RUN_LOG_SECTION_HEADER = "Angular-momentum screen integration (plot_angular_momentum.py)"

# Matches Logging::write_run_log's own section-header convention: a title line immediately followed
# by a dash-underline exactly as long as the title (see run_log.cpp/run_log.txt's own "Frequency
# spectrum\n------------------\n"-style sections, and this script's own appended ones below).
_SECTION_BOUNDARY_RE = re.compile(r'(?m)^([^\n]+)\n(-+)\n')


def _find_run_log_sections(content):
    """
    Returns [(start_index, title), ...] for every section header in run_log.txt (see
    _SECTION_BOUNDARY_RE) -- used to splice out only this script's own previously-appended section
    (identified by its exact header text) without disturbing any other section, including a
    *different* label_suffix of this same script (e.g. the '-- incident beam' variant) that may
    follow it. A naive "truncate everything after the first marker match" (this function's own
    previous implementation) would silently delete such a later, differently-labeled section too --
    confirmed by running this script on radiation_field.dat, then incident_field.dat, then
    radiation_field.dat again, which wiped the middle (incident) section outright.
    """
    return [(m.start(), m.group(1)) for m in _SECTION_BOUNDARY_RE.finditer(content)
            if len(m.group(2)) == len(m.group(1))]


def append_integration_to_run_log(run_dir, detector_type, integrated, label_suffix=""):
    """
    Appends (or, on a re-run, replaces) an S_z/L_z/J_z summary table in the run's own run_log.txt --
    theory/numerical_calculation_of_angular_momentum.md's 'Surface Integration' step, integrated over
    the whole screen, one row per frequency. Only this script's own previously-appended section
    (identified by _RUN_LOG_SECTION_HEADER + label_suffix, via _find_run_log_sections) is ever
    replaced -- in place, regardless of what other sections precede or follow it. `label_suffix`
    (e.g. " -- incident beam", set when run on incident_field.dat) keeps that table in its own section
    rather than overwriting the actual scattered-radiation run's own logged section, and also triggers
    _UNCALIBRATED_MAGNITUDE_CAVEAT being written into that section.
    """
    header = _RUN_LOG_SECTION_HEADER + label_suffix
    log_path = os.path.join(run_dir, "run_log.txt")
    if not os.path.exists(log_path):
        raise FileNotFoundError(
            f"'{log_path}' not found -- this run predates run_log.txt (Logging::write_run_log); "
            "re-run the solver to regenerate it.")

    with open(log_path) as f:
        content = f.read()

    sections = _find_run_log_sections(content)
    for i, (start, title) in enumerate(sections):
        if title == header:
            end = sections[i + 1][0] if i + 1 < len(sections) else len(content)
            content = content[:start] + content[end:]
            break

    lines = [content.rstrip("\n"), "", header, "-" * len(header)]
    if label_suffix:
        lines += [f"  {_UNCALIBRATED_MAGNITUDE_CAVEAT}", ""]
    lines += [
        "  Rough sum_p Delta_A_p * (dQ_z/domega)_p screen integration "
        "(theory/numerical_calculation_of_angular_momentum.md,",
        f"  'Surface Integration' step), {detector_type} detector, trapezoidal quadrature in "
        f"{'(x, y)' if detector_type == 'rectangular' else '(r^2, phi)'}, area weights in atomic units.",
        f"  {'i_omega':>8} {'omega/omega_1':>14} {'S_z [a.u.]':>16} {'L_z [a.u.]':>16} {'J_z [a.u.]':>16}",
    ]
    for _, row in integrated.iterrows():
        lines.append(f"  {int(row['i_omega']):>8d} {row['omega']:>14.6f} {row['Sz_total']:>16.6e} "
                     f"{row['Lz_total']:>16.6e} {row['Jz_total']:>16.6e}")

    with open(log_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Appended angular-momentum screen integration to {log_path}")


def plot_angular_momentum(radiation_filepath):
    """
    Plots dS_z/d(omega/omega_1), dL_z/d(omega/omega_1), and dJ_z/d(omega/omega_1) (one single-panel
    diverging-colormap heatmap each, since all three are real-valued) over the detector screen, one
    figure set per configured frequency, then appends the screen-integrated totals to run_log.txt.

    `radiation_filepath` may point at either a run's radiation_field.dat or its incident_field.dat
    sibling, detected by filename -- see this module's own doc comment for how each is handled.
    """
    is_incident = os.path.basename(radiation_filepath) == "incident_field.dat"
    title_suffix = " (incident beam)" if is_incident else ""
    if is_incident:
        print(_UNCALIBRATED_MAGNITUDE_CAVEAT)

    result, detector_type, config_path = compute_angular_momentum_density(radiation_filepath)

    _, _, x_label, y_label = get_screen_coordinates(radiation_filepath, config_path)
    detector_geometry_label = get_detector_geometry_label(detector_type, config_path)

    aspect = 'equal'
    if detector_type == 'rectangular':
        Nx = int(read_config_value('rectangular_detector_Nx', config_path)[0])
        Ny = int(read_config_value('rectangular_detector_Ny', config_path)[0])
        grid_shape = (Nx, Ny)
        x, y = get_rectangular_cell_edges(config_path)
        axes_unit = read_config_value('rectangular_detector_x_min', config_path)[1]
    else:
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        grid_shape = (N_R, N_phi)
        x, y = get_circular_cell_edges(config_path)
        axes_unit = read_config_value('circular_detector_R_min', config_path)[1]

    w0 = get_laser_lg_w0_in_axes_units(radiation_filepath, axes_unit) if axes_unit else None

    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder", MODULE_NAME,
                           "incident" if is_incident else "emitted")
    os.makedirs(png_dir, exist_ok=True)

    for i_omega, subset in result.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_ratio = subset['omega'].iloc[0]

        for column, tag, symbol_label in QUANTITIES:
            values = subset[column].to_numpy()
            abs_max = np.abs(values).max() if values.size else None

            fig, ax = plt.subplots(figsize=(6.5, 5.5), layout='constrained')
            sc = ax.pcolormesh(x, y, values.reshape(grid_shape), cmap='RdBu_r', vmin=-abs_max, vmax=abs_max)
            ax.set_aspect(aspect, adjustable='box')
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(symbol_label)
            ax.grid(True, alpha=0.25)
            fig.colorbar(sc, ax=ax, shrink=0.85)
            add_w0_secondary_axes(ax, w0)

            fig.suptitle(f"{tag.capitalize()} angular momentum density{title_suffix}, $\\omega$ index "
                         f"{i_omega} ($\\omega/\\omega_1$={omega_ratio:.4g}), {detector_type} detector, "
                         f"{detector_geometry_label}", fontsize=11)

            output_name = f"angular_momentum_{tag}_omega{i_omega}.png"
            output_img = os.path.join(png_dir, output_name)
            plt.savefig(output_img, dpi=200, bbox_inches='tight')
            print(f"Successfully saved plot to {output_img}")
            plt.close(fig)

    run_dir = os.path.dirname(radiation_filepath)
    integrated = integrate_angular_momentum_density(result, detector_type, config_path)
    append_integration_to_run_log(run_dir, detector_type, integrated,
                                  label_suffix=" -- incident beam" if is_incident else "")


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

    plot_angular_momentum(input_file)
