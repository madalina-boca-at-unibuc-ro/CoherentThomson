"""
Computes and plots the z-component of the spin (SAM), orbital (OAM), and total (TAM) angular
momentum spectral densities, their z-directed fluxes, AND the electromagnetic energy spectral
density and its z-directed flux (Poynting vector), from a run's coherently-summed Faraday tensor
(radiation_field.dat) or the incident beam's own analytic field (incident_field.dat, via
--incident) -- see theory/numerical_calculation_of_angular_momentum.md ("Applications for the
angular momentum density", sections 1-2 for AM density, 3-4 for AM flux, 5-6 for energy density/
flux) for the derivation. Formerly plot_angular_momentum.py; renamed plot_observables.py once the
energy density/flux quantities (sections 5-6 of the theory doc, added after the angular-momentum
sections were already implemented) were folded into this same script rather than split into a
second one, since they share every piece of this script's machinery (screen geometry, incident/
emitted split, run_log.txt integration/table-splicing) with nothing angular-momentum-specific in
any of it. Supersedes plot_spin_angular_momentum.py, which only computed S_z and (per that doc's
own correction) was missing the epsilon_0 = 1/(4*pi) prefactor -- see
Core::PhysUtils::AtomicUnits::epsilon_0 (phys_utils.hpp) for the same constant on the C++ side.

Angular momentum density (accumulated angular momentum per unit screen area):
$$\\frac{d\\mathcal{S}_z}{d\\omega} = \\frac{4\\epsilon_0}{\\omega}\\mathrm{Im}[\\tilde E_x^*\\tilde E_y]$$
$$\\frac{d\\mathcal{L}_z}{d\\omega} = \\frac{2\\epsilon_0}{\\omega}\\sum_i\\mathrm{Im}\\left[\\tilde E_i^*\\,\\hat L_z\\tilde E_i\\right]$$
$$\\frac{d\\mathcal{J}_z}{d\\omega} = \\frac{d\\mathcal{L}_z}{d\\omega} + \\frac{d\\mathcal{S}_z}{d\\omega}$$

Angular momentum flux (angular momentum flow per unit screen area, through the screen along its
own normal Oz -- a genuinely different physical quantity from the density above, not the same
thing in different units, see CLAUDE.md's density-vs-flux distinction from this project's earlier,
since-deleted angular-momentum scripts):
$$\\frac{d\\Sigma_{zz}}{d\\omega} = \\frac{2\\epsilon_0 c^2}{\\omega}\\mathrm{Im}\\left[\\tilde B_x^*\\tilde E_x + \\tilde B_y^*\\tilde E_y - \\tilde B_z^*\\tilde E_z\\right]$$
$$\\frac{d\\Lambda_{zz}}{d\\omega} = \\frac{2\\epsilon_0 c^2}{\\omega}\\mathrm{Im}\\left[\\tilde B_y^*\\,\\hat L_z\\tilde E_x - \\tilde B_x^*\\,\\hat L_z\\tilde E_y + \\tilde B_z^*\\tilde E_z\\right]$$
$$\\frac{d(\\Sigma_{zz}+\\Lambda_{zz})}{d\\omega} \\;\\text{(total flux)}$$

where $\\hat L_z = x\\partial_y - y\\partial_x$ (rectangular) $= \\partial_\\phi$ (circular) -- the
*same* operator for both density and flux, just combined with different field components afterward
(E with itself for density, B with E for flux); see _lz_operator.

Energy density and its z-directed flux (Poynting vector) -- theory doc sections 5-6. Unlike every
angular-momentum quantity above, these are bilinear directly in E/B (never in the vector potential
A), so their spectral-density formula ([Eq. (bilinear-spectral-density)] applied to E/B themselves)
carries **no explicit 1/omega prefactor** -- don't add one by analogy with S_z/L_z/Sigma_zz/
Lambda_zz above, that would be wrong:
$$\\frac{du}{d\\omega} = \\epsilon_0\\left(|\\tilde E_x|^2+|\\tilde E_y|^2+|\\tilde E_z|^2\\right) + \\epsilon_0 c^2\\left(|\\tilde B_x|^2+|\\tilde B_y|^2+|\\tilde B_z|^2\\right)$$
$$\\frac{dP_z}{d\\omega} = 2\\epsilon_0 c^2\\,\\mathrm{Re}\\left[\\tilde E_x\\tilde B_y^* - \\tilde E_y\\tilde B_x^*\\right]$$

**Always uses the total (physical) Faraday tensor, LR+SR+BR summed** -- there is no
long-range/short-range/boundary decomposition of any of these observables (that split only has
meaning for the derivation's intermediate F^{mu nu} terms, see
theory/FT_Faraday_tensor-direct_and_simplified_forms.md); only the sum is the actual radiated field
E/B appearing in the formulas above. BR is treated as zero if the file predates the boundary term;
for incident_field.dat, SR/BR are already identically zero so the sum reduces to LR alone.

Ex/Ey/Ez/Bx/By/Bz are recovered directly from the stored F^{mu nu} columns (E_i=c*F^{i0},
Bx=F^{32}, By=F^{13}, Bz=F^{21} -- the same convention faraday_frame_utils.py's
extract_rotated_faraday_fields documents).

**Restricted to rectangular and circular detectors** (theory doc's own "Target Screens"
instruction): a spherical detector's (theta, phi) grid has no flat local (x, y) plane for the OAM
operator's x*d/dy - y*d/dx (rectangular) or d/dphi (circular) construction to act on, so this script
rejects detector_type=spherical outright -- unlike the old, now-superseded
plot_spin_angular_momentum.py, which did support spherical for S_z alone (no spatial derivative
needed there). u/P_z need no such operator and could in principle support a spherical detector, but
this script computes every QUANTITIES entry together in one pass, so they're gated by the same
rectangular/circular check as the angular-momentum quantities rather than special-cased around it.

**OAM spatial derivatives**: rectangular uses numpy.gradient in x/y (central differences, taking the
grid spacing into account, per the theory doc's own numerical-implementation instructions); circular
uses an exact periodic centered difference in phi (numpy.gradient has no periodic-boundary mode, so
this is implemented by hand -- see _periodic_phi_derivative), which is what x*d/dy - y*d/dx reduces
to exactly at fixed r on a polar grid. Both are finite-difference approximations, not exact for a
rapidly-varying near-field phase (the same caveat this project documents elsewhere for
np.gradient-based derivatives, e.g. CLAUDE.md's Fresnel-number bullet) -- a coarse grid close to the
beam (small Fresnel number) may under-resolve the true spatial variation. Note the L_z density/flux
formulas are invariant under an isotropic rescaling of x/y (x*d/dy-y*d/dx is the rotation generator,
scale-free under x->c*x, y->c*y), so it does not matter that get_screen_coordinates returns x/y in
the config's own length unit (e.g. 'lambda') rather than atomic units -- only the screen-integrated
totals below need an explicit unit conversion, since an area DOES scale with the unit choice.

**Frequency convention**: neither .dat file exports raw atomic-unit omega, only the dimensionless
omega/omega_1 ratio (see CLAUDE.md's radiation_field.dat 'omega' column bullet) -- recovering true
fundamental_frequency in Python would mean re-parsing it out of run_log.txt's human-readable dump.
So the 1/omega prefactor in S_z/L_z/Sigma_zz/Lambda_zz above uses that same ratio, making those
plotted quantities d.../d(omega/omega_1) rather than d.../domega in raw atomic units -- consistent
with every other frequency-axis quantity this project already plots in omega/omega_1 units, and
differing from the true d.../domega only by the constant factor omega_1 (fundamental_frequency),
which does not affect the spatial pattern at any given omega index or the relative comparison
across harmonics. u/P_z have no such prefactor to rescale (see above), so this distinction doesn't
apply to them -- they're already true du/domega, dP_z/domega in atomic units, at whatever raw
frequency the omega/omega_1 index happens to correspond to.

**Screen integration**: appends a summary table (one column per quantity below) to the run's own
run_log.txt (theory doc's "Surface Integration" step), replacing any such table from a previous run
of this script on the same run_log.txt so re-running doesn't pile up duplicates. Trapezoidal
quadrature in (x, y) for rectangular, (r^2, phi) for circular (matching CircularDetector's own
equal-area radial spacing, so a constant-r^2-step trapezoidal rule is exact-in-area for any N_R) --
area weights are converted to atomic-unit length^2 (get_screen_area_weights_au), unlike the density/
flux values themselves.

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
    "applied -- only RATIOS between the columns below (S_z/L_z/J_z density, Sigma_zz/Lambda_zz/"
    "Flux_tot flux, and P_z/u) are physically meaningful; the absolute values are on an arbitrary scale and are "
    "NOT comparable to a real scattered-radiation run's own values."
)

# (result column, output-filename tag, figure title, math-symbol panel title). Angular-momentum
# density (S_z/L_z/J_z), angular-momentum flux (Sigma_zz/Lambda_zz/Flux_tot), and energy
# density/flux (u/P_z) are three deliberately different physical quantities (accumulated angular
# momentum, its flow through the screen, and the analogous energy pair), not the same thing in
# different units -- kept as one list purely because every quantity here shares the same
# plotting/integration machinery.
QUANTITIES = (
    ('dSz', 'spin', 'Spin angular momentum density', r'$d\mathcal{S}_z/d(\omega/\omega_1)$'),
    ('dLz', 'orbital', 'Orbital angular momentum density', r'$d\mathcal{L}_z/d(\omega/\omega_1)$'),
    ('dJz', 'total', 'Total angular momentum density', r'$d\mathcal{J}_z/d(\omega/\omega_1)$'),
    ('dSigmazz', 'flux_spin', 'Spin angular momentum flux', r'$d\Sigma_{zz}/d(\omega/\omega_1)$'),
    ('dLambdazz', 'flux_orbital', 'Orbital angular momentum flux', r'$d\Lambda_{zz}/d(\omega/\omega_1)$'),
    ('dFluxTotal', 'flux_total', 'Total angular momentum flux',
     r'$d(\Sigma_{zz}{+}\Lambda_{zz})/d(\omega/\omega_1)$'),
    ('du', 'energy_density', 'Electromagnetic energy density', r'$du/d(\omega/\omega_1)$'),
    ('dPz', 'energy_flux', 'Energy flux (Poynting vector)', r'$dP_z/d(\omega/\omega_1)$'),
)

# Short labels for run_log.txt's table columns (same order/columns as QUANTITIES).
_RUN_LOG_LABELS = {
    'dSz': 'S_z', 'dLz': 'L_z', 'dJz': 'J_z',
    'dSigmazz': 'Sigma_zz', 'dLambdazz': 'Lambda_zz', 'dFluxTotal': 'Flux_tot',
    'du': 'u', 'dPz': 'P_z',
}


def _complex_column(data, prefix, mu, nu):
    return (data[f'{prefix}_F{mu}{nu}_re'] + 1j * data[f'{prefix}_F{mu}{nu}_im']).to_numpy()


def _get_total_field(data, index_pairs, scale=1.0):
    """
    Generic total-(LR+SR+BR) field extractor: one complex numpy array per (mu, nu) pair in
    `index_pairs`, each scaled by `scale`. Shared by _get_E_total (E_i=c*F^{i0}, scale=C_LIGHT) and
    _get_B_total (B_i=F^{jk} directly, scale=1) -- see this module's own doc comment for why only
    the total, never a single range, is used here.
    """
    def component(mu, nu):
        value = _complex_column(data, 'LR', mu, nu) + _complex_column(data, 'SR', mu, nu)
        if {f'BR_F{mu}{nu}_re', f'BR_F{mu}{nu}_im'}.issubset(data.columns):
            value = value + _complex_column(data, 'BR', mu, nu)
        return scale * value

    return tuple(component(mu, nu) for mu, nu in index_pairs)


def _get_E_total(data):
    """Returns (Ex, Ey, Ez) from the total Faraday tensor. E_i = c*F^{i0}."""
    return _get_total_field(data, ((1, 0), (2, 0), (3, 0)), scale=C_LIGHT)


def _get_B_total(data):
    """Returns (Bx, By, Bz) from the total Faraday tensor. Bx=F^{32}, By=F^{13}, Bz=F^{21}."""
    return _get_total_field(data, ((3, 2), (1, 3), (2, 1)), scale=1.0)


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


def _lz_operator(detector_type, E, x_grid=None, y_grid=None, x_vals=None, y_vals=None):
    """
    hat{L}_z E = x*dE/dy - y*dE/dx (rectangular, via numpy.gradient) or dE/dphi (circular, via
    _periodic_phi_derivative) -- the one operator both the OAM density (Eq. primary-definition-of-Lz-
    rectangular/-circular) and OAM flux (Eq. primary-definition-of-Lambda-zz-rectangular/-circular)
    formulas need, just combined with different field components afterward by their respective
    callers. A size-1 rectangular axis has no defined derivative (a single-row/column detector has no
    spatial variation to differentiate along that axis) -- treated as zero, matching numpy.gradient's
    own minimum-2-points requirement. `x_grid`/`y_grid`/`x_vals`/`y_vals` are ignored for circular.
    """
    if detector_type == 'rectangular':
        dE_dx = np.gradient(E, x_vals, axis=0) if E.shape[0] > 1 else np.zeros_like(E)
        dE_dy = np.gradient(E, y_vals, axis=1) if E.shape[1] > 1 else np.zeros_like(E)
        return x_grid * dE_dy - y_grid * dE_dx
    return _periodic_phi_derivative(E)


def compute_observables(radiation_filepath):
    """
    Returns (result, detector_type, config_path). `result` has one row per radiation_field.dat row
    with omega/omega_1 != 0 (i_omega, omega, i_screen, dSz, dLz, dJz, dSigmazz, dLambdazz,
    dFluxTotal, du, dPz), computed from the physical total field (LR+SR+BR, see
    _get_E_total/_get_B_total).
    """
    try:
        data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    except Exception as e:
        print(f"Error reading file '{radiation_filepath}': {e}")
        sys.exit(1)

    required = {f'{p}_F{m}{n}_{part}' for p in ('LR', 'SR')
                for m, n in ((1, 0), (2, 0), (3, 0), (3, 2), (1, 3), (2, 1))
                for part in ('re', 'im')}
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
            "theory/numerical_calculation_of_angular_momentum.md's angular-momentum calculation is "
            "only implemented for those two (its own 'Target Screens' instruction).")

    x_grid = y_grid = x_vals = y_vals = None
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
        Bx, By, Bz = _get_B_total(subset)

        dSz = (4.0 * EPSILON_0 / omega_ratio) * np.imag(np.conj(Ex) * Ey)
        dSigmazz = (2.0 * EPSILON_0 * C_LIGHT ** 2 / omega_ratio) * np.imag(
            np.conj(Bx) * Ex + np.conj(By) * Ey - np.conj(Bz) * Ez)

        # No 1/omega_ratio factor here -- see this module's own doc comment: u/P_z are bilinear
        # directly in E/B, never in the vector potential A, so their spectral density carries no
        # such prefactor (unlike dSz/dSigmazz above).
        du = (EPSILON_0 * (np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2)
              + EPSILON_0 * C_LIGHT ** 2 * (np.abs(Bx) ** 2 + np.abs(By) ** 2 + np.abs(Bz) ** 2))
        dPz = 2.0 * EPSILON_0 * C_LIGHT ** 2 * np.real(Ex * np.conj(By) - Ey * np.conj(Bx))

        Ex_grid, Ey_grid, Ez_grid = (a.reshape(grid_shape) for a in (Ex, Ey, Ez))
        Bx_grid, By_grid, Bz_grid = (a.reshape(grid_shape) for a in (Bx, By, Bz))
        LEx = _lz_operator(detector_type, Ex_grid, x_grid, y_grid, x_vals, y_vals)
        LEy = _lz_operator(detector_type, Ey_grid, x_grid, y_grid, x_vals, y_vals)
        LEz = _lz_operator(detector_type, Ez_grid, x_grid, y_grid, x_vals, y_vals)

        dLz = (2.0 * EPSILON_0 / omega_ratio) * np.imag(
            np.conj(Ex_grid) * LEx + np.conj(Ey_grid) * LEy + np.conj(Ez_grid) * LEz).ravel()
        dLambdazz = (2.0 * EPSILON_0 * C_LIGHT ** 2 / omega_ratio) * np.imag(
            np.conj(By_grid) * LEx - np.conj(Bx_grid) * LEy + np.conj(Bz_grid) * Ez_grid).ravel()

        rows.append(pd.DataFrame({
            'i_omega': i_omega, 'omega': omega_ratio, 'i_screen': subset['i_screen'].to_numpy(),
            'dSz': dSz, 'dLz': dLz, 'dJz': dSz + dLz,
            'dSigmazz': dSigmazz, 'dLambdazz': dLambdazz, 'dFluxTotal': dSigmazz + dLambdazz,
            'du': du, 'dPz': dPz,
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


def integrate_observables(result, detector_type, config_path):
    """
    Screen integration (theory doc's 'Surface Integration' step): Q(omega) = sum_p Delta A_p *
    (dQ/domega)_p, for every quantity in QUANTITIES (density and flux alike -- the formula is
    identical regardless of what physical quantity is being integrated). Returns a DataFrame with
    one row per unique i_omega.
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
        row = {'i_omega': i_omega, 'omega': subset['omega'].iloc[0]}
        for column, *_ in QUANTITIES:
            row[column] = np.sum(area_weights * subset[column].to_numpy())
        rows.append(row)
    return pd.DataFrame(rows).sort_values('i_omega').reset_index(drop=True)


_RUN_LOG_SECTION_HEADER = "Screen-integrated observables (plot_observables.py)"

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
    Appends (or, on a re-run, replaces) a density/flux summary table (one column per QUANTITIES
    entry) in the run's own run_log.txt -- theory/numerical_calculation_of_angular_momentum.md's
    'Surface Integration' step, integrated over the whole screen, one row per frequency. Only this
    script's own previously-appended section (identified by _RUN_LOG_SECTION_HEADER + label_suffix,
    via _find_run_log_sections) is ever replaced -- in place, regardless of what other sections
    precede or follow it. `label_suffix` (e.g. " -- incident beam", set when run on
    incident_field.dat) keeps that table in its own section rather than overwriting the actual
    scattered-radiation run's own logged section, and also triggers _UNCALIBRATED_MAGNITUDE_CAVEAT
    being written into that section.
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
        "  Rough sum_p Delta_A_p * (dQ/domega)_p screen integration "
        "(theory/numerical_calculation_of_angular_momentum.md,",
        f"  'Surface Integration' step), {detector_type} detector, trapezoidal quadrature in "
        f"{'(x, y)' if detector_type == 'rectangular' else '(r^2, phi)'}, area weights in atomic units.",
    ]
    header_cells = ['i_omega', 'omega/omega_1'] + [f'{_RUN_LOG_LABELS[c]} [a.u.]' for c, *_ in QUANTITIES]
    lines.append(f"  {header_cells[0]:>8} {header_cells[1]:>14} "
                 + ' '.join(f'{cell:>16}' for cell in header_cells[2:]))
    for _, row in integrated.iterrows():
        cells = [f"{int(row['i_omega']):>8d}", f"{row['omega']:>14.6f}"]
        cells += [f"{row[c]:>16.6e}" for c, *_ in QUANTITIES]
        lines.append('  ' + ' '.join(cells))

    with open(log_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Appended screen-integrated observables to {log_path}")


def plot_observables(radiation_filepath):
    """
    Plots every quantity in QUANTITIES (S_z/L_z/J_z angular-momentum density,
    Sigma_zz/Lambda_zz/Flux_tot angular-momentum flux, u energy density, P_z energy flux -- one
    single-panel diverging-colormap heatmap each, since all eight are real-valued) over the
    detector screen, one figure set per configured frequency, then appends the screen-integrated
    totals to run_log.txt.

    `radiation_filepath` may point at either a run's radiation_field.dat or its incident_field.dat
    sibling, detected by filename -- see this module's own doc comment for how each is handled.
    """
    is_incident = os.path.basename(radiation_filepath) == "incident_field.dat"
    title_suffix = " (incident beam)" if is_incident else ""
    if is_incident:
        print(_UNCALIBRATED_MAGNITUDE_CAVEAT)

    result, detector_type, config_path = compute_observables(radiation_filepath)

    x_centers, y_centers, x_label, y_label = get_screen_coordinates(radiation_filepath, config_path)
    detector_geometry_label = get_detector_geometry_label(detector_type, config_path)

    aspect = 'equal'
    if detector_type == 'rectangular':
        Nx = int(read_config_value('rectangular_detector_Nx', config_path)[0])
        Ny = int(read_config_value('rectangular_detector_Ny', config_path)[0])
        grid_shape = (Nx, Ny)
        # Flat shading needs cell *corners* (one more point per axis than the data) -- coarse
        # rectangular grids already render as reasonably-shaped axis-aligned squares this way, so
        # there's no strong reason to switch this one to gouraud too (see the circular branch below).
        x, y = get_rectangular_cell_edges(config_path)
        shading = 'flat'
        axes_unit = read_config_value('rectangular_detector_x_min', config_path)[1]
    else:
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        grid_shape = (N_R, N_phi)
        # Flat-shaded quads over a coarse (r, phi) grid mapped to Cartesian x/y render as visible
        # pie-slice facets, worst at large radius where each azimuthal wedge spans a wide arc length
        # (confirmed on a real 50-electron/64-phi-point run). Gouraud shading fixes this by linearly
        # interpolating color across each cell from its four corner *values* instead of flat-filling
        # it -- needs cell *centers* (matplotlib requires X/Y/C all the same shape for shading=
        # 'gouraud', unlike flat shading's one-bigger edges array) and the values at those centers,
        # both already available from get_screen_coordinates/compute_observables; no extra data
        # or a different (edges) array needed, unlike the flat-shaded branch above. Preferred over a
        # Delaunay/tricontourf re-triangulation of the point cloud (an alternative considered and
        # rejected): this still uses the detector's own known grid connectivity, so it can't
        # fabricate data across the R_min>0 center hole the way a blind triangulation would, and
        # isn't at risk of a degenerate-triangulation crash at R_min=0, where the whole inner ring
        # collapses to N_phi coincident points at the origin -- gouraud just draws a degenerate
        # (zero-area) quad there, same as it always has for any repeated vertex. Every quantity here
        # (unlike plot_field.py's Phase panel) is a plain real-valued scalar with no branch-cut wrap,
        # so gouraud is safe for all of them, not just a subset.
        x = x_centers.reshape(grid_shape)
        y = y_centers.reshape(grid_shape)
        shading = 'gouraud'
        axes_unit = read_config_value('circular_detector_R_min', config_path)[1]

    w0 = get_laser_lg_w0_in_axes_units(radiation_filepath, axes_unit) if axes_unit else None

    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder", MODULE_NAME,
                           "incident" if is_incident else "emitted")
    os.makedirs(png_dir, exist_ok=True)

    for i_omega, subset in result.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_ratio = subset['omega'].iloc[0]

        for column, tag, name, symbol_label in QUANTITIES:
            values = subset[column].to_numpy()
            abs_max = np.abs(values).max() if values.size else None

            fig, ax = plt.subplots(figsize=(6.5, 5.5), layout='constrained')
            sc = ax.pcolormesh(x, y, values.reshape(grid_shape), cmap='RdBu_r', vmin=-abs_max, vmax=abs_max,
                              shading=shading)
            ax.set_aspect(aspect, adjustable='box')
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(symbol_label)
            ax.grid(True, alpha=0.25)
            fig.colorbar(sc, ax=ax, shrink=0.85)
            add_w0_secondary_axes(ax, w0)

            fig.suptitle(f"{name}{title_suffix}, $\\omega$ index {i_omega} "
                         f"($\\omega/\\omega_1$={omega_ratio:.4g}), {detector_type} detector, "
                         f"{detector_geometry_label}", fontsize=11)

            output_name = f"observable_{tag}_omega{i_omega}.png"
            output_img = os.path.join(png_dir, output_name)
            plt.savefig(output_img, dpi=200, bbox_inches='tight')
            print(f"Successfully saved plot to {output_img}")
            plt.close(fig)

    run_dir = os.path.dirname(radiation_filepath)
    integrated = integrate_observables(result, detector_type, config_path)
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

    plot_observables(input_file)
