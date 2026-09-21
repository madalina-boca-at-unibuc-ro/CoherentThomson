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

**Frequency convention (FIXED -- see CLAUDE.md's "missing 1/omega" bug note)**: neither .dat file
exports raw atomic-unit omega, only the dimensionless omega/omega_1 ratio (see CLAUDE.md's
radiation_field.dat 'omega' column bullet). An earlier version of this script used that ratio
directly as the divisor in S_z/L_z/Sigma_zz/Lambda_zz's 1/omega prefactor, making those four (and
J_z/Flux_tot, built from them) d.../d(omega/omega_1) = omega_1 * d.../domega_true -- while u/P_z
(no 1/omega prefactor at all, per the theory doc) stayed true d.../domega_true, unscaled. That put
the two families on different frequency bases, off by the constant factor omega_1
(fundamental_frequency): comparing e.g. J_z/u directly gave a bare dimensionless number
(coincidentally close to the topological charge m at the fundamental, since omega_1/omega=1 there)
instead of a quantity with units of inverse time. Fixed by reading the run's own raw
fundamental_frequency (atomic units, inverse time) back out of run_log.txt
(_read_fundamental_frequency_au -- Logging::write_run_log already prints
`sim_par.fundamental_frequency * c`, i.e. the internal omega/c-convention value converted back to a
genuine angular frequency) and multiplying it onto the exported omega/omega_1 ratio to recover the
actual raw omega (atomic units) used as the 1/omega divisor everywhere below -- so all eight
QUANTITIES are now genuine d.../domega on the same basis, and e.g. int_dJz_domega/int_du_domega now
comes out numerically close to m/omega_1 (inverse time), not the bare integer m.

**Screen integration**: appends a summary table (one column per quantity below) to the run's own
run_log.txt (theory doc's "Surface Integration" step), replacing any such table from a previous run
of this script on the same run_log.txt so re-running doesn't pile up duplicates. Trapezoidal
quadrature in (x, y) for rectangular, (r^2, phi) for circular (matching CircularDetector's own
equal-area radial spacing, so a constant-r^2-step trapezoidal rule is exact-in-area for any N_R) --
area weights are converted to atomic-unit length^2 (get_screen_area_weights_au), unlike the density/
flux values themselves.

**Flux/density ratios**: alongside that table, also prints (stdout) and appends to run_log.txt a
second table of each screen-integrated flux total divided by its corresponding density total
(Sigma_zz/S_z, Lambda_zz/L_z, Flux_tot/J_z, P_z/u -- see _RATIO_PAIRS), each expected to land close to
C_LIGHT: flux equals density times transport speed, for radiation propagating at c (same relationship
as the Poynting vector to the EM energy density, S=c*u) -- see CLAUDE.md's own confirmation of this on
a real run (~3e-5 relative agreement for all four pairs). A pair whose density integrates close to
zero (e.g. S_z can nearly cancel for some polarization/mode combinations) gives a noise-dominated,
not-necessarily-close-to-C_LIGHT ratio for that one pair -- expected, not itself a sign of a bug.

**OAM-to-energy ratio**: a third table (_ENERGY_NORMALIZED_PAIRS, Lambda_zz/P_z -- the orbital
angular-momentum *flux* over the energy *flux*/Poynting vector, not the density pair L_z/u this used
before) divides the screen-integrated orbital-AM-flux total by the screen-integrated energy-flux
total, expected to land close to m/omega_N (inverse time, atomic units; m = the beam's own
topological charge, i.e. config.cfg's laser_lg_l; omega_N = that row's actual raw angular frequency,
(omega/omega_1)*fundamental_frequency_au, fundamental_frequency_au read back out of run_log.txt by
_read_fundamental_frequency_au) -- the angular-momentum/energy = m/omega relation for radiation
carrying m units of angular momentum per photon of energy omega (hbar=1). Since flux = c * density
for both the numerator and denominator (see _RATIO_PAIRS' own C_LIGHT check), Lambda_zz/P_z and
L_z/u are algebraically identical (the c cancels) -- the flux pair is used here because Lambda_zz/P_z
draws on the B-field bilinears (never used by L_z/u, which is pure-E), so agreement between this
ratio and the theoretical m/omega_N is an independent check on the flux formulas too, not just a
restatement of the density one. The theoretical m/omega_N value is appended alongside the numerical
ratio (not just quoted once in the caption) so every row can be checked directly, since omega_N
scales with harmonic index N and a single caption value would only be valid for N=1 (see CLAUDE.md's
"missing 1/omega" bug-fix note for why this needs the run's actual raw omega_1, not the dimensionless
omega/omega_1 ratio the .dat files export).

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
    ('dSz', 'spin', 'Spin angular momentum density', r'$d\mathcal{S}_z/d\omega$'),
    ('dLz', 'orbital', 'Orbital angular momentum density', r'$d\mathcal{L}_z/d\omega$'),
    ('dJz', 'total', 'Total angular momentum density', r'$d\mathcal{J}_z/d\omega$'),
    ('dSigmazz', 'flux_spin', 'Spin angular momentum flux', r'$d\Sigma_{zz}/d\omega$'),
    ('dLambdazz', 'flux_orbital', 'Orbital angular momentum flux', r'$d\Lambda_{zz}/d\omega$'),
    ('dFluxTotal', 'flux_total', 'Total angular momentum flux',
     r'$d(\Sigma_{zz}{+}\Lambda_{zz})/d\omega$'),
    ('du', 'energy_density', 'Electromagnetic energy density', r'$du/d\omega$'),
    ('dPz', 'energy_flux', 'Energy flux (Poynting vector)', r'$dP_z/d\omega$'),
)

# Short labels for run_log.txt's table columns (same order/columns as QUANTITIES).
_RUN_LOG_LABELS = {
    'dSz': 'S_z', 'dLz': 'L_z', 'dJz': 'J_z',
    'dSigmazz': 'Sigma_zz', 'dLambdazz': 'Lambda_zz', 'dFluxTotal': 'Flux_tot',
    'du': 'u', 'dPz': 'P_z',
}

# (flux column, density column, ratio label). Each pair's flux/density ratio is expected to equal
# C_LIGHT: a flux (crossing the screen per unit time/area) equals its density (sitting on/near the
# screen per unit area) times the transport speed, same relationship as the Poynting vector to the EM
# energy density (S=c*u) -- see CLAUDE.md's "Flux_total_quantity / Density_total_quantity == C_LIGHT"
# verification note, which found this to ~3e-5 relative on a real run for all four pairs below.
_RATIO_PAIRS = (
    ('dSigmazz', 'dSz', 'Sigma_zz/S_z'),
    ('dLambdazz', 'dLz', 'Lambda_zz/L_z'),
    ('dFluxTotal', 'dJz', 'Flux_tot/J_z'),
    ('dPz', 'du', 'P_z/u'),
)

# (numerator column, denominator column, ratio label). Lambda_zz/P_z is the OAM-to-energy ratio: for
# radiation carrying m units of angular momentum per photon of energy omega (hbar=1, atomic units),
# int_dLambdazz_domega/int_dPz_domega is expected to land close to m/omega_N (inverse time, atomic
# units, omega_N the row's actual raw angular frequency -- see this module's own "OAM-to-energy
# ratio" doc comment for why the flux pair is used instead of the algebraically-identical density
# pair L_z/u). See CLAUDE.md's "missing 1/omega" bug-fix note, whose own int_dJz_domega/int_du_domega
# check (the sum of S_z/u and L_z/u) motivated this fix in the first place. Kept in its own
# tuple/table rather than folded into _RATIO_PAIRS above: that table's caption specifically expects
# C_LIGHT, which does not apply here.
_ENERGY_NORMALIZED_PAIRS = (
    ('dLambdazz', 'dPz', 'Lambda_zz/P_z'),
)


# Matches Logging::write_run_log's own write_kv_au format ("  <label, left-padded to 28> <value, "
# scientific> a.u.\n") for the specific "fundamental_frequency" line.
_FUNDAMENTAL_FREQUENCY_RE = re.compile(r'(?m)^\s*fundamental_frequency\s+([+-]?[\d.]+[eE][+-]?\d+)\s+a\.u\.')


def _read_fundamental_frequency_au(run_dir):
    """
    Returns the run's own raw fundamental_frequency (omega_1, atomic units, inverse time), parsed
    out of run_log.txt -- Logging::write_run_log prints `sim_par.fundamental_frequency * c`, i.e.
    the internal omega/c-convention value already converted back to a genuine angular frequency
    (run_log.cpp), not the omega/c value itself. Needed to convert radiation_field.dat/
    incident_field.dat's exported dimensionless omega/omega_1 ratio back into a real omega for the
    1/omega prefactor in dSz/dLz/dSigmazz/dLambdazz below -- see this module's own "Frequency
    convention" doc comment for the bug this fixes.
    """
    log_path = os.path.join(run_dir, "run_log.txt")
    if not os.path.exists(log_path):
        raise FileNotFoundError(
            f"'{log_path}' not found -- this run predates run_log.txt (Logging::write_run_log); "
            "re-run the solver to regenerate it.")
    with open(log_path) as f:
        content = f.read()
    match = _FUNDAMENTAL_FREQUENCY_RE.search(content)
    if not match:
        raise ValueError(f"Could not find a 'fundamental_frequency' line in '{log_path}'.")
    return float(match.group(1))


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

    fundamental_frequency_au = _read_fundamental_frequency_au(os.path.dirname(radiation_filepath))

    detector_type, _ = read_config_value('detector_type', config_path)
    if detector_type not in ('rectangular', 'circular'):
        raise ValueError(
            f"detector_type '{detector_type}' is not a flat rectangular/circular screen -- "
            "theory/numerical_calculation_of_angular_momentum.md's angular-momentum calculation is "
            "only implemented for those two (its own 'Target Screens' instruction).")

    # The L_z operator (x*d/dy - y*d/dx / d/dphi, see _lz_operator) and every observable this script
    # computes implicitly assume the screen's own local x/y axes coincide with the simulation frame's
    # x/y -- i.e. the screen normal is along the laser's Oz axis, forward or backward. The laser
    # itself is now always fixed along Oz (see CLAUDE.md), but detector_direction_theta/phi can still
    # point the screen anywhere -- so that assumption isn't automatically satisfied and has to be
    # checked explicitly, rather than silently computing a wrong L_z for a tilted screen. Checked via
    # the direction vector's own x/y components (Core::MathUtils::create_unit_light_like_vector's
    # convention: dir = (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))) rather than theta alone,
    # so it's robust to any equivalent (theta, phi) parametrization of +-Oz.
    theta_s, theta_unit = read_config_value('detector_direction_theta', config_path)
    phi_s, phi_unit = read_config_value('detector_direction_phi', config_path)
    theta = float(theta_s) * convert_unit_to_number(theta_unit, config_path)
    phi = float(phi_s) * convert_unit_to_number(phi_unit, config_path)
    dir_x = np.sin(theta) * np.cos(phi)
    dir_y = np.sin(theta) * np.sin(phi)
    transverse_tol = 1e-9
    if abs(dir_x) > transverse_tol or abs(dir_y) > transverse_tol:
        raise ValueError(
            f"detector_direction_theta/phi = ({theta_s} {theta_unit}, {phi_s} {phi_unit}) does not point "
            "along (0,0,1) or (0,0,-1) in the simulation frame -- every observable this script computes "
            "(S_z/L_z/J_z, their fluxes, u/P_z) assumes the screen is transverse to the laser's own "
            "propagation axis (Oz), which a tilted detector_direction violates. Not yet supported here; "
            "use detector_direction_theta = 0.0 pi or 1.0 pi (any detector_direction_phi).")

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
        omega_au = omega_ratio * fundamental_frequency_au

        Ex, Ey, Ez = _get_E_total(subset)
        Bx, By, Bz = _get_B_total(subset)

        dSz = (4.0 * EPSILON_0 / omega_au) * np.imag(np.conj(Ex) * Ey)
        dSigmazz = (2.0 * EPSILON_0 * C_LIGHT ** 2 / omega_au) * np.imag(
            np.conj(Bx) * Ex + np.conj(By) * Ey - np.conj(Bz) * Ez)

        # No 1/omega factor here -- see this module's own doc comment: u/P_z are bilinear
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

        dLz = (2.0 * EPSILON_0 / omega_au) * np.imag(
            np.conj(Ex_grid) * LEx + np.conj(Ey_grid) * LEy + np.conj(Ez_grid) * LEz).ravel()
        dLambdazz = (2.0 * EPSILON_0 * C_LIGHT ** 2 / omega_au) * np.imag(
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


def compute_ratio_table(integrated, pairs):
    """
    Returns a DataFrame (one row per i_omega) of numerator/denominator ratios for each
    (numerator column, denominator column, label) entry in `pairs` -- called with _RATIO_PAIRS
    (flux/density, expected ~ C_LIGHT) and _ENERGY_NORMALIZED_PAIRS (OAM/energy, expected ~
    m/omega_N). A near-zero denominator (e.g. S_z for a near-cancelling spin term at high
    topological charge) gives a noise-dominated ratio for that one pair -- expected, not a bug.
    """
    rows = []
    for _, row in integrated.iterrows():
        entry = {'i_omega': row['i_omega'], 'omega': row['omega']}
        for num_col, den_col, label in pairs:
            with np.errstate(divide='ignore', invalid='ignore'):
                entry[label] = row[num_col] / row[den_col]
        rows.append(entry)
    return pd.DataFrame(rows)


def add_theoretical_oam_energy_ratio(ratios, fundamental_frequency_au, m):
    """
    Adds an 'm/omega_N (theory)' column to an _ENERGY_NORMALIZED_PAIRS ratio table: the theoretical
    OAM/energy ratio m/omega_N for that row's own actual raw angular frequency omega_N =
    (omega/omega_1)*fundamental_frequency_au -- computed per row (not once from omega_1 alone) since
    omega_N scales with harmonic index N, so the theoretical target is only m/omega_1 at the
    fundamental itself. Meant to sit next to the numerical Lambda_zz/P_z column so the two can be
    compared directly row by row, rather than only quoted once in a caption.
    """
    ratios = ratios.copy()
    omega_N_au = ratios['omega'].to_numpy() * fundamental_frequency_au
    with np.errstate(divide='ignore', invalid='ignore'):
        ratios['m/omega_N (theory)'] = m / omega_N_au
    return ratios


_RUN_LOG_SECTION_HEADER = "Screen-integrated observables (plot_observables.py)"
_RUN_LOG_RATIO_HEADER = "Flux/density ratios (plot_observables.py)"
_RUN_LOG_ENERGY_RATIO_HEADER = "OAM-to-energy ratios (plot_observables.py)"

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


def _write_run_log_section(run_dir, header, body_lines):
    """
    Splices `body_lines` into run_log.txt under a dash-underlined `header` (matching
    Logging::write_run_log's own section convention), replacing this same header's previous section
    in place if a prior run of this script already wrote one (via _find_run_log_sections), regardless
    of what other sections precede or follow it -- so re-running doesn't pile up duplicates and can't
    disturb an unrelated section (e.g. a differently-`label_suffix`d one, see _find_run_log_sections'
    own doc comment for the bug this guards against). Shared by append_integration_to_run_log and
    append_ratios_to_run_log, which previously duplicated this splice logic.
    """
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

    lines = [content.rstrip("\n"), "", header, "-" * len(header)] + body_lines
    with open(log_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Appended '{header}' section to {log_path}")


def append_integration_to_run_log(run_dir, detector_type, integrated, label_suffix=""):
    """
    Appends (or, on a re-run, replaces) a density/flux summary table (one column per QUANTITIES
    entry) in the run's own run_log.txt -- theory/numerical_calculation_of_angular_momentum.md's
    'Surface Integration' step, integrated over the whole screen, one row per frequency.
    `label_suffix` (e.g. " -- incident beam", set when run on incident_field.dat) keeps that table in
    its own section rather than overwriting the actual scattered-radiation run's own logged section,
    and also triggers _UNCALIBRATED_MAGNITUDE_CAVEAT being written into that section.
    """
    body = []
    if label_suffix:
        body += [f"  {_UNCALIBRATED_MAGNITUDE_CAVEAT}", ""]
    body += [
        "  Rough sum_p Delta_A_p * (dQ/domega)_p screen integration "
        "(theory/numerical_calculation_of_angular_momentum.md,",
        f"  'Surface Integration' step), {detector_type} detector, trapezoidal quadrature in "
        f"{'(x, y)' if detector_type == 'rectangular' else '(r^2, phi)'}, area weights in atomic units.",
    ]
    header_cells = ['i_omega', 'omega/omega_1'] + [f'{_RUN_LOG_LABELS[c]} [a.u.]' for c, *_ in QUANTITIES]
    body.append(f"  {header_cells[0]:>8} {header_cells[1]:>14} "
                + ' '.join(f'{cell:>16}' for cell in header_cells[2:]))
    for _, row in integrated.iterrows():
        cells = [f"{int(row['i_omega']):>8d}", f"{row['omega']:>14.6f}"]
        cells += [f"{row[c]:>16.6e}" for c, *_ in QUANTITIES]
        body.append('  ' + ' '.join(cells))

    _write_run_log_section(run_dir, _RUN_LOG_SECTION_HEADER + label_suffix, body)


def _append_ratio_table_to_run_log(run_dir, ratios, pairs, header, caption_lines, label_suffix="",
                                   extra_value_columns=()):
    """
    Appends (or, on a re-run, replaces) a numerator/denominator ratio table (`pairs`, as computed by
    compute_ratio_table) in the run's own run_log.txt, one row per frequency. Shared by
    append_ratios_to_run_log (_RATIO_PAIRS, expected ~ C_LIGHT) and append_energy_ratios_to_run_log
    (_ENERGY_NORMALIZED_PAIRS, expected ~ m/omega_N). `extra_value_columns` (column name already
    present in `ratios`, e.g. from add_theoretical_oam_energy_ratio) are printed as extra columns
    as-is, after the `pairs` ratio columns -- not computed as a ratio of two other columns themselves.
    """
    body = list(caption_lines)
    ratio_labels = [label for _, _, label in pairs]
    header_cells = ['i_omega', 'omega/omega_1'] + ratio_labels + list(extra_value_columns)
    body.append(f"  {header_cells[0]:>8} {header_cells[1]:>14} "
                + ' '.join(f'{cell:>16}' for cell in header_cells[2:]))
    for _, row in ratios.iterrows():
        cells = [f"{int(row['i_omega']):>8d}", f"{row['omega']:>14.6f}"]
        cells += [f"{row[label]:>16.6e}" for _, _, label in pairs]
        cells += [f"{row[col]:>16.6e}" for col in extra_value_columns]
        body.append('  ' + ' '.join(cells))

    _write_run_log_section(run_dir, header + label_suffix, body)


def append_ratios_to_run_log(run_dir, ratios, label_suffix=""):
    """
    Appends (or, on a re-run, replaces) the flux/density ratio table (_RATIO_PAIRS,
    compute_ratio_table) in the run's own run_log.txt, one row per frequency -- each column
    expected to land close to C_LIGHT, the sanity check documented in CLAUDE.md's
    "Flux_total_quantity / Density_total_quantity == C_LIGHT" note.
    """
    caption_lines = [
        f"  Each column is flux_total/density_total for one quantity pair, expected ~ C_LIGHT "
        f"({C_LIGHT:.6f} a.u.) since flux = density * transport speed c, for radiation propagating "
        "at c (same relationship as the Poynting vector to the EM energy density, S=c*u).",
    ]
    _append_ratio_table_to_run_log(run_dir, ratios, _RATIO_PAIRS, _RUN_LOG_RATIO_HEADER,
                                    caption_lines, label_suffix=label_suffix)


def append_energy_ratios_to_run_log(run_dir, ratios, m, label_suffix=""):
    """
    Appends (or, on a re-run, replaces) the OAM-to-energy ratio table (_ENERGY_NORMALIZED_PAIRS,
    compute_ratio_table, plus the theoretical 'm/omega_N (theory)' column from
    add_theoretical_oam_energy_ratio) in the run's own run_log.txt, one row per frequency -- the
    numerical Lambda_zz/P_z column is expected to match the theoretical column, both ~ m/omega_N
    (inverse time, atomic units) for radiation carrying m units of angular momentum per photon of
    energy omega (hbar=1 units); see CLAUDE.md's "missing 1/omega" bug-fix note.
    """
    caption_lines = [
        f"  m (topological charge, config.cfg's laser_lg_l) = {m}. Each row's Lambda_zz/P_z is the "
        "angular-momentum flux integrated over the screen divided by the energy flux integrated over "
        "the screen; 'm/omega_N (theory)' is the theoretical prediction m/omega_N for that same row's "
        "own actual raw angular frequency omega_N (not just omega_1) -- the two columns are expected "
        "to agree, for radiation carrying m units of angular momentum per photon of energy omega "
        "(hbar=1).",
    ]
    _append_ratio_table_to_run_log(run_dir, ratios, _ENERGY_NORMALIZED_PAIRS, _RUN_LOG_ENERGY_RATIO_HEADER,
                                    caption_lines, label_suffix=label_suffix,
                                    extra_value_columns=('m/omega_N (theory)',))


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
    label_suffix = " -- incident beam" if is_incident else ""
    integrated = integrate_observables(result, detector_type, config_path)
    append_integration_to_run_log(run_dir, detector_type, integrated, label_suffix=label_suffix)

    ratios = compute_ratio_table(integrated, _RATIO_PAIRS)
    print(f"Flux/density ratios (expect ~ C_LIGHT = {C_LIGHT:.6f} a.u.):")
    print(ratios.to_string(index=False))
    append_ratios_to_run_log(run_dir, ratios, label_suffix=label_suffix)

    fundamental_frequency_au = _read_fundamental_frequency_au(run_dir)
    m = int(read_config_value('laser_lg_l', config_path)[0])
    energy_ratios = compute_ratio_table(integrated, _ENERGY_NORMALIZED_PAIRS)
    energy_ratios = add_theoretical_oam_energy_ratio(energy_ratios, fundamental_frequency_au, m)
    print(f"OAM-to-energy ratios (m={m}; numerical Lambda_zz/P_z column expected to match the "
          "theoretical m/omega_N column):")
    print(energy_ratios.to_string(index=False))
    append_energy_ratios_to_run_log(run_dir, energy_ratios, m, label_suffix=label_suffix)


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
