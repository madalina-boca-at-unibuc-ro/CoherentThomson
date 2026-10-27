"""
Laguerre-Gauss decomposition of one run's emitted field on its own detector screen (Stage 2 of
theory/LG_mode_decomposition.md; spec theory/lg_screen_decomposition_simple.md; implementation notes and
conventions theory/lg_decomposition_implementation.md).

For every frequency, each selected field component E (default: both lab-frame helicity components
E_+ = (Ex - i Ey)/sqrt(2), E_- = (Ex + i Ey)/sqrt(2)) is fitted, for each trial waist w0, by weighted complex
least squares

    min_c  sum_i A_i |E_i - sum_{p,l} c_{p,l} u_{p,l}(r_i, phi_i, s; w0)|^2

with the unit-norm LG modes u_{p,l} (waist at the origin, s = the screen distance along its normal) and the
screen's trapezoidal area weights A_i. Reported per (frequency, component, w0): the complex coefficients,
eps = sqrt(sum A|E - E_rec|^2 / sum A|E|^2), captured = sum|c|^2 / sum A|E|^2 (= 1 - eps^2 when the modes are
orthonormal on the screen) and the condition number of the weighted design matrix.

The azimuthal index is fixed from angular-momentum conservation about lab +z: J = N(m + sigma) with m, sigma the
laser's l and spin about lab +z (both times sign(n_laser,z)), and the helicity-s component (s = +-1 about lab +z)
carries l = J - s. Everything (grid azimuth, field components, l) is in the lab frame, so a screen facing -Oz needs
no special case beyond mapping its local grid to lab coordinates.

For electrons at rest and an LG laser, each helicity component is also compared with a parameter-free
coherent-source prediction (coherent_source_prediction: the laser mode to the N-th power, Fresnel-propagated, times
the single-electron angular factor theta^|dl|), reported as a fidelity. Harmonics N >= 2 of a structured laser are
not single LG modes, so this, not the dominant-mode w0, is their test.

    python3 py_scripts/radiation/lg_decomposition.py [--incident] [--harmonic N] [--component plus minus Ex Ey]
        [--l L] [--p-min P] [--p-max P] [--w0 V1,V2,... [--w0-unit U]] [--n-w0 K] [run_folder]
    python3 py_scripts/radiation/lg_decomposition.py --self-test

--incident fits incident_field.dat (the laser's own mode in its own frame, on the plane s = the detector distance
downstream of the waist, incident_config.cfg; s = 0, the waist, for runs before 2026-10-27): a pure LG_p^l laser
must come out as the single coefficient (p, l = laser_lg_l) in E_sigma with eps ~ 0, and its series waist must be
the laser's w0 at z0 = 0.
"""
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd
import scipy.special
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file, field_config_path, is_incident_screen_config
from plot_field import read_config_value, read_length, get_screen_coordinates, DETECTOR_AXES_UNIT
from faraday_frame_utils import BOHR_RADIUS_M, C_LIGHT, convert_unit_to_number, get_detector_local_rotation
from plot_observables import (
    _get_E_total,
    _read_fundamental_frequency_au,
    _read_laser_direction,
    _read_laser_spin,
    _write_run_log_section,
    get_screen_area_weights_au,
)

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

_RUN_LOG_HEADER = "LG decomposition (lg_decomposition.py)"
_COMPONENTS = ('plus', 'minus', 'Ex', 'Ey')
# Default trial-waist scan: _N_W0_DEFAULT points log-spaced over [1/_W0_SCAN_FACTOR, _W0_SCAN_FACTOR] x the
# second-moment estimate (_estimate_w0).
_W0_SCAN_FACTOR = 2.0
_N_W0_DEFAULT = 41
# Optional waist-position scan (--scan-z0): _N_Z0_DEFAULT points over +-_Z0_SCAN_ZR Rayleigh ranges of the w0
# scan's central waist.
_N_Z0_DEFAULT = 41
_Z0_SCAN_ZR = 3.0
# Weighted design matrix condition number above which the coefficients are flagged as not interpretable.
_COND_WARN = 1e8
# Components with a smaller share of the transverse power (E_+ vs E_-) get no coefficient figure: their fit is noise.
_MIN_HELICITY_FRACTION_PLOT = 1e-4
# Tolerance for treating the laser spin as exactly 0 or +-1 (otherwise J is not an integer).
_SPIN_TOL = 1e-9
_TRANSVERSE_TOL = 1e-9
# Grid step relative to the beam width w(s) at the reported w0 above which the fit is flagged as under-resolved:
# circular grids compare the r^2 step with w^2 (the trapezoid's error is ~(2h/w^2)^2/12), rectangular ones the
# larger of dx, dy with w.
_MAX_STEP_OVER_WIDTH = {'circular': 0.25, 'rectangular': 0.5}
# A fit is interpretable (modes orthonormal on the screen) when |captured - (1 - eps^2)| is below this.
_CAPTURED_TOL = 0.05
# Below this |t| (= |c_1/c_0|/sqrt(1+|l|)) the summary doesn't compare t with c_2/c_1: both are residual-sized.
_MIN_T_FOR_CHECK = 0.05
# Coherent-source prediction: source radius cap (laser waists) without a beam radius, radial grid (at least
# _SOURCE_MIN_POINTS, else _SOURCE_POINTS_PER_CYCLE per cycle of the Bessel/Fresnel phase), screen radii per
# Bessel chunk, and the prediction's power fraction whose radius bounds the "outside" share.
_SOURCE_RMAX_W0 = 8.0
_SOURCE_MIN_POINTS = 4000
_SOURCE_POINTS_PER_CYCLE = 40
_FRESNEL_CHUNK = 256
_PREDICTION_CONTAINED = 0.999


# --------------------------------------------------------------------------------------------------------------
# LG basis
# --------------------------------------------------------------------------------------------------------------

def generalized_laguerre(p, alpha, x):
    """L_p^alpha(x) by the three-term recurrence (exact for integer p >= 0)."""
    L_prev = np.ones_like(x)
    if p == 0:
        return L_prev
    L = 1.0 + alpha - x
    for k in range(1, p):
        L_prev, L = L, ((2 * k + 1 + alpha - x) * L - (k + alpha) * L_prev) / (k + 1)
    return L


def lg_mode(p, l, r, phi, s, w0, k):
    """
    Unit-norm LG mode (theory/lg_screen_decomposition_simple.md), exp(-i omega t) convention, waist at s = 0:

        u = sqrt(2 p!/(pi (p+|l|)!)) / w (sqrt2 r/w)^|l| L_p^|l|(2r^2/w^2) e^{-r^2/w^2}
            e^{i l phi} e^{i k r^2/(2 R_c)} e^{-i (2p+|l|+1) psi}

    so that int |u|^2 dA = 1 on the infinite plane at any s. The carrier e^{iks} is left out: it is the same at
    every screen point (s is the screen distance), so it only shifts the phase of every coefficient by k s.
    """
    a = abs(l)
    z_R = 0.5 * k * w0 ** 2
    w = w0 * math.sqrt(1.0 + (s / z_R) ** 2)
    inv_Rc = s / (s ** 2 + z_R ** 2)
    gouy = math.atan2(s, z_R)
    norm = math.exp(0.5 * (math.log(2.0) + math.lgamma(p + 1) - math.log(math.pi) - math.lgamma(p + a + 1)))
    rho_sq = (r / w) ** 2
    amplitude = (norm / w) * (2.0 * rho_sq) ** (0.5 * a) * generalized_laguerre(p, a, 2.0 * rho_sq) * np.exp(-rho_sq)
    phase = l * phi + 0.5 * k * r ** 2 * inv_Rc - (2 * p + a + 1) * gouy
    return amplitude * np.exp(1j * phase)


def build_basis(modes, r, phi, s, w0, k):
    """Design matrix B (n_points, n_modes): one lg_mode column per (p, l) in `modes`."""
    return np.stack([lg_mode(p, l, r, phi, s, w0, k) for p, l in modes], axis=1)


def fit_modes(target, basis, weights):
    """
    Weighted complex least squares (sqrt(A_i)-scaled rows, numpy's SVD-based lstsq; the Gram matrix is never
    formed). Returns a dict with the coefficients, eps, captured, cond, rank and a status string.
    """
    norm_sq = np.sum(weights * np.abs(target) ** 2)
    result = {'coeffs': np.full(basis.shape[1], np.nan + 0j), 'eps': np.nan, 'captured': np.nan,
              'cond': np.nan, 'rank': 0, 'status': 'ok'}
    if not (np.isfinite(norm_sq) and norm_sq > 0.0):
        result['status'] = 'zero-norm component: eps undefined'
        return result

    sqrt_w = np.sqrt(weights)
    coeffs, _, rank, singular_values = np.linalg.lstsq(sqrt_w[:, None] * basis, sqrt_w * target, rcond=None)
    residual = target - basis @ coeffs
    result.update(coeffs=coeffs, rank=int(rank),
                  eps=math.sqrt(np.sum(weights * np.abs(residual) ** 2) / norm_sq),
                  captured=float(np.sum(np.abs(coeffs) ** 2) / norm_sq),
                  cond=float(singular_values[0] / singular_values[-1]) if singular_values[-1] > 0 else np.inf)
    if rank < basis.shape[1]:
        result['status'] = f'rank deficient ({rank}/{basis.shape[1]})'
    elif result['cond'] > _COND_WARN:
        result['status'] = f'ill-conditioned (cond > {_COND_WARN:.0e})'
    return result


def series_waist(coeffs, modes, w_b, k):
    """
    Waist and waist position implied by the first coefficient ratio, assuming the target is one LG_0^l beam with
    its own complex beam parameter (theory/lg_decomposition_implementation.md, "Series waist"). On the basis waist
    plane (s = 0) such a beam is r^|l| e^{-alpha r^2} e^{il phi}; the Laguerre generating function gives
    c_{p+1}/c_p = t sqrt((p+1+|l|)/(p+1)) with t = (alpha w_b^2 - 1)/(alpha w_b^2 + 1). Inverting t from c_1/c_0,
    q = -i z_Rb (1 - t)/(1 + t) (z_Rb = k w_b^2/2) is the beam parameter at s = 0, so the waist sits at
    z0 = -Re q along the propagation direction and w0 = sqrt(2 (-Im q)/k).

    Returns None unless `modes` has a single l and contains p = 0 and 1; otherwise a dict with t (from c_1/c_0),
    t2 (the same quantity from c_2/c_1, NaN without p = 2: equal to t for a pure LG_0^l beam), w0 and z0.
    """
    if len({l for _, l in modes}) != 1:
        return None
    a = abs(modes[0][1])
    c = {p: coeff for (p, _), coeff in zip(modes, coeffs)}
    if 0 not in c or 1 not in c or not (abs(c[0]) > 0.0 and np.isfinite(c[0])):
        return None
    t = c[1] / c[0] / math.sqrt(1 + a)
    t2 = c[2] / c[1] / math.sqrt((2 + a) / 2.0) if 2 in c and abs(c[1]) > 0.0 else np.nan
    u = (1.0 - t) / (1.0 + t)
    if not u.real > 0.0:                 # |t| >= 1: not a convergent single-beam series
        return {'t': t, 't2': t2, 'w0': np.nan, 'z0': np.nan}
    z_Rb = 0.5 * k * w_b ** 2
    return {'t': t, 't2': t2, 'w0': w_b * math.sqrt(u.real), 'z0': -z_Rb * u.imag}


# --------------------------------------------------------------------------------------------------------------
# Coherent-source prediction (electrons at rest at the waist)
# --------------------------------------------------------------------------------------------------------------

def fresnel_propagate_radial(source, rr, L, rho, s, k):
    """
    Fresnel field at distance s of the source S(r) e^{iL phi} on the plane s = 0, sampled on the radial grid `rr`
    (trapezoid), on the radii `rho`:

        E(rho) = (k/(i s)) (-i)^L e^{ik rho^2/(2s)} int S(r) e^{ik r^2/(2s)} J_L(k rho r/s) r dr,

    times e^{iL phi}, with the carrier e^{iks} left out as in lg_mode (the self-test checks the two agree).
    """
    w_r = np.full(len(rr), rr[1] - rr[0]); w_r[[0, -1]] *= 0.5
    integrand = source * np.exp(0.5j * k * rr ** 2 / s) * rr * w_r
    E = np.empty(len(rho), dtype=complex)
    for start in range(0, len(rho), _FRESNEL_CHUNK):        # bounded memory: chunk x len(rr) Bessel values
        chunk = rho[start:start + _FRESNEL_CHUNK]
        E[start:start + _FRESNEL_CHUNK] = scipy.special.jv(L, np.outer(chunk, k * rr / s)) @ integrand
    return (k / (1j * s)) * (-1j) ** L * np.exp(0.5j * k * rho ** 2 / s) * E


def coherent_source_prediction(r, phi, s, k, N, l_target, m, laser_p, laser_w0, R_beam):
    """
    Screen field of harmonic N from electrons at rest at the waist of an LG_p^m laser (m about lab +z), in the
    weak-field limit (theory/lg_decomposition_implementation.md, "Coherent-source prediction"): the coherent
    source is the laser mode to the N-th power over the beam cylinder, u_{p,m}(r, 0)^N e^{iNm phi} for r < R_beam,
    Fresnel-propagated to s, times each electron's own near-axis angular factor (r/s)^|dl| e^{i dl phi} that
    takes the source's N m to the fitted component's l_target (dl = l_target - N m). No free parameter; only the
    overall complex scale is undetermined. Returns (field on the screen points, dl).
    """
    dl = l_target - N * m
    r_max = min(R_beam, _SOURCE_RMAX_W0 * laser_w0) if R_beam else _SOURCE_RMAX_W0 * laser_w0
    # Enough points for the Bessel argument k rho r/s and the Fresnel phase k r^2/(2s) across the source.
    cycles = (k * r_max * r.max() / s + 0.5 * k * r_max ** 2 / s) / (2.0 * np.pi)
    rr = np.linspace(0.0, r_max, max(_SOURCE_MIN_POINTS, int(_SOURCE_POINTS_PER_CYCLE * cycles)))
    source = lg_mode(laser_p, m, rr, 0.0, 0.0, laser_w0, k).real ** N       # phi = 0, s = 0: real
    rho, index = np.unique(r, return_inverse=True)
    radial = fresnel_propagate_radial(source, rr, N * m, rho, s, k)[index]
    return radial * np.exp(1j * N * m * phi) * (r / s) ** abs(dl) * np.exp(1j * dl * phi), dl


def prediction_fidelity(pred, target, weights):
    """|<pred|E>|^2/(||pred||^2 ||E||^2): the share of the field's power that the prediction explains after the best
    complex scale (1 - fidelity is the residual eps^2 of that one-parameter fit)."""
    overlap = np.sum(weights * np.conj(pred) * target)
    return abs(overlap) ** 2 / (np.sum(weights * np.abs(pred) ** 2) * np.sum(weights * np.abs(target) ** 2))


# --------------------------------------------------------------------------------------------------------------
# Field components and angular-momentum bookkeeping
# --------------------------------------------------------------------------------------------------------------

def helicity_components(Ex, Ey):
    """(E_+, E_-): amplitudes on e_+- = (x + i y)/sqrt2 (lab x, y), i.e. spin +-1 about lab +z for
    exp(-i omega t) -- S_z ~ Im(Ex* Ey) = |E_+|^2 - |E_-|^2."""
    return (Ex - 1j * Ey) / math.sqrt(2.0), (Ex + 1j * Ey) / math.sqrt(2.0)


def component_field(name, Ex, Ey):
    E_plus, E_minus = helicity_components(Ex, Ey)
    return {'plus': E_plus, 'minus': E_minus, 'Ex': Ex, 'Ey': Ey}[name]


def component_l_values(name, J, l_override=None):
    """Azimuthal indices fitted for a component: J - s for helicity s; Ex/Ey mix both helicities, so both."""
    if l_override is not None:
        return [l_override]
    return {'plus': [J - 1], 'minus': [J + 1], 'Ex': [J - 1, J + 1], 'Ey': [J - 1, J + 1]}[name]


def total_angular_momentum(N, m, sigma):
    """J = N(m + sigma) per emitted photon (m, sigma about lab +z); None if sigma is not 0 or +-1."""
    for s in (-1, 0, 1):
        if abs(sigma - s) < _SPIN_TOL:
            return N * (m + s)
    return None


def azimuthal_spectrum(E_grid, radial_weights):
    """
    Circular grid only: power fraction per azimuthal index l, from an FFT over the N_phi - 1 unique phi columns
    (phi = 2pi duplicates phi = 0) of each ring, rings weighted by their trapezoidal r^2 weights. Returns
    (l_values, fractions).
    """
    unique = E_grid[:, :-1]
    M = unique.shape[1]
    a_l = np.fft.fft(unique, axis=1) / M          # a_l(r) of E(r, phi) = sum_l a_l(r) e^{i l phi}
    power = np.sum(radial_weights[:, None] * np.abs(a_l) ** 2, axis=0)
    l_values = np.round(np.fft.fftfreq(M, 1.0 / M)).astype(int)
    total = power.sum()
    return l_values, (power / total if total > 0 else power * np.nan)


def _estimate_w0(target, weights, r, l, s, k):
    """
    Starting waist for the default scan: the screen width of a p = 0 LG mode with the same second moment,
    w(s)^2 = 2<r^2>/(|l|+1), inverted for w0 on the far branch (s > z_R) as in plot_observables'
    estimate_gaussian_waist; s = 0 (incident field) gives w0 = w(s) directly.
    """
    power = weights * np.abs(target) ** 2
    w_s_sq = 2.0 * np.sum(power * r ** 2) / np.sum(power) / (abs(l) + 1)
    if s == 0.0:
        return math.sqrt(w_s_sq)
    B = 2.0 * s / k                    # = lambda s / pi
    disc = w_s_sq ** 2 - 4.0 * B ** 2
    if disc < 0.0:                     # narrower than any Gaussian can be at s: take the minimum-width waist
        return math.sqrt(B)
    return math.sqrt(2.0 * B ** 2 / (w_s_sq + math.sqrt(disc)))


def _dominant_share(coeffs):
    """Largest single |c|^2 over sum |c|^2 (1 = the whole fitted field in one mode)."""
    c_sq = np.abs(coeffs) ** 2
    total = np.sum(c_sq)
    return float(c_sq.max() / total) if np.isfinite(total) and total > 0.0 else -np.inf


def _choose_w0(fits):
    """
    The reported (w0, z0, fit) of a scan (items of `fits`; z0 = 0 unless the waist position is scanned): the
    trial at which one mode holds the largest share of sum|c|^2, among the interpretable fits
    (|captured - (1 - eps^2)| < _CAPTURED_TOL, i.e. the modes are orthonormal on the screen; all fits if none
    is). For a single LG_p^l beam that is its own waist (and waist position), whatever p and l. The smallest-eps
    w0 is not used: with several p free and a w0-independent residual, eps(w0) is flat and its minimum is set by
    noise.
    """
    interpretable = [item for item in fits
                     if abs(item[-1]['captured'] - (1.0 - item[-1]['eps'] ** 2)) < _CAPTURED_TOL]
    return max(interpretable or fits, key=lambda item: _dominant_share(item[-1]['coeffs']))


def _radius_containing(fraction, power, r):
    """Smallest screen radius containing `fraction` of the (area-weighted) power."""
    order = np.argsort(r)
    cumulative = np.cumsum(power[order])
    return r[order][min(np.searchsorted(cumulative, fraction * cumulative[-1]), len(r) - 1)]


# --------------------------------------------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------------------------------------------

def screen_lab_xy(x_local, y_local, local_rotation):
    """Lab-frame transverse coordinates of local screen points (z_local = 0 plane; the screen faces +-Oz)."""
    lab = local_rotation @ np.vstack([x_local, y_local, np.zeros_like(x_local)])
    return lab[0], lab[1]


def _grid_step_au(detector_type, config_path):
    """The screen grid's resolution in atomic units: the r^2 step (equal-area rings) of a circular grid, the larger
    of dx, dy of a rectangular one. Compared with w^2 and w respectively (_MAX_STEP_OVER_WIDTH)."""
    if detector_type == 'circular':
        R_min, R_max = (read_length(f'circular_detector_{key}', config_path, 'a.u.') for key in ('R_min', 'R_max'))
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        return (R_max ** 2 - R_min ** 2) / max(N_R - 1, 1)
    steps = []
    for axis in ('x', 'y'):
        lo, hi = (read_length(f'rectangular_detector_{axis}_{end}', config_path, 'a.u.') for end in ('min', 'max'))
        steps.append(abs(hi - lo) / max(int(read_config_value(f'rectangular_detector_N{axis}', config_path)[0]) - 1, 1))
    return max(steps)


def _check_screen_along_z(local_rotation, config_path):
    normal = local_rotation[:, 2]
    if abs(normal[0]) > _TRANSVERSE_TOL or abs(normal[1]) > _TRANSVERSE_TOL:
        theta = read_config_value('detector_direction_theta', config_path)
        raise ValueError(f"detector_direction_theta = {' '.join(theta)}: the screen does not face +-Oz. The LG "
                         "decomposition needs a screen transverse to lab z (same restriction as plot_observables).")


def _source_model(run_dir, config_path):
    """
    Laser and beam parameters for coherent_source_prediction, or a string saying why it doesn't apply: it needs an
    LG laser and electrons at rest (average_px/py/pz = 0). The beam radius is the one run_log.txt records (it may
    come from the LG table, beam_cylinder_radius < 0); None if absent, and the source is then capped at
    _SOURCE_RMAX_W0 laser waists.
    """
    if read_config_value('laser_type', config_path)[0] != 'laguerre_gauss':
        return "needs laser_type laguerre_gauss"
    if any(float(read_config_value(f'average_p{axis}', config_path)[0]) != 0.0 for axis in 'xyz'):
        return "electrons not at rest (average_p != 0): not implemented"
    R_beam = None
    with open(os.path.join(run_dir, 'run_log.txt')) as run_log:
        for line in run_log:
            if line.strip().startswith('cylinder radius'):
                R_beam = float(line.split()[2])
                break
    return {'p': int(read_config_value('laser_lg_p', config_path)[0]),
            'w0': read_length('laser_lg_w0', config_path, 'a.u.'), 'R_beam': R_beam}


def _parse_w0_values(w0_spec, unit, config_path):
    """'10,20,30' in `unit` (any convert_unit_to_number unit: lambda, w0, a.u.) -> atomic units."""
    values = [float(t) * convert_unit_to_number(unit, config_path) for t in w0_spec.split(',') if t.strip()]
    if not values or any(v <= 0 for v in values):
        raise ValueError(f"--w0 needs one or more positive comma-separated values, got '{w0_spec}'")
    return values


def _parse_z0_values(z0_spec, unit, config_path):
    """'-1000,0,1000' in `unit` -> atomic units (any sign: the basis waist's position along the propagation)."""
    values = [float(t) * convert_unit_to_number(unit, config_path) for t in z0_spec.split(',') if t.strip()]
    if not values:
        raise ValueError(f"--z0 needs one or more comma-separated values, got '{z0_spec}'")
    return values


# --------------------------------------------------------------------------------------------------------------
# Run analysis
# --------------------------------------------------------------------------------------------------------------

def decompose_run(field_filepath, components, harmonic=None, l_override=None, p_min=0, p_max=3, w0_spec=None,
                  w0_unit='lambda', n_w0=_N_W0_DEFAULT, scan_z0=False, z0_spec=None, z0_unit='lambda',
                  n_z0=_N_Z0_DEFAULT):
    """
    Fits every selected component at every (selected) frequency for every trial w0 and, with scan_z0 (or an
    explicit z0_spec), every trial waist position z0 of the basis (s -> s - z0; default z0 = 0, the origin).
    Returns (table, summary, context): `table` one row per (i_omega, component, w0, z0); `summary` one row per
    (i_omega, component) with the reported (w0, z0) (_choose_w0), the smallest-eps one, and the helicity/azimuthal
    diagnostics; `context` the run-level parameters.
    """
    run_dir = os.path.dirname(field_filepath)
    is_incident = os.path.basename(field_filepath) == "incident_field.dat"
    config_path = field_config_path(field_filepath)    # incident_config.cfg for a recent incident_field.dat
    data = pd.read_csv(field_filepath, sep=' ', comment='#')

    detector_type, _ = read_config_value('detector_type', config_path)
    if detector_type not in ('rectangular', 'circular'):
        raise ValueError(f"detector_type '{detector_type}': only flat rectangular/circular screens are supported.")

    # Screen points in atomic units, in the frame whose +z the angular momentum refers to: the lab frame for the
    # emitted field, the laser's own frame (grid and components) for incident_field.dat.
    x_local, y_local, _, _ = get_screen_coordinates(field_filepath, config_path)
    axes_to_au = convert_unit_to_number(DETECTOR_AXES_UNIT, config_path)
    x_local, y_local = x_local * axes_to_au, y_local * axes_to_au
    m = int(read_config_value('laser_lg_l', config_path)[0])
    sigma = _read_laser_spin(config_path)
    if is_incident:
        x, y = x_local, y_local
        # The incident plane: the detector distance downstream of the waist (incident_config.cfg), or the waist
        # itself for runs before 2026-10-27.
        s = (abs(read_length(f'{detector_type}_detector_distance', config_path, 'a.u.'))
             if is_incident_screen_config(config_path) else 0.0)
    else:
        local_rotation = get_detector_local_rotation(config_path)
        _check_screen_along_z(local_rotation, config_path)
        x, y = screen_lab_xy(x_local, y_local, local_rotation)
        s = abs(read_length(f'{detector_type}_detector_distance', config_path, 'a.u.'))
        laser_dir, laser_label = _read_laser_direction(config_path)
        if abs(laser_dir[0]) > _TRANSVERSE_TOL or abs(laser_dir[1]) > _TRANSVERSE_TOL:
            print(f"WARNING: {laser_label} is not along +-Oz: J = N(m + sigma) about lab z does not apply.")
        n_sign = int(np.sign(laser_dir[2]))
        m *= n_sign
        sigma *= n_sign
    r = np.hypot(x, y)
    phi = np.arctan2(y, x)
    weights = get_screen_area_weights_au(detector_type, config_path)

    radial_weights = None
    if detector_type == 'circular':
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        grid_shape = (N_R, N_phi)
        radial_weights = weights.reshape(grid_shape).sum(axis=1)     # ring areas, (1/2) w_{r^2} 2 pi
        # Lab azimuth increases with the local grid's phi index on a +Oz screen and decreases on a -Oz one
        # (local y = -lab y), so the FFT index is l_lab * det(R_perp).
        phi_orientation = 1 if is_incident else int(np.sign(np.linalg.det(local_rotation[:2, :2])))

    grid_step = _grid_step_au(detector_type, config_path)
    source_model = None if is_incident else _source_model(run_dir, config_path)
    fundamental_frequency_au = _read_fundamental_frequency_au(run_dir)
    lambda_L = convert_unit_to_number('lambda', config_path)
    p_values = list(range(p_min, p_max + 1))

    table_rows, summary_rows = [], []
    for i_omega, subset in data.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_ratio = subset['omega'].iloc[0]
        if len(subset) != len(weights):
            raise ValueError(f"i_omega={i_omega}: {len(subset)} screen points, {len(weights)} expected from the config.")
        N = 1 if is_incident else int(round(omega_ratio))
        if harmonic is not None and N != harmonic:
            continue
        if not is_incident and abs(omega_ratio - N) > 1e-6:
            print(f"WARNING: i_omega={i_omega}: omega/omega_1 = {omega_ratio:.6g} is not a harmonic; using N = {N}.")
        k = omega_ratio * fundamental_frequency_au / C_LIGHT
        J = total_angular_momentum(N, m, sigma)
        if J is None and l_override is None:
            raise ValueError(f"laser spin sigma = {sigma:+.6g} is neither 0 nor +-1, so J = N(m + sigma) is not an "
                             "integer; pass --l explicitly.")

        Ex, Ey, _ = _get_E_total(subset)
        E_plus, E_minus = helicity_components(Ex, Ey)
        P_plus, P_minus = (np.sum(weights * np.abs(E) ** 2) for E in (E_plus, E_minus))
        helicity_fraction = {'plus': P_plus / (P_plus + P_minus), 'minus': P_minus / (P_plus + P_minus)}

        for name in components:
            target = component_field(name, Ex, Ey)
            l_values = component_l_values(name, J, l_override)
            modes = [(p, l) for l in l_values for p in p_values]
            power = weights * np.abs(target) ** 2
            summary = {'i_omega': i_omega, 'omega': omega_ratio, 'N': N, 'component': name, 'J': J,
                       'modes': modes, 'helicity_fraction': helicity_fraction.get(name, np.nan), 'fit': None,
                       'sph_phase': np.nan, 'azimuthal': None, 'prediction': None}
            summary_rows.append(summary)
            if not np.sum(power) > 0.0:
                continue                   # e.g. the opposite helicity of the incident beam: nothing to fit

            if w0_spec:
                w0_list = _parse_w0_values(w0_spec, w0_unit, config_path)
            else:
                w0_est = _estimate_w0(target, weights, r, l_values[0], s, k)
                w0_list = list(w0_est * np.geomspace(1.0 / _W0_SCAN_FACTOR, _W0_SCAN_FACTOR, n_w0))
            if z0_spec:
                z0_list = _parse_z0_values(z0_spec, z0_unit, config_path)
            elif scan_z0:
                # +-_Z0_SCAN_ZR Rayleigh ranges of the scan's central waist, odd count so that z0 = 0 is included.
                z_R_centre = 0.5 * k * math.sqrt(w0_list[0] * w0_list[-1]) ** 2
                z0_list = list(np.linspace(-_Z0_SCAN_ZR * z_R_centre, _Z0_SCAN_ZR * z_R_centre, 2 * (n_z0 // 2) + 1))
            else:
                z0_list = [0.0]

            fits = []
            for z0 in z0_list:
                for w0 in w0_list:
                    fit = fit_modes(target, build_basis(modes, r, phi, s - z0, w0, k), weights)
                    fits.append((w0, z0, fit))
                    row = {'i_omega': i_omega, 'omega': omega_ratio, 'N': N, 'component': name,
                           'l': '/'.join(str(l) for l in l_values), 'w0_au': w0, 'w0_lambdaL': w0 / lambda_L,
                           'z0_au': z0, 'z0_lambdaL': z0 / lambda_L, 'eps': fit['eps'], 'captured': fit['captured'],
                           'cond': fit['cond'], 'rank': fit['rank'], 'norm_sq': np.sum(power),
                           'dominant_share': _dominant_share(fit['coeffs'])}
                    for (p, l), c in zip(modes, fit['coeffs']):
                        row[f'c_p{p}_l{l}_re'], row[f'c_p{p}_l{l}_im'] = c.real, c.imag
                    series = series_waist(fit['coeffs'], modes, w0, k)
                    if series is not None:
                        row['w0_series_lambdaL'] = series['w0'] / lambda_L
                        row['z0_series_lambdaL'] = (z0 + series['z0']) / lambda_L
                    table_rows.append(row)

            w0_min_eps, z0_min_eps, fit_min_eps = min(fits, key=lambda item: item[-1]['eps'])
            best_w0, best_z0, best = _choose_w0(fits)
            series = series_waist(best['coeffs'], modes, best_w0, k)
            if series is not None:
                series['z0'] += best_z0          # relative to the basis waist plane -> from the origin
            summary.update(w0_best=best_w0, z0_best=best_z0, w0_list=w0_list, z0_list=z0_list, fit=best,
                           w0_min_eps=w0_min_eps, z0_min_eps=z0_min_eps, fit_min_eps=fit_min_eps,
                           edge_of_scan=len(w0_list) > 1 and best_w0 in (w0_list[0], w0_list[-1]),
                           z0_edge_of_scan=len(z0_list) > 1 and best_z0 in (z0_list[0], z0_list[-1]),
                           series=series)
            s_waist = s - best_z0
            w_s = best_w0 * math.sqrt(1.0 + (s_waist / (0.5 * k * best_w0 ** 2)) ** 2)
            summary['step_over_width'] = grid_step / (w_s ** 2 if detector_type == 'circular' else w_s)
            if s > 0.0:
                # Exact spherical phase k(sqrt(s^2+r^2) - s) minus the paraxial k r^2/(2s) the LG modes assume.
                r99 = _radius_containing(0.99, power, r)
                summary['sph_phase'] = k * r99 ** 4 / (8.0 * s ** 3)
            if radial_weights is not None and name in ('plus', 'minus'):
                l_fft, fractions = azimuthal_spectrum(target.reshape(grid_shape), radial_weights)
                summary['azimuthal'] = (phi_orientation * l_fft, fractions)
            if isinstance(source_model, str):
                summary['prediction'] = source_model                 # why it doesn't apply
            elif source_model is not None and name in ('plus', 'minus') and len(l_values) == 1 and s > 0.0:
                pred, dl = coherent_source_prediction(r, phi, s, k, N, l_values[0], m, source_model['p'],
                                                      source_model['w0'], source_model['R_beam'])
                pred_power = weights * np.abs(pred) ** 2
                r_contained = _radius_containing(_PREDICTION_CONTAINED, pred_power, r)
                summary['prediction'] = {'fidelity': prediction_fidelity(pred, target, weights), 'dl': dl,
                                         'r_contained': r_contained,
                                         'outside': power[r > r_contained].sum() / power.sum()}

    context = {'run_dir': run_dir, 'config_path': config_path, 'is_incident': is_incident, 'm': m,
               'sigma': sigma, 's': s, 'lambda_L': lambda_L, 'detector_type': detector_type}
    return pd.DataFrame(table_rows), summary_rows, context


def _summary_lines(summary_rows, context, table_name):
    lambda_L = context['lambda_L']
    lines = [
        "  theory/lg_decomposition_implementation.md: weighted complex least squares of one field component on "
        "unit-norm LG modes",
        "  u_{p,l}(r, phi, s; w0) (waist at the origin, s = screen distance; incident: laser frame), one w0 "
        "for all modes;",
        "  reported at the trial w0 where one mode holds the largest share of sum|c|^2 (the beam's own waist for a "
        "single LG_p^l; the",
        "  smallest-eps w0 is a basis label only). plus/minus = (Ex -+ i Ey)/sqrt2, spin +-1 about lab +z; fitted "
        "l = J - spin,",
        f"  J = N(m + sigma), m = {context['m']}, sigma = {context['sigma']:+.6g} (about lab +z). eps = relative "
        "area-weighted residual;",
        "  captured = sum|c|^2/sum A|E|^2 (= 1 - eps^2 if the modes are orthonormal on the screen). Coefficients "
        "carry the common phase e^{iks}.",
        f"  Lengths in lambda_L = {lambda_L * BOHR_RADIUS_M * 1e9:.6g} nm; s = {context['s'] / lambda_L:.6g} lambda_L. "
        f"Full w0 scan: {table_name}.",
    ]
    for row in summary_rows:
        fit = row['fit']
        lines.append("")
        lines.append(f"  i_omega {row['i_omega']} (omega/omega_1 = {row['omega']:.6f}, N = {row['N']}), component "
                     f"{row['component']}, J = {row['J']}, helicity power fraction = {row['helicity_fraction']:.4g}")
        if fit is None:
            lines.append("    zero power in this component: not fitted")
            continue
        dominant = int(np.argmax(np.abs(fit['coeffs'])))
        p_dom, l_dom = row['modes'][dominant]
        if len(row['z0_list']) > 1:
            lines.append(f"    basis waist position z0 = {row['z0_best'] / lambda_L:.6g} lambda_L along propagation "
                         f"(scan {row['z0_list'][0] / lambda_L:.4g}..{row['z0_list'][-1] / lambda_L:.4g}, "
                         f"{len(row['z0_list'])} values{'; AT SCAN EDGE' if row['z0_edge_of_scan'] else ''}), with")
        lines.append(f"    w0 = {row['w0_best'] / lambda_L:.6g} lambda_L, where (p,l) = ({p_dom},{l_dom}) holds "
                     f"{_dominant_share(fit['coeffs']):.4g} of sum|c|^2 (scan {row['w0_list'][0] / lambda_L:.4g}"
                     f"..{row['w0_list'][-1] / lambda_L:.4g}, {len(row['w0_list'])} values"
                     f"{'; AT SCAN EDGE' if row['edge_of_scan'] else ''}): eps = {fit['eps']:.4e}, captured = "
                     f"{fit['captured']:.6f}, cond = {fit['cond']:.3e}, status: {fit['status']}")
        lines.append(f"    smallest eps = {row['fit_min_eps']['eps']:.4e} at w0 = {row['w0_min_eps'] / lambda_L:.6g} "
                     f"lambda_L" + (f", z0 = {row['z0_min_eps'] / lambda_L:.6g} lambda_L" if len(row['z0_list']) > 1
                                    else "") + " (a basis label, not the beam's waist)")
        step_ratio, step_max = row['step_over_width'], _MAX_STEP_OVER_WIDTH[context['detector_type']]
        step_label = 'r^2 step / w(s)^2' if context['detector_type'] == 'circular' else 'grid step / w(s)'
        if step_ratio > step_max:
            lines.append(f"    WARNING: UNDER-RESOLVED, {step_label} = {step_ratio:.3g} > {step_max:g} at the reported "
                         "w0: "
                         "the grid doesn't sample the beam, nothing below is meaningful (a finer or smaller screen "
                         "is needed)")
        else:
            lines.append(f"    resolution: {step_label} = {step_ratio:.3g} (flagged above {step_max:g})")
        series = row.get('series')
        if series is not None:
            t, t2 = series['t'], series['t2']
            # Near the beam's own waist c_1, c_2 -> 0 and c_2/c_1 is a ratio of residuals; the w0_series columns
            # across the scan are then the single-beam check.
            check = (f"t from c_2/c_1 = {abs(t2):.4f} e^{{i {np.angle(t2) / np.pi:+.3f} pi}} (equal for a single beam)"
                     if abs(t) >= _MIN_T_FOR_CHECK else
                     f"|t| < {_MIN_T_FOR_CHECK:g}: basis ~ beam, c_2/c_1 not meaningful (single-beam check: "
                     "w0_series_lambdaL across the scan)")
            lines.append(f"    series waist (single LG_0^l beam, from c_1/c_0): w0 = {series['w0'] / lambda_L:.6g} "
                         f"lambda_L, waist at z0 = {series['z0'] / lambda_L:.6g} lambda_L along propagation; "
                         f"t = {abs(t):.4f} e^{{i {np.angle(t) / np.pi:+.3f} pi}}, {check}")
        if np.isfinite(row['sph_phase']):
            lines.append(f"    spherical - paraxial phase at the 99% power radius: {row['sph_phase']:.3e} rad")
        if row['azimuthal'] is not None:
            l_fft, fractions = row['azimuthal']
            top = np.argsort(fractions)[::-1][:3]
            expected = component_l_values(row['component'], row['J'])[0] if row['J'] is not None else None
            at_expected = fractions[l_fft == expected].sum() if expected is not None else np.nan
            lines.append(f"    azimuthal power: l = {expected}: {at_expected:.4g}; largest: "
                         + ', '.join(f"l={l_fft[i]}: {fractions[i]:.4g}" for i in top))
        prediction = row.get('prediction')
        if isinstance(prediction, str):
            lines.append(f"    coherent-source prediction: {prediction}")
        elif prediction is not None:
            lines.append(f"    coherent-source prediction (rest, weak field; laser mode^N, angular factor "
                         f"theta^{abs(prediction['dl'])}): fidelity = {prediction['fidelity']:.4f}, measured power "
                         f"outside its {100 * _PREDICTION_CONTAINED:g}% radius "
                         f"({prediction['r_contained'] / lambda_L:.4g} lambda_L) = {prediction['outside']:.4g}")
        norm_sq = np.sum(np.abs(fit['coeffs']) ** 2)
        for (p, l), c in zip(row['modes'], fit['coeffs']):
            lines.append(f"    p={p:2d} l={l:+3d}  c = {c.real:+.6e} {c.imag:+.6e}i  |c| = {abs(c):.6e}  "
                         f"|c|^2/sum|c|^2 = {abs(c) ** 2 / norm_sq if norm_sq > 0 else np.nan:.4g}")
    return lines


def _plot(table, summary_rows, context):
    png_dir = os.path.join(context['run_dir'], "png_folder", MODULE_NAME, "lg_decomposition",
                           "incident" if context['is_incident'] else "emitted")
    os.makedirs(png_dir, exist_ok=True)
    for i_omega, subset in table.groupby('i_omega'):
        fig, (ax_eps, ax_c) = plt.subplots(1, 2, figsize=(11, 4.2), layout='constrained')
        for row in (r for r in summary_rows if r['i_omega'] == i_omega and r['fit'] is not None):
            curve = subset[(subset['component'] == row['component']) & (subset['z0_au'] == row['z0_best'])]
            label = f"{row['component']} (l = {curve['l'].iloc[0]})"
            ax_eps.semilogx(curve['w0_lambdaL'], curve['eps'], 'o-', ms=3, label=label)
            c_sq = np.abs(row['fit']['coeffs']) ** 2
            ticks = np.arange(len(row['modes']))
            ax_c.plot(ticks, c_sq / c_sq.sum() if c_sq.sum() > 0 else c_sq, 'o-', label=label)
            ax_c.set_xticks(ticks, [f"{p},{l}" for p, l in row['modes']], rotation=45)   # same p order for all
        ax_eps.set_xlabel(r'trial $w_0$ [$\lambda_L$]')
        ax_eps.set_ylabel(r'$\epsilon(w_0)$')
        ax_eps.grid(True, which='both', alpha=0.3)
        ax_eps.legend()
        ax_c.set_xlabel('(p, l)')
        ax_c.set_ylabel(r'$|c_{p,l}|^2/\sum|c|^2$ at reported $w_0$')
        ax_c.set_yscale('log')
        ax_c.grid(True, alpha=0.3)
        fig.suptitle(f"LG decomposition, $\\omega$ index {i_omega} ($\\omega/\\omega_1$ = {subset['omega'].iloc[0]:.4g})"
                     + (" (incident beam)" if context['is_incident'] else ""))
        output = os.path.join(png_dir, f"lg_decomposition_omega{i_omega}.png")
        fig.savefig(output, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Successfully saved plot to {output}")

    for row in summary_rows:
        if row['fit'] is None:
            continue
        if row['helicity_fraction'] < _MIN_HELICITY_FRACTION_PLOT:     # NaN for Ex/Ey: always plotted
            print(f"Skipping the coefficient plot of i_omega={row['i_omega']} {row['component']}: helicity power "
                  f"fraction {row['helicity_fraction']:.3g} < {_MIN_HELICITY_FRACTION_PLOT:g}.")
            continue
        _plot_coefficients(table, row, context, png_dir)
        if len(row['z0_list']) > 1:
            _plot_w0_z0_map(table, row, context, png_dir)


def _plot_w0_z0_map(table, row, context, png_dir):
    """
    z0 scan only: the dominant mode's share of sum|c|^2 and eps over the (w0, z0) grid of one component at one
    frequency, the reported point marked. Grid cells that are not interpretable (|captured - (1 - eps^2)| >=
    _CAPTURED_TOL: the modes are not orthonormal on the screen) are masked.
    """
    scan = table[(table['i_omega'] == row['i_omega']) & (table['component'] == row['component'])]
    w0 = np.sort(scan['w0_lambdaL'].unique())
    z0 = np.sort(scan['z0_lambdaL'].unique())
    grids = {}
    for name in ('dominant_share', 'eps', 'captured'):
        grids[name] = scan.pivot_table(index='z0_lambdaL', columns='w0_lambdaL', values=name).loc[z0, w0].to_numpy()
    bad = np.abs(grids['captured'] - (1.0 - grids['eps'] ** 2)) >= _CAPTURED_TOL

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), layout='constrained')
    for ax, name, label in ((axes[0], 'dominant_share', r'dominant mode: $\max|c|^2/\sum|c|^2$'),
                            (axes[1], 'eps', r'$\epsilon$')):
        values = np.ma.masked_where(bad, grids[name])
        mesh = ax.pcolormesh(w0, z0, values, shading='nearest', cmap='viridis' if name == 'dominant_share'
                             else 'magma_r')
        fig.colorbar(mesh, ax=ax, label=label)
        ax.plot(row['w0_best'] / context['lambda_L'], row['z0_best'] / context['lambda_L'], 'r+', ms=12, mew=2)
        ax.set_xscale('log')
        ax.set_xlabel(r'basis $w_0$ [$\lambda_L$]')
        ax.set_ylabel(r'basis waist position $z_0$ [$\lambda_L$]')
    fig.suptitle(f"LG basis scan, {row['component']}, $\\omega$ index {row['i_omega']} ($\\omega/\\omega_1$ = "
                 f"{row['omega']:.4g}); red +: reported (gray: not interpretable)"
                 + (" (incident beam)" if context['is_incident'] else ""))
    for ax in axes:
        ax.set_facecolor('0.85')
    output = os.path.join(png_dir, f"lg_w0_z0_map_omega{row['i_omega']}_{row['component']}.png")
    fig.savefig(output, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Successfully saved plot to {output}")


def _plot_coefficients(table, row, context, png_dir):
    """
    One component at one frequency, the whole coefficient set: |c|^2/||E||^2 vs w0 per mode (with eps^2 and
    captured), the Argand diagram of c/||E|| at the reported w0 (_choose_w0) and arg(c/c_dominant) vs w0.
    ||E||^2 = sum A|E|^2, so |c|^2/||E||^2 is the mode's power fraction when the modes are orthonormal on the
    screen. Phases are relative to the mode that dominates at the reported w0, which removes the arbitrary common
    e^{iks}. With a z0 scan, the w0 curves are at the reported z0.
    """
    curve = table[(table['i_omega'] == row['i_omega']) & (table['component'] == row['component'])
                  & (table['z0_au'] == row['z0_best'])]
    curve = curve.sort_values('w0_au')
    w0 = curve['w0_lambdaL'].to_numpy()
    norm = np.sqrt(curve['norm_sq'].to_numpy())
    coeffs = np.stack([(curve[f'c_p{p}_l{l}_re'] + 1j * curve[f'c_p{p}_l{l}_im']).to_numpy() / norm
                       for p, l in row['modes']], axis=1)                  # (n_w0, n_modes), c/||E||
    best_index = int(np.argmin(np.abs(w0 - row['w0_best'] / context['lambda_L'])))
    dominant = int(np.argmax(np.abs(coeffs[best_index])))
    rotated_best = coeffs[best_index] * np.exp(-1j * np.angle(coeffs[best_index, dominant]))
    single_l = len({l for _, l in row['modes']}) == 1
    labels = [f"p={p}" if single_l else f"p={p}, l={l}" for p, l in row['modes']]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, (ax_pow, ax_arg, ax_phase) = plt.subplots(1, 3, figsize=(16, 4.8), layout='constrained')
    for j, label in enumerate(labels):
        color = colors[j % len(colors)]
        ax_pow.semilogx(w0, np.abs(coeffs[:, j]) ** 2, 'o-', ms=2.5, color=color, label=label)
        # Phases are noise where the mode carries almost nothing; mask below 1e-6 of the power.
        phase = np.angle(coeffs[:, j] / coeffs[:, dominant])
        phase[np.abs(coeffs[:, j]) ** 2 < 1e-6] = np.nan
        ax_phase.semilogx(w0, phase / np.pi, 'o', ms=2.5, color=color, label=label)
        z = rotated_best[j]
        ax_arg.plot([0, z.real], [0, z.imag], '-', color=color, lw=1.5)
        ax_arg.plot(z.real, z.imag, 'o', color=color, label=label)
        ax_arg.annotate(label, (z.real, z.imag), textcoords='offset points', xytext=(5, 5), fontsize=8)
    ax_pow.semilogx(w0, curve['eps'] ** 2, 'k--', lw=1.2, label=r'$\epsilon^2$')
    ax_pow.semilogx(w0, curve['captured'], 'k:', lw=1.2, label=r'captured $=\sum|c|^2/\|E\|^2$')
    for ax in (ax_pow, ax_phase):
        ax.axvline(w0[best_index], color='gray', lw=0.8, alpha=0.7)
        ax.set_xlabel(r'trial $w_0$ [$\lambda_L$]')
        ax.grid(True, which='both', alpha=0.3)
    ax_pow.set_ylabel(r'$|c_{p,l}|^2/\|E\|^2$')
    ax_pow.set_yscale('log')
    ax_pow.set_ylim(bottom=max(1e-8, ax_pow.get_ylim()[0]))
    ax_pow.legend(fontsize=8)
    ax_pow.set_title('power fraction per mode')

    limit = 1.1 * max(np.abs(rotated_best).max(), 1e-30)
    circle = np.exp(1j * np.linspace(0, 2 * np.pi, 200))
    ax_arg.plot(circle.real, circle.imag, color='gray', lw=0.6, alpha=0.5)         # |c| = ||E||: the whole field
    ax_arg.axhline(0, color='gray', lw=0.6)
    ax_arg.axvline(0, color='gray', lw=0.6)
    ax_arg.set_xlim(-limit, limit)
    ax_arg.set_ylim(-limit, limit)
    ax_arg.set_aspect('equal')
    ax_arg.set_xlabel(r'Re $c/\|E\|$')
    ax_arg.set_ylabel(r'Im $c/\|E\|$')
    ax_arg.grid(True, alpha=0.3)
    ax_arg.set_title(f"at reported $w_0$ = {w0[best_index]:.4g} $\\lambda_L$, rotated: {labels[dominant]} real")

    ax_phase.set_ylabel(rf'arg$(c/c_{{\rm {labels[dominant].replace(" ", "")}}})/\pi$')
    ax_phase.set_ylim(-1.05, 1.05)
    ax_phase.set_title('relative phase')

    fig.suptitle(f"LG coefficients, {row['component']} (l = {curve['l'].iloc[0]}), $\\omega$ index {row['i_omega']} "
                 f"($\\omega/\\omega_1$ = {row['omega']:.4g})" + (" (incident beam)" if context['is_incident'] else ""))
    output = os.path.join(png_dir, f"lg_coefficients_omega{row['i_omega']}_{row['component']}.png")
    fig.savefig(output, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Successfully saved plot to {output}")


def run(field_filepath, args):
    table, summary_rows, context = decompose_run(
        field_filepath, args.component, harmonic=args.harmonic, l_override=args.l, p_min=args.p_min,
        p_max=args.p_max, w0_spec=args.w0, w0_unit=args.w0_unit, n_w0=args.n_w0, scan_z0=args.scan_z0,
        z0_spec=args.z0, z0_unit=args.z0_unit, n_z0=args.n_z0)
    if table.empty:
        print("Nothing to fit (no frequency matched --harmonic?).")
        return

    suffix = "_incident" if context['is_incident'] else ""
    table_path = os.path.join(context['run_dir'], f"lg_decomposition{suffix}.dat")
    table.to_csv(table_path, sep=' ', index=False, float_format='%.10e')
    print(f"Wrote the full scan to {table_path}")

    lines = _summary_lines(summary_rows, context, os.path.basename(table_path))
    print(_RUN_LOG_HEADER + ":")
    print('\n'.join(lines))
    header = _RUN_LOG_HEADER + (" -- incident beam" if context['is_incident'] else "")
    _write_run_log_section(context['run_dir'], header, lines)
    _plot(table, summary_rows, context)


# --------------------------------------------------------------------------------------------------------------
# Self-test (no run needed)
# --------------------------------------------------------------------------------------------------------------

def self_test():
    """
    1. Normalization/orthogonality of lg_mode on a fine synthetic circular grid, at the waist and at s = 3 z_R.
    2. Recovery of known coefficients (l != 0, both helicities populated) on a +Oz and a -Oz screen, with the field
       synthesized directly in lab coordinates (x_lab, y_lab) = (x_loc, +-y_loc), through the same
       screen_lab_xy/helicity_components/fit path the run analysis uses.
    Returns True if every check passes.
    """
    ok = True

    def check(label, value, tol):
        nonlocal ok
        passed = abs(value) <= tol
        ok &= passed
        print(f"  [{'PASS' if passed else 'FAIL'}] {label}: {value:.3e} (tol {tol:.0e})")

    k, w0 = 2.0 * np.pi, 3.0          # wavelength 1, w0 = 3 wavelengths
    z_R = 0.5 * k * w0 ** 2

    def circular_grid(R_max, N_R, N_phi):
        r_vals = np.sqrt(np.linspace(0.0, R_max ** 2, N_R))
        phi_vals = np.linspace(0.0, 2.0 * np.pi, N_phi)
        r_g, phi_g = np.meshgrid(r_vals, phi_vals, indexing='ij')
        w_r = np.full(N_R, R_max ** 2 / (N_R - 1)); w_r[[0, -1]] *= 0.5
        w_p = np.full(N_phi, 2.0 * np.pi / (N_phi - 1)); w_p[[0, -1]] *= 0.5
        return r_g.ravel(), phi_g.ravel(), (0.5 * w_r[:, None] * w_p[None, :]).ravel()

    # The trapezoid in u = r^2 is O(h^2) for l = 0, 1 (|u|^2 has a nonzero u-slope at u = 0; ~ (2h/w^2)^2/12),
    # spectrally accurate for |l| >= 2, so this grid is much finer in r than a typical run's.
    print("1. Normalization and orthogonality (int over a screen of radius 6 w(s), 6000 rings):")
    for s in (0.0, 3.0 * z_R):
        w_s = w0 * math.sqrt(1.0 + (s / z_R) ** 2)
        r, phi, A = circular_grid(6.0 * w_s, N_R=6000, N_phi=129)
        for p, l in ((0, 0), (1, 0), (2, 3), (0, -2), (4, 1)):
            u = lg_mode(p, l, r, phi, s, w0, k)
            check(f"s = {s / z_R:g} z_R, (p,l) = ({p},{l}): int|u|^2 - 1", np.sum(A * np.abs(u) ** 2) - 1.0, 1e-4)
        for (p1, l1), (p2, l2) in (((0, 2), (1, 2)), ((0, 2), (0, -2)), ((1, 1), (3, 1))):
            overlap = np.sum(A * np.conj(lg_mode(p1, l1, r, phi, s, w0, k)) * lg_mode(p2, l2, r, phi, s, w0, k))
            check(f"s = {s / z_R:g} z_R, <{p1},{l1}|{p2},{l2}>", abs(overlap), 1e-4)

    print("2. Coefficient recovery on +Oz and -Oz screens (J = 3: E_+ has l = 2, E_- has l = 4):")
    s = 3.0 * z_R
    w_s = w0 * math.sqrt(1.0 + (s / z_R) ** 2)
    r_loc, phi_loc, A = circular_grid(6.0 * w_s, N_R=256, N_phi=65)
    x_loc, y_loc = r_loc * np.cos(phi_loc), r_loc * np.sin(phi_loc)
    c_plus = {(0, 2): 0.8 - 0.1j, (1, 2): 0.3j, (2, 2): -0.05}
    c_minus = {(0, 4): 0.2 + 0.1j}
    for label, rotation, y_sign in (("+Oz", np.eye(3), 1.0), ("-Oz", np.diag([1.0, -1.0, -1.0]), -1.0)):
        x_lab, y_lab = x_loc, y_sign * y_loc              # independent of screen_lab_xy
        r_lab, phi_lab = np.hypot(x_lab, y_lab), np.arctan2(y_lab, x_lab)
        E_p = sum(c * lg_mode(p, l, r_lab, phi_lab, s, w0, k) for (p, l), c in c_plus.items())
        E_m = sum(c * lg_mode(p, l, r_lab, phi_lab, s, w0, k) for (p, l), c in c_minus.items())
        Ex, Ey = (E_p + E_m) / math.sqrt(2.0), 1j * (E_p - E_m) / math.sqrt(2.0)

        x, y = screen_lab_xy(x_loc, y_loc, rotation)
        r, phi = np.hypot(x, y), np.arctan2(y, x)
        E_plus, E_minus = helicity_components(Ex, Ey)
        J = total_angular_momentum(1, 2, 1.0)
        for name, target, truth in (('plus', E_plus, c_plus), ('minus', E_minus, c_minus)):
            modes = [(p, l) for l in component_l_values(name, J) for p in range(4)]
            fit = fit_modes(target, build_basis(modes, r, phi, s, w0, k), A)
            error = max(abs(c - truth.get(mode, 0.0)) for mode, c in zip(modes, fit['coeffs']))
            check(f"{label} screen, {name}: max |c - c_true|", error, 1e-4)
            check(f"{label} screen, {name}: eps", fit['eps'], 1e-4)
        modes = [(p, l) for l in component_l_values('Ex', J) for p in range(3)]
        fit = fit_modes(Ex, build_basis(modes, r, phi, s, w0, k), A)
        check(f"{label} screen, Ex with l = J-1, J+1: eps", fit['eps'], 1e-4)
        modes = [(p, 2) for p in range(3)]
        fit = fit_modes(E_plus, build_basis(modes, r, phi, s, 1.5 * w0, k), A)
        print(f"  (info) {label} screen, plus at a wrong w0 = 1.5 w0_true, p <= 2: eps = {fit['eps']:.3e}")

    # A beam not in the basis (other waist, waist displaced by z0) on a wide, fine screen, where the modes are
    # orthonormal and the fitted c_0, c_1 are its projections.
    print("3. Series waist: one LG_0^l beam (w0 = 1.3 or 0.8 x the basis waist, waist at z0) fitted with p = 0..8:")
    s = 3.0 * z_R
    for l, w_ratio, z0_ratio in ((0, 1.3, 0.0), (0, 0.8, 0.3), (2, 1.3, -0.4)):
        w_true, z0 = w_ratio * w0, z0_ratio * z_R
        z_R_true = 0.5 * k * w_true ** 2
        w_s = max(w_true * math.sqrt(1.0 + ((s - z0) / z_R_true) ** 2), w0 * math.sqrt(1.0 + (s / z_R) ** 2))
        r, phi, A = circular_grid(8.0 * w_s, N_R=3000, N_phi=17)
        target = lg_mode(0, l, r, phi, s - z0, w_true, k)
        modes = [(p, l) for p in range(9)]
        fit = fit_modes(target, build_basis(modes, r, phi, s, w0, k), A)
        series = series_waist(fit['coeffs'], modes, w0, k)
        tag = f"l = {l}, w0 = {w_ratio:g} w_b, z0 = {z0_ratio:g} z_R"
        check(f"{tag}: w0_series/w0_true - 1", series['w0'] / w_true - 1.0, 1e-3)
        check(f"{tag}: (z0_series - z0)/z_R", (series['z0'] - z0) / z_R, 1e-3)
        check(f"{tag}: |t2 - t|", abs(series['t2'] - series['t']), 1e-3)

    # The (w0, z0) scan and _choose_w0 on a displaced beam: the true pair is on the grid and must be picked, with
    # all of sum|c|^2 in its own mode; at z0 = 0 no basis waist does that.
    print("4. Waist-position scan: LG_1^2 with w0 = 1.3 w_b at z0 = 0.5 z_R, grid containing the true (w0, z0):")
    l, w_true, z0_true = 2, 1.3 * w0, 0.5 * z_R
    r, phi, A = circular_grid(8.0 * w_true * math.sqrt(1.0 + ((s - z0_true) / (0.5 * k * w_true ** 2)) ** 2),
                              N_R=3000, N_phi=17)
    target = lg_mode(1, l, r, phi, s - z0_true, w_true, k)
    modes = [(p, l) for p in range(4)]
    w0_grid, z0_grid = w_true * np.array([0.8, 0.9, 1.0, 1.1, 1.25]), z0_true * np.array([-1.0, 0.0, 0.5, 1.0, 1.5])
    fits = [(wb, zb, fit_modes(target, build_basis(modes, r, phi, s - zb, wb, k), A)) for zb in z0_grid
            for wb in w0_grid]
    w_pick, z_pick, fit = _choose_w0(fits)
    check("picked w0/w0_true - 1", w_pick / w_true - 1.0, 1e-12)
    check("picked (z0 - z0_true)/z_R", (z_pick - z0_true) / z_R, 1e-12)
    check("dominant share - 1 at the picked point", _dominant_share(fit['coeffs']) - 1.0, 1e-6)
    best_at_origin = max(_dominant_share(f['coeffs']) for _, zb, f in fits if zb == 0.0)
    print(f"  (info) best dominant share with the basis waist at the origin: {best_at_origin:.4f}")

    # The Fresnel propagation behind coherent_source_prediction against the analytic modes (amplitude and phase,
    # Gouy and curvature included), and the N = 1, dl = 0 prediction, which is the laser mode itself.
    print("5. Coherent-source prediction: Fresnel-propagated LG modes vs lg_mode at s = 3 z_R:")
    s = 3.0 * z_R
    w_s = w0 * math.sqrt(1.0 + (s / z_R) ** 2)
    r, phi, A = circular_grid(5.0 * w_s, N_R=400, N_phi=33)
    rr = np.linspace(0.0, 8.0 * w0, 6000)
    for p, l in ((0, 0), (1, -2), (2, 3)):
        exact = lg_mode(p, l, r, phi, s, w0, k)
        rho, index = np.unique(r, return_inverse=True)
        fresnel = (fresnel_propagate_radial(lg_mode(p, l, rr, 0.0, 0.0, w0, k).real, rr, l, rho, s, k)[index]
                   * np.exp(1j * l * phi))
        check(f"(p,l) = ({p},{l}): max |E_Fresnel - u| / max |u|", np.abs(fresnel - exact).max() / np.abs(exact).max(),
              1e-6)
    pred, dl = coherent_source_prediction(r, phi, s, k, 1, -2, -2, 1, w0, None)
    check("N = 1 prediction of LG_1^-2: 1 - fidelity", 1.0 - prediction_fidelity(pred, lg_mode(1, -2, r, phi, s, w0, k),
                                                                                A), 1e-9)
    print("SELF-TEST " + ("PASSED" if ok else "FAILED"))
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('run_folder', nargs='?', help="run folder (default: the latest run)")
    parser.add_argument('--incident', action='store_true', help="fit incident_field.dat instead")
    parser.add_argument('--harmonic', type=int, help="only this harmonic N (default: every frequency)")
    parser.add_argument('--component', nargs='+', choices=_COMPONENTS, default=['plus', 'minus'])
    parser.add_argument('--l', type=int, help="override the azimuthal index (default: J - spin from the config)")
    parser.add_argument('--p-min', type=int, default=0)
    parser.add_argument('--p-max', type=int, default=3)
    parser.add_argument('--w0', help="comma-separated trial waists, e.g. 20,25,30 (default: a log scan around a "
                                     "second-moment estimate)")
    parser.add_argument('--w0-unit', default='lambda', help="unit of --w0: lambda (laser), w0 (laser waist), a.u.")
    parser.add_argument('--n-w0', type=int, default=_N_W0_DEFAULT, help="points in the default w0 scan")
    parser.add_argument('--scan-z0', action='store_true', help="also scan the basis waist position z0 (default: "
                        f"{_N_Z0_DEFAULT} values over +-{_Z0_SCAN_ZR:g} Rayleigh ranges; default off: z0 = 0)")
    parser.add_argument('--z0', help="comma-separated basis waist positions along propagation (implies a z0 scan)")
    parser.add_argument('--z0-unit', default='lambda', help="unit of --z0: lambda (laser), w0 (laser waist), a.u.")
    parser.add_argument('--n-z0', type=int, default=_N_Z0_DEFAULT, help="points in the default z0 scan (odd)")
    parser.add_argument('--self-test', action='store_true', help="run the synthetic checks and exit")
    args = parser.parse_args()

    if args.self_test:
        sys.exit(0 if self_test() else 1)
    if not 0 <= args.p_min <= args.p_max:
        parser.error("need 0 <= --p-min <= --p-max")

    field_filename = "incident_field.dat" if args.incident else "radiation_field.dat"
    if args.run_folder:
        field_filepath = os.path.join(args.run_folder, field_filename)
        if not os.path.exists(field_filepath):
            print(f"Error: '{field_filepath}' not found")
            sys.exit(1)
    else:
        field_filepath = find_latest_output_file(field_filename)
        print(f"No folder given; using latest run: {field_filepath}")
    run(field_filepath, args)


if __name__ == "__main__":
    main()
