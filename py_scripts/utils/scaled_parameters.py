#!/usr/bin/env python3
"""
Scaled parameters of a finished run (theory/specs/11_scaled_parameters.md, with its §10 amendments): every scaled
quantity the paper uses, computed from the run folder's run_log.txt (and config.cfg for what the log lacks). Plots and
figure scripts take their scaled units (theta_*, Z_*, ...) from here and nowhere else.

Leaf module: imports only utils.run_output_utils (never plot_field or plot_observables), so any script can use it.

    sp = scaled_parameters(run_dir)     # ScaledParameters; per-harmonic rows in sp.harmonic_rows
    sp.warnings                         # the §5 checks
    to_scaled_angles(x, y, sp)          # screen coordinates (lambda) -> (theta_x/theta_*, theta_y/theta_*)
    dominant_component(sp)              # 'plus' or 'minus'

Command line:
    python3 py_scripts/utils/scaled_parameters.py [run_folder ...] [--write-log] [--expect L|H]
                                                  [--table md|csv] [--output FILE]
One run: prints every block. Several runs (or --table): one row per run and harmonic, Markdown or CSV.
"""
import dataclasses
import math
import os
import re
import sys

import numpy as np
import scipy.special

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from utils.run_output_utils import find_run_log_sections, write_run_log_section  # noqa: E402

HARTREE_EV = 27.211386245988
ELECTRON_REST_MEV = 0.51099895
RUN_LOG_HEADER = "Scaled parameters (scaled_parameters.py)"
SCREEN_GROUPS = {'L': 2.81, 'H': 2.85}          # Z/Z_* of the run plan's screens (amendment 4: Z/Z_* only)
_NUM = r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?'


# ------------------------------------------------------------------------------------------------------------------
# The LG mode's peak intensity (one copy; plot_spectrum_map.peak_a2_fraction wraps it)
# ------------------------------------------------------------------------------------------------------------------

def lg_peak_intensity(p, l):
    """kappa = max over rho of |u|^2 for the C++ LaguerreGaussLaser amplitude with its custom Npn: 2 on axis for (0, 0),
    0.541 for (0, 2). |a|^2 peaks at kappa a0^2 (cycle average kappa a0^2/2)."""
    n = abs(int(l))
    p = int(p)
    Npn = np.sqrt(2.0) / math.factorial(n) * np.sqrt(math.factorial(p + n) / math.factorial(p))
    u = np.linspace(0.0, 10.0 + 4.0 * (p + n), 200001)  # rho^2 / w0^2
    amplitude = Npn * (2.0 * u) ** (n / 2.0) * scipy.special.hyp1f1(-p, n + 1, 2.0 * u) * np.exp(-u)
    return float(np.max(amplitude ** 2))


# ------------------------------------------------------------------------------------------------------------------
# run_log.txt parsing
# ------------------------------------------------------------------------------------------------------------------

def _sections(content):
    """{title: text} of run_log.txt's top-level sections."""
    marks = find_run_log_sections(content)
    return {title: content[start:(marks[i + 1][0] if i + 1 < len(marks) else len(content))]
            for i, (start, title) in enumerate(marks)}


def _line(text, label):
    """The value part of the first top-level ('  label  ...', two-space indent) line with this label, or None."""
    m = re.search(rf'(?m)^  {re.escape(label)}\s+(.*)$', text or '')
    return m.group(1).strip() if m else None


def _first_number(text, label, unit=None):
    """First number on the line; with `unit`, the number just before that unit word (e.g. 'lambda', 'a.u.')."""
    value = _line(text, label)
    if value is None:
        return None
    if unit is not None:
        m = re.search(rf'({_NUM})\s+{re.escape(unit)}(?:\s|$)', value)
        return float(m.group(1)) if m else None
    m = re.search(_NUM, value)
    return float(m.group(0)) if m else None


def _vector(text, label):
    value = _line(text, label)
    m = re.search(r'\(([^)]*)\)', value or '')
    return [float(x) for x in m.group(1).split(',')] if m else None


# ------------------------------------------------------------------------------------------------------------------
# The quantities
# ------------------------------------------------------------------------------------------------------------------

@dataclasses.dataclass
class ScaledParameters:
    run_dir: str
    warnings: list = dataclasses.field(default_factory=list)
    notes: list = dataclasses.field(default_factory=list)
    head_on: bool = False
    # physics parameters
    a0: float = None
    N_c: float = None
    a0sq_Nc: float = None
    p: int = None
    l: int = None
    sigma: float = None
    m: int = None
    aligned: bool = None
    w0_over_lambda: float = None
    R_over_w0: float = None
    d_over_w0: float = None
    n_electrons: int = None
    D: float = None
    gamma: float = None
    beta: float = None
    A_ref: float = None
    kappa: float = None
    # screen and frequency scales
    k1_over_kL: float = None
    lambda_1: float = None              # in lambda_L
    Z_star: float = None                # in lambda_L
    theta_star: float = None            # rad
    Z_over_Zstar: float = None
    theta_max: float = None             # rad
    theta_max_over_thetastar: float = None
    screen: str = None                  # 'circular', 'rectangular', 'line', 'spherical'
    frequency_mode: str = None
    harmonics: list = None
    harmonic_rows: list = None
    # validity diagnostics
    epsilon: float = None
    xi: float = None
    fresnel_phase_error: float = None
    paraxial_angle: float = None
    # physical block
    lambda_L_nm: float = None
    lambda_1_nm: float = None
    photon_energy_eV: float = None
    w0_um: float = None
    R_um: float = None
    sigma_field_fs: float = None
    fwhm_intensity_fs: float = None
    t_int_fs: float = None
    kinetic_energy_MeV: float = None


def _warn(sp, message):
    sp.warnings.append(message)


def _frequencies(section):
    """[(omega/omega_1, N or None)] from the frequencies_list lines of a harmonic or omega_list run, or
    ('dense', first, last) for a dense spectrum."""
    rows = re.findall(rf'(?m)^\s+\[\d+\]\s+\S+\s+\(({_NUM}) omega_1\)(?:\s+N = (\d+))?\s*$', section or '')
    if rows:
        return [(float(x), int(n) if n else None) for x, n in rows]
    first = re.search(rf'(?m)^\s+first = \S+\s+\(({_NUM}) omega_1\)', section or '')
    last = re.search(rf'(?m)^\s+last\s+= \S+\s+\(({_NUM}) omega_1\)', section or '')
    if first and last:
        return ('dense', float(first.group(1)), float(last.group(1)))
    return None


def scaled_parameters(run_dir):
    """Every quantity of spec 11 §4 (+ §10) for one finished run; None and a warning for what is missing."""
    sp = ScaledParameters(run_dir=os.path.abspath(run_dir))
    log_path = os.path.join(run_dir, 'run_log.txt')
    config_path = os.path.join(run_dir, 'config.cfg')
    with open(log_path) as f:
        content = f.read()
    sec = _sections(content)
    laser, beam, det, freq = (sec.get(k) for k in ('Laser', 'Electron beam', 'Detector', 'Frequency spectrum'))
    for name, text in (('Laser', laser), ('Electron beam', beam), ('Detector', det), ('Frequency spectrum', freq)):
        if text is None:
            _warn(sp, f"run_log.txt has no '{name}' section")

    # ---- laser ----
    n_L = _vector(laser, 'n')
    sp.a0 = _first_number(laser, 'a0')
    zeta = []
    for i in (1, 2):
        v = _vector(laser, f'zeta_{i}')
        zeta.append(complex(v[0], v[1]) if v else None)
    lp = _line(laser, 'LG p, l')
    if lp:
        sp.p, sp.l = (int(x) for x in lp.split(','))
    w0 = _first_number(laser, 'LG w0', 'lambda')            # in lambda_L
    sp.w0_over_lambda = w0
    lambda_au = _first_number(laser, 'lambda', 'a.u.')
    sp.lambda_L_nm = _first_number(laser, 'lambda (nm)')
    period_fs = _first_number(laser, 'period T (fs)')
    flat_cycles = _first_number(laser, 'flat_duration (cycles)')
    wing_phase = _first_number(laser, 'wing_sigma (phase)')
    if flat_cycles not in (None, 0.0):
        _warn(sp, f'envelope: flat_duration = {flat_cycles} cycles, not a pure Gaussian: N_c = None')
    elif wing_phase is not None:
        sp.N_c = wing_phase / (2.0 * math.pi)

    # ---- beam ----
    p_vec = _vector(beam, 'average_p (px,py,pz)')
    sp.n_electrons = int(_first_number(beam, 'particle count')) if _line(beam, 'particle count') else None
    R_sheet = _first_number(beam, 'cylinder radius', 'lambda')
    d = _first_number(beam, 'lattice spacing', 'lambda')
    if d is None:
        sp.notes.append('no lattice spacing in the log (random beam or older run): no replica check')

    # ---- geometry check (head-on, screen normal +z) ----
    det_dir = _vector(det, 'direction (theta,phi)')
    sp.head_on = bool(n_L and p_vec and det_dir is not None
                      and n_L[2] < -1 + 1e-9 and abs(p_vec[0]) < 1e-12 and abs(p_vec[1]) < 1e-12 and p_vec[2] >= 0
                      and abs(det_dir[0]) < 1e-12)
    if not sp.head_on:
        _warn(sp, 'geometry is not head-on (laser along -z, beam along +z, screen normal +z): no scaled block')

    # ---- helicity and AM about lab +z ----
    if n_L and all(z is not None for z in zeta):
        s_nz = math.copysign(1.0, n_L[2])
        norm = abs(zeta[0]) ** 2 + abs(zeta[1]) ** 2
        sp.sigma = 2.0 * (np.conj(zeta[0]) * zeta[1]).imag / norm * s_nz
        if sp.l is not None:
            sp.m = int(round(sp.l * s_nz))
            sp.aligned = None if sp.m == 0 or abs(sp.sigma) < 1e-9 else (math.copysign(1, sp.m) == math.copysign(1, sp.sigma))

    # ---- Doppler factor, dressing ----
    if p_vec and n_L:
        p0 = math.sqrt(1.0 + sum(x * x for x in p_vec))
        sp.D = p0 - sum(a * b for a, b in zip(n_L, p_vec))
        sp.gamma = p0
        sp.beta = math.sqrt(sum(x * x for x in p_vec)) / p0
        sp.kinetic_energy_MeV = (p0 - 1.0) * ELECTRON_REST_MEV
    if sp.a0 is not None:
        sp.A_ref = sp.a0 ** 2 / 2.0
        if sp.N_c is not None:
            sp.a0sq_Nc = sp.a0 ** 2 * sp.N_c
    if sp.p is not None:
        sp.kappa = lg_peak_intensity(sp.p, sp.l)
    else:
        laser_type = _line(laser, 'type')
        sp.kappa = 1.0 if laser_type and 'plane' in laser_type else None

    # ---- frequency scales and the check against the log ----
    if sp.D is not None and sp.A_ref is not None:
        sp.k1_over_kL = sp.D ** 2 / (1.0 + sp.A_ref)
        sp.lambda_1 = 1.0 / sp.k1_over_kL
        if sp.lambda_L_nm is not None:
            sp.lambda_1_nm = sp.lambda_L_nm * sp.lambda_1
        log_ratio = _first_number(freq, 'fundamental / laser omega')
        if log_ratio is not None and abs(sp.k1_over_kL / log_ratio - 1.0) > 1e-6:
            _warn(sp, f'frequency: D^2/(1 + a0^2/2) = {sp.k1_over_kL:.7g} but the log has fundamental / laser omega '
                      f'= {log_ratio:.7g} (A_ref or geometry mismatch)')
        omega1_au = _first_number(freq, 'fundamental_frequency', 'a.u.')
        if omega1_au is not None:
            sp.photon_energy_eV = omega1_au * HARTREE_EV
    if sp.k1_over_kL is not None and w0 is not None:
        k1w0 = 2.0 * math.pi * sp.k1_over_kL * w0
        sp.Z_star = 0.5 * k1w0 * w0                 # k_1 w0^2 / 2, in lambda_L
        sp.theta_star = 2.0 / k1w0
    if sp.D is not None and w0 is not None:
        sp.epsilon = 2.0 / (sp.D * 2.0 * math.pi * w0)
    if sp.beta is not None and sp.N_c is not None and w0 is not None:
        z_R = math.pi * w0 ** 2                     # k_L w0^2 / 2 in lambda_L
        sp.xi = sp.beta / (1.0 + sp.beta) * sp.N_c / z_R

    # ---- screen (amendment 1: lines) ----
    Z = _first_number(det, 'distance', 'lambda')
    det_type = _line(det, 'type') or ''
    grid = re.match(r'(\d+) x (\d+)', _line(det, 'grid') or '')
    n1, n2 = (int(grid.group(1)), int(grid.group(2))) if grid else (None, None)
    rho_max = None
    if 'Circular' in det_type:
        rho_max = _first_number(det, 'R_max', 'lambda')
        sp.screen = 'line' if n2 == 1 else 'circular'
    elif 'Rectangular' in det_type:
        xs = [_first_number(det, k, 'lambda') for k in ('x_min', 'x_max', 'y_min', 'y_max')]
        if None not in xs:
            x_min, x_max, y_min, y_max = xs
            if n1 == 1 or n2 == 1:
                sp.screen = 'line'
                rho_max = max(math.hypot(x, y) for x in (x_min, x_max) for y in (y_min, y_max))
            else:
                sp.screen = 'rectangular'
                rho_max = min(x_max, -x_min, y_max, -y_min)
    else:
        sp.screen = 'spherical'
        sp.notes.append('spherical screen: no screen quantities')
    if Z is not None and rho_max is not None:
        sp.theta_max = math.atan(rho_max / Z)
        sp.paraxial_angle = sp.theta_max
        if sp.Z_star:
            sp.Z_over_Zstar = Z / sp.Z_star
            sp.theta_max_over_thetastar = sp.theta_max / sp.theta_star
    if R_sheet is not None and w0:
        sp.R_over_w0 = R_sheet / w0
    if d is not None and w0:
        sp.d_over_w0 = d / w0

    # ---- harmonics (amendment 2) ----
    freqs = _frequencies(freq)
    band = None
    if sp.kappa is not None and sp.A_ref is not None:
        band = ((1.0 + sp.A_ref) / (1.0 + sp.kappa * sp.A_ref), 1.0 + sp.A_ref)
    omega_max_N = {}
    if isinstance(freqs, tuple):
        sp.frequency_mode = 'dense'
        _, first, last = freqs
        if band is None:
            Ns = range(max(1, int(round(first))), int(round(last)) + 1)
            omega_max_N = {N: min(last, N + 0.5) for N in Ns}
        else:
            low, high = band
            for N in range(1, int(last / low) + 2):
                if N * low <= last and N * high >= first:
                    omega_max_N[N] = min(last, 0.5 * (N * high + (N + 1) * low))
    elif freqs:
        sp.frequency_mode = 'omega_list' if any(abs(x - round(x)) > 1e-9 for x, _ in freqs) else 'harmonics'
        for x, N in freqs:
            N = N if N is not None else int(round(x))
            omega_max_N[N] = max(omega_max_N.get(N, 0.0), x)
    else:
        _warn(sp, 'no frequencies_list in the log')
    sp.harmonics = sorted(omega_max_N)

    rows = []
    for N in sp.harmonics:
        w_max = omega_max_N[N]
        row = dict(N=N, omega_max=w_max)
        if sp.m is not None and sp.sigma is not None and abs(abs(sp.sigma) - 1) < 1e-9:
            row['J_N'] = N * (sp.m + int(round(sp.sigma)))
            row['L_N'] = row['J_N'] - int(round(sp.sigma))
            row['component'] = 'plus' if sp.sigma > 0 else 'minus'
        if sp.N_c:
            row['linewidth'] = 1.0 / (N * sp.N_c)
        if band is not None:
            row['band_low'], row['band_high'] = N * band[0], N * band[1]
            if sp.N_c and sp.a0 is not None:
                row['N_fringes'] = math.sqrt(math.pi / 2.0) * N * sp.N_c * sp.kappa * sp.a0 ** 2 / 2.0
        if d is not None and w0 and sp.theta_max_over_thetastar:
            row['replica_over_thetastar'] = math.pi * w0 / (w_max * d)
            row['aliasing_margin'] = row['replica_over_thetastar'] / sp.theta_max_over_thetastar
            if row['aliasing_margin'] <= 1.25:
                _warn(sp, f'aliasing: N = {N}: lattice replica at {row["replica_over_thetastar"]:.3g} theta_*, '
                          f'margin {row["aliasing_margin"]:.3g} <= 1.25 over the screen')
        if sp.theta_max_over_thetastar and w0 and sp.k1_over_kL:
            row['angular_detuning'] = sp.theta_max_over_thetastar ** 2 / (sp.k1_over_kL * (2 * math.pi * w0) ** 2)
            if 'linewidth' in row:
                row['angular_detuning_over_linewidth'] = row['angular_detuning'] / row['linewidth']
        rows.append(row)
    sp.harmonic_rows = rows

    if sp.harmonics and sp.k1_over_kL and Z and rho_max is not None:
        k_max = 2.0 * math.pi * sp.k1_over_kL * max(omega_max_N.values())
        sp.fresnel_phase_error = k_max * (rho_max + (R_sheet or 0.0)) ** 4 / (8.0 * Z ** 3)

    # ---- physical block ----
    if w0 is not None and lambda_au is not None and sp.lambda_L_nm is not None:
        sp.w0_um = w0 * sp.lambda_L_nm * 1e-3
        if R_sheet is not None:
            sp.R_um = R_sheet * sp.lambda_L_nm * 1e-3
    if sp.N_c and period_fs:
        sp.sigma_field_fs = sp.N_c * period_fs
        sp.fwhm_intensity_fs = math.sqrt(2.0 * math.log(2.0)) * sp.N_c * period_fs
        sp.t_int_fs = 0.59 * sp.N_c * period_fs

    # ---- remaining checks (§5) ----
    if sp.theta_max_over_thetastar is not None and sp.theta_max_over_thetastar < 5:
        sp.notes.append(f'containment: theta_max = {sp.theta_max_over_thetastar:.3g} theta_* < 5: the screen may clip')
    if sp.paraxial_angle is not None and sp.paraxial_angle >= 0.1:
        sp.notes.append(f'paraxial: theta_max = {sp.paraxial_angle:.3g} rad >= 0.1')
    if sp.fresnel_phase_error is not None and sp.fresnel_phase_error >= 0.1:
        sp.notes.append(f'Fresnel: next phase term {sp.fresnel_phase_error:.3g} rad >= 0.1 over the screen')
    if not sp.head_on:
        for name in ('Z_star', 'theta_star', 'Z_over_Zstar', 'theta_max_over_thetastar'):
            setattr(sp, name, None)
    return sp


def check_screen_group(sp, group):
    """--expect L|H (amendment 4): Z/Z_* within 1% of the group's value."""
    target = SCREEN_GROUPS[group]
    if sp.Z_over_Zstar is None or abs(sp.Z_over_Zstar / target - 1.0) > 0.01:
        _warn(sp, f'screen group {group}: Z/Z_* = {sp.Z_over_Zstar} is not within 1% of {target}')


def to_scaled_angles(x, y, sp, Z=None):
    """Screen coordinates (in lambda_L, as plot_field.get_screen_coordinates returns them) -> (theta_x/theta_*,
    theta_y/theta_*), with theta = x/Z as plot_field's angular axes. Z defaults to the run's screen distance."""
    if Z is None:
        with open(os.path.join(sp.run_dir, 'run_log.txt')) as f:
            Z = _first_number(_sections(f.read()).get('Detector'), 'distance', 'lambda')
    return np.asarray(x) / Z / sp.theta_star, np.asarray(y) / Z / sp.theta_star


def dominant_component(sp):
    """'plus' or 'minus': E_sigma, sigma the laser spin about lab +z."""
    if sp.sigma is None or abs(sp.sigma) < 1e-9:
        raise ValueError(f'laser spin sigma = {sp.sigma}: no dominant helicity')
    return 'plus' if sp.sigma > 0 else 'minus'


# ------------------------------------------------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------------------------------------------------

_BLOCKS = [
    ('Physics parameters', ['a0', 'N_c', 'a0sq_Nc', 'p', 'l', 'sigma', 'm', 'aligned', 'w0_over_lambda', 'R_over_w0',
                            'd_over_w0', 'n_electrons', 'D', 'gamma', 'beta', 'A_ref', 'kappa']),
    ('Screen and frequency scales', ['k1_over_kL', 'lambda_1', 'Z_star', 'theta_star', 'Z_over_Zstar', 'theta_max',
                                     'theta_max_over_thetastar', 'screen', 'frequency_mode', 'harmonics']),
    ('Validity diagnostics', ['epsilon', 'xi', 'fresnel_phase_error', 'paraxial_angle']),
    ('Physical block', ['lambda_L_nm', 'lambda_1_nm', 'photon_energy_eV', 'w0_um', 'R_um', 'sigma_field_fs',
                        'fwhm_intensity_fs', 't_int_fs', 'kinetic_energy_MeV']),
]
_UNITS = {'lambda_1': 'lambda_L', 'Z_star': 'lambda_L', 'theta_star': 'rad', 'theta_max': 'rad',
          'paraxial_angle': 'rad', 'fresnel_phase_error': 'rad'}
_ROW_KEYS = ['N', 'omega_max', 'J_N', 'L_N', 'component', 'linewidth', 'band_low', 'band_high', 'N_fringes',
             'replica_over_thetastar', 'aliasing_margin', 'angular_detuning', 'angular_detuning_over_linewidth']


def _fmt(v):
    if v is None:
        return '-'
    if isinstance(v, bool) or isinstance(v, (int, np.integer)) or isinstance(v, str) or isinstance(v, list):
        return str(v)
    return f'{v:.6g}'


def report_lines(sp):
    lines = [f'  run                         {sp.run_dir}']
    for title, keys in _BLOCKS:
        lines.append(f'  {title}:')
        for k in keys:
            unit = f'  {_UNITS[k]}' if k in _UNITS else ''
            lines.append(f'    {k:<34}{_fmt(getattr(sp, k))}{unit}')
    lines.append('  Per harmonic:')
    lines.append('    ' + ' '.join(f'{k:>12}' for k in _ROW_KEYS))
    for row in sp.harmonic_rows or []:
        lines.append('    ' + ' '.join(f'{_fmt(row.get(k)):>12}' for k in _ROW_KEYS))
    lines.append('  Checks:')
    for w in sp.warnings:
        lines.append(f'    WARNING: {w}')
    for n in sp.notes:
        lines.append(f'    note: {n}')
    if not sp.warnings and not sp.notes:
        lines.append('    all passed')
    return lines


_TABLE_KEYS = ['run', 'a0', 'N_c', 'a0sq_Nc', 'p', 'l', 'sigma', 'D', 'w0_over_lambda', 'R_over_w0', 'd_over_w0',
               'n_electrons', 'Z_over_Zstar', 'theta_max_over_thetastar', 'lambda_L_nm', 'lambda_1_nm']


def table(sps, fmt='md'):
    keys = _TABLE_KEYS + _ROW_KEYS
    rows = []
    for sp in sps:
        base = {k: getattr(sp, k, None) for k in _TABLE_KEYS if k != 'run'}
        base['run'] = os.path.basename(sp.run_dir.rstrip('/'))
        for row in sp.harmonic_rows or [{}]:
            rows.append({**base, **row})
    if fmt == 'csv':
        return [','.join(keys)] + [','.join(_fmt(r.get(k)) for k in keys) for r in rows]
    return (['| ' + ' | '.join(keys) + ' |', '|' + '---|' * len(keys)]
            + ['| ' + ' | '.join(_fmt(r.get(k)) for k in keys) + ' |' for r in rows])


def main(argv):
    args = list(argv)
    opts = {}
    for flag in ('--expect', '--table', '--output'):
        if flag in args:
            i = args.index(flag)
            opts[flag] = args[i + 1]
            del args[i:i + 2]
    write_log = '--write-log' in args
    if write_log:
        args.remove('--write-log')
    if not args:
        print(__doc__)
        return 1
    sps = []
    for run in args:
        sp = scaled_parameters(run)
        if '--expect' in opts:
            check_screen_group(sp, opts['--expect'])
        sps.append(sp)
        if write_log:
            write_run_log_section(run, RUN_LOG_HEADER, report_lines(sp))
    if len(sps) == 1 and '--table' not in opts:
        out = report_lines(sps[0])
    else:
        out = table(sps, opts.get('--table', 'md'))
    text = '\n'.join(out) + '\n'
    if '--output' in opts:
        with open(opts['--output'], 'w') as f:
            f.write(text)
        print(f"Wrote {opts['--output']}")
    else:
        print(text, end='')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
