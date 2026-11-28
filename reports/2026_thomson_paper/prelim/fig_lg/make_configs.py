#!/usr/bin/env python3
"""Generate the Sec. IV.D (LG content) and IV.E (ponderomotive lens) configs: LG (0, 2) and LG (2, 2), a0 = 0.1 and
0.25, 10 um, gamma = 13.52, N_c = 30, circular polarization, harmonics N = 1 and 2 (decided 2026-11-27; the same four
cases serve both sections, the paper / supplement split is chosen from the results).

Two stages, as in IV.B / IV.C (spectrum first):
  1. line: dense spectrum on a radial screen line (as ../fig_spectrum), one run per case and harmonic. At a0 = 0.25
     omega_{N,ref} lies below the line by ~1.5%, more than the 1.25% (N = 1) / 0.9% (N = 2) line width, so the maps are
     imaged at the peak of the screen-integrated spectrum (PEAKS below, from ../fig_maps/peak_frequencies.py-style
     analysis of these runs). The a0 = 0.1 lines are cheap and give the supplement spectra and the peak offsets.
  2. map: full circular screen. a0 = 0.1: one run per mode with N = 1, 2 at omega_{N,ref} (offset 0.25%, a fifth of
     the line width; agreed 2026-11-27), so it can run before the lines. a0 = 0.25: one run per harmonic at its peak,
     written only once PEAKS has the value.

Shared physics: config/reference/spectrum_line_lattice_2d.cfg (polar lattice, beam radius from the LG table: 2.37 w0
for (0, 2), 3.22 w0 for (2, 2); long_distance_direct_approx, trajectory_NT 100, cutoff 6). Each map config is its line
config (stage 1, harmonic 1) with only the screen, the frequencies and the lattice spacing changed.

Scaled choices (manuscript Sec. III.B; runs_sec4_plan.md screen L, for a0^2 N_c <~ 2):
  Z/Z_* = 2.81 at both a0; Z_* = k_{1,ref} w0^2 / 2 dressed with a0^2/2. theta_* Z_* = w0.
  Screen radius: weak-field Fresnel estimate (source u^N) puts 99.99% of the flux within 2.8-3.1 theta_* for (0, 2)
  and 3.7 / 4.8 theta_* (N = 1 / 2) for (2, 2); lines and maps go to 8 theta_* at a0 = 0.1 and the a0 = 0.25 lines to
  16 theta_* (the lens widens the emission; its map radius is set from the line runs).
  Lattice spacing: the largest of 1, 0.5, 0.25 lambda with the grating replica pi w0 / ((omega/omega_1)_max d) beyond
  REPLICA_MARGIN x theta_max.

Don't edit the configs by hand: change this script and regenerate.
Usage: python3 make_configs.py   (writes configs/<id>.cfg and prints a summary with cost estimates)
"""
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
BASE = os.path.join(REPO, 'config', 'reference', 'spectrum_line_lattice_2d.cfg')
OUT = os.path.join(HERE, 'configs')
OUTPUT_ROOT = 'output/paper_2026/prelim/fig_lg'

W0_OVER_LAMBDA = 10.0
CUTOFF = 6
N_C = 30.0
CYCLES = 2 * CUTOFF * N_C                 # trajectory span in cycles (360)
FREQ, PZ = 0.004556, 13.483               # 10 um, gamma = 13.520 (runs_sec4_plan.md common settings)
Z_OVER_ZS = 2.81
REPLICA_MARGIN = 1.5
SPACINGS = (1.0, 0.5, 0.25)
LG_RADIUS_W0 = {(0, 2): 2.37, (2, 2): 3.22}   # IoUtils::lg_beam_radius_in_w0 (beam_cylinder_radius -1)
MODES = [(0, 2), (2, 2)]
# a0 tag -> a0, line theta_max / theta_*, line points, frequency step, per N: line window (omega/omega_1)
AMPLITUDES = {
    '0p1': dict(a0=0.1, theta_max=8.0, Nx=128, step=0.0005, N={1: (0.90, 1.10), 2: (1.90, 2.10)}),
    '0p25': dict(a0=0.25, theta_max=16.0, Nx=192, step=0.0005, N={1: (0.94, 1.10), 2: (1.88, 2.20)}),
}
# map screen: circular, equal-area rings; N_PHI - 1 distinct azimuths (> 2 |J| + 1 = 13 for J = -6).
# a0 = 0.25 radius from the stage-1 lines (dev_notes 2026-11-28): LG (2, 2) has 0.4% (N = 1) and 1.9% (N = 2) of the
# flux beyond 8 theta_*, 1.5-2e-4 beyond 12; so 16 theta_*, with N_R doubled to keep the ring density.
MAP = {'0p1': dict(theta_max=8.0, N_R=512), '0p25': dict(theta_max=16.0, N_R=1024)}
N_PHI = 65
# (mode tag, a0 tag, N) -> omega_peak (omega/omega_1) of the screen-integrated spectrum of the stage-1 line run
# (parabolic refinement of the 0.0005 grid; dev_notes 2026-11-28). Stage 2: one map per entry, one harmonic per run.
PEAKS = {
    ('LG02', '0p1', 1): 1.0034, ('LG02', '0p1', 2): 2.0055,
    ('LG02', '0p25', 1): 1.0234, ('LG02', '0p25', 2): 2.0408,
    ('LG22', '0p1', 1): 1.0035, ('LG22', '0p1', 2): 2.0042,
    ('LG22', '0p25', 1): 1.0269, ('LG22', '0p25', 2): 2.0429,
}
# costs: dense line ~(6.4 + 0.5 N_omega) us, map ~2.5 us per electron, screen point and cycle (fig_spectrum / fig_maps)
OVERRIDES = {
    'require_clean_build': 'true',
    'num_threads': '0',
    'laser_wing_sigma': f'{N_C:.1f} cycles_adim',
    'laser_wing_sigma_cutoff': str(CUTOFF),
    'laser_flat_duration': '0.0 cycles_adim',
    'laser_delay': '0.0 cycles_adim',
    'laser_frequency': f'{FREQ}',
    'average_pz': f'{PZ} mc',
    'laser_lg_l': '2',
    'laser_lg_w0': f'{W0_OVER_LAMBDA:.1f} lambda',
    'detector_type': 'rectangular',
    'rectangular_detector_x_min': '0.0 w0',
    'rectangular_detector_y_min': '0.0 w0',
    'rectangular_detector_y_max': '0.0 w0',
    'rectangular_detector_Ny': '1',
    'dense_frequency_spectrum': 'true',
}
OPTIONAL_KEYS = {'omega_list'}   # absent from the base: appended


def doppler(pz):
    gamma = math.sqrt(1.0 + pz * pz)
    return gamma + pz, gamma


def set_key(lines, key, value):
    pattern = re.compile(rf'^(\s*){re.escape(key)}(\s+)(\S+(?:\s+[A-Za-z_]\w*)?)(\s*(#.*)?)$')
    for i, line in enumerate(lines):
        m = pattern.match(line.rstrip('\n'))
        if m:
            comment = f'  {m.group(5)}' if m.group(5) else ''
            lines[i] = f'{key:<23} {value}{comment}\n'
            return
    if key not in OPTIONAL_KEYS:
        raise KeyError(f'{key} not found in {BASE}')
    lines.append(f'{key:<23} {value}\n')


def spacing(omega_max, theta_max):
    return next(s for s in SPACINGS if math.pi * W0_OVER_LAMBDA / (omega_max * s) >= REPLICA_MARGIN * theta_max)


def sites(mode, d):
    return math.pi * (LG_RADIUS_W0[mode] * W0_OVER_LAMBDA / d) ** 2


def write(run_id, lines, header):
    with open(os.path.join(OUT, f'{run_id}.cfg'), 'w') as f:
        f.write(header + ''.join(lines))


def main():
    os.makedirs(OUT, exist_ok=True)
    with open(BASE) as f:
        base = f.readlines()
    D, gamma = doppler(PZ)
    total = 0.0
    print(f'{"id":<32}{"Z [w0]":>9}{"R [w0]":>8}{"d [lambda]":>11}{"electrons":>10}{"points":>8}{"N_omega":>8}'
          f'{"CPU h":>8}')
    for p, l in MODES:
        mode_tag = f'LG{p}{l}'
        for a_tag, amp in AMPLITUDES.items():
            a0 = amp['a0']
            Z = Z_OVER_ZS * D * D / (1.0 + a0 * a0 / 2.0) * math.pi * W0_OVER_LAMBDA   # in w0
            line_cfg = {}
            # stage 1: line spectra
            for N, (wmin, wmax) in amp['N'].items():
                run_id = f'10um_{mode_tag}_a0_{a_tag}_N{N}_line'
                x_max = amp['theta_max'] * Z_OVER_ZS
                d = spacing(wmax, amp['theta_max'])
                n_omega = int(round((wmax - wmin) / amp['step'])) + 1
                lines = list(base)
                values = dict(OVERRIDES)
                values.update({
                    'output_folder': f'{OUTPUT_ROOT}/{run_id}',
                    'laser_a0': f'{a0}',
                    'laser_lg_p': str(p),
                    'rectangular_detector_distance': f'{Z:.1f} w0',
                    'rectangular_detector_x_max': f'{x_max:.3f} w0',
                    'rectangular_detector_Nx': str(amp['Nx']),
                    'omega_min': f'{wmin} first_harmonic_frequency',
                    'omega_max': f'{wmax} first_harmonic_frequency',
                    'N_omega': str(n_omega),
                    'beam_lattice_spacing': f'{d} lambda',
                })
                for key, value in values.items():
                    set_key(lines, key, value)
                line_cfg[N] = lines
                header = (f'# GENERATED by reports/2026_thomson_paper/prelim/fig_lg/make_configs.py: do not edit.\n'
                          f'# {run_id}: Sec. IV.D/E line spectrum. LG ({p}, {l}), a0 = {a0}, N_c = {N_C:g}, cutoff '
                          f'{CUTOFF}, N = {N}; D = {D:.4f} (gamma = {gamma:.4f}); Z/Z_* = {Z_OVER_ZS}, line to '
                          f'{amp["theta_max"]:g} theta_*.\n')
                write(run_id, lines, header)
                n_e = sites((p, l), d)
                cpu = n_e * amp['Nx'] * CYCLES * (6.4 + 0.5 * n_omega) * 1e-6 / 3600
                total += cpu
                print(f'{run_id:<32}{Z:>9.1f}{x_max:>8.2f}{d:>11}{n_e:>10.0f}{amp["Nx"]:>8}{n_omega:>8}{cpu:>8.1f}')
            # stage 2: maps
            theta_max, n_r = MAP[a_tag]['theta_max'], MAP[a_tag]['N_R']
            jobs = []
            if a_tag == '0p1':   # stage 1: both harmonics at omega_N,ref in one run
                jobs.append((f'10um_{mode_tag}_a0_{a_tag}_map', 2.0, {'N_harmonics': '2', 'N_harmonics_min': '1'},
                             'N = 1, 2 at omega_N,ref'))
            for N in (1, 2):     # stage 2: one harmonic at its spectral peak
                peak = PEAKS.get((mode_tag, a_tag, N))
                if peak is None:
                    print(f'{f"10um_{mode_tag}_a0_{a_tag}_N{N}_map":<32}  (after the line runs: PEAKS)')
                    continue
                jobs.append((f'10um_{mode_tag}_a0_{a_tag}_N{N}_map', peak,
                             {'omega_list': f'{N}:{peak} first_harmonic_frequency'}, f'N = {N} at the peak {peak}'))
            for run_id, omega_max, freq_keys, what in jobs:
                R = theta_max * Z_OVER_ZS
                d = spacing(omega_max, theta_max)
                lines = list(line_cfg[1])
                values = {
                    'output_folder': f'{OUTPUT_ROOT}/{run_id}',
                    'dense_frequency_spectrum': 'false',
                    'detector_type': 'circular',
                    'circular_detector_R_min': '0.0 lambda',
                    'circular_detector_R_max': f'{R:.3f} w0',
                    'circular_detector_N_R': str(n_r),
                    'circular_detector_N_phi': str(N_PHI),
                    'circular_detector_distance': f'{Z:.1f} w0',
                    'beam_lattice_spacing': f'{d} lambda',
                }
                values.update(freq_keys)
                for key, value in values.items():
                    set_key(lines, key, value)
                header = (f'# GENERATED by reports/2026_thomson_paper/prelim/fig_lg/make_configs.py: do not edit.\n'
                          f'# {run_id}: Sec. IV.D/E field map. LG ({p}, {l}), a0 = {a0}, N_c = {N_C:g}, {what}; '
                          f'circular screen Z/Z_* = {Z_OVER_ZS}, R = {theta_max:g} theta_*, {n_r} x {N_PHI}; '
                          f'd = {d} lambda.\n')
                write(run_id, lines, header)
                n_e = sites((p, l), d)
                n_freq = 2 if 'N_harmonics' in freq_keys else 1
                cpu = n_e * n_r * N_PHI * CYCLES * 2.5e-6 * n_freq / 3600
                total += cpu
                print(f'{run_id:<32}{Z:>9.1f}{R:>8.2f}{d:>11}{n_e:>10.0f}{n_r * N_PHI:>8}{n_freq:>8}{cpu:>8.1f}')
    print(f'total {total:.0f} CPU h (estimate)')


if __name__ == '__main__':
    main()
