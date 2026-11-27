#!/usr/bin/env python3
"""Generate the preliminary configs for the Sec. IV.B spectrum figure: 10 um and its 800 nm scaled twin, N = 1 and 2, for
the cases in CASES (main text: a0 = 0.1 with N_c = 30, a0 = 1 with N_c = 8; supplement: a0 = 0.5 with N_c = 30; a0 = 1
with N_c = 30 is kept as run, its 16 theta_* line was too short).

Base: config/reference/spectrum_line_lattice_2d.cfg (LG (0, 2), circular, polar lattice, dense spectrum on a radial
line, long_distance_direct_approx). Only the keys in OVERRIDES / per-run values change; the generated files are
complete configs and are not edited by hand (edit this script and regenerate).

Scaled choices (theory: manuscript Sec. III.B; run plan ~/Dropbox/pub/ms/2026-Thomson/runs_sec4_plan.md):
  Z_* = k_{1,ref} w0^2 / 2, theta_* = 2 / (k_{1,ref} w0), k_{1,ref} dressed with a0^2/2 (N = 1 scale for all N);
  Z/Z_* fixed per a0 and equal at both wavelengths (2.81 at a0 = 0.1, 2.85 at a0 = 1: screens L and H of the plan);
  line from 0 to theta_max (per case), so x_max = (theta_max/theta_*)(Z/Z_*) w0.
  Lattice spacing d from the aliasing rule theta_rep / theta_* = pi w0 / ((omega/omega_1)_max d) > theta_max/theta_*;
  for the N_c = 8 and a0 = 0.5 cases with the replica beyond 1.5 theta_max, so its skirt stays off the line.
  Line lengths from probe runs (2026-11-26, line to 80 theta_*): at a0 = 1, N_c = 8 and a0 = 0.5, N_c = 30 all but
  ~1e-3 of the flux is within 16 theta_*, so those lines go to 24 theta_*.

Usage: python3 make_configs.py   (writes configs/<id>.cfg next to this script and prints a summary)
"""
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
BASE = os.path.join(REPO, 'config', 'reference', 'spectrum_line_lattice_2d.cfg')
OUT = os.path.join(HERE, 'configs')
OUTPUT_ROOT = 'output/paper_2026/prelim/fig_spectrum'

W0_OVER_LAMBDA = 10.0
CUTOFF = 6
# laser_frequency (a.u.) and average_pz (mc): the run plan's common settings
LASERS = {'10um': (0.004556, 13.483), '800nm': (0.05695, 3.753)}   # gamma 13.520, 3.884: lambda_1 = 13.71 nm undressed at both
# per case (id tag): a0, N_c, Z/Z_*, theta_max/theta_*, number of line points, frequency step (omega/omega_1),
# and per N: (omega_min, omega_max) in omega/omega_1 and lattice spacing in lambda
CASES = {
    'a0_0p1': dict(a0=0.1, N_c=30.0, Z_over_Zs=2.81, theta_max=8.0, Nx=128, step=0.0005,
                   N={1: ((0.90, 1.10), 1.0), 2: ((1.90, 2.10), 1.0)}),
    'a0_1p0': dict(a0=1.0, N_c=30.0, Z_over_Zs=2.85, theta_max=16.0, Nx=192, step=0.001,
                   N={1: ((1.08, 1.60), 1.0), 2: ((2.16, 3.20), 0.5)}),
    'a0_1p0_Nc8': dict(a0=1.0, N_c=8.0, Z_over_Zs=2.85, theta_max=24.0, Nx=288, step=0.002,
                       N={1: ((1.08, 1.60), 0.5), 2: ((2.16, 3.20), 0.25)}),
    'a0_0p5_Nc30': dict(a0=0.5, N_c=30.0, Z_over_Zs=2.85, theta_max=24.0, Nx=288, step=0.001,
                        N={1: ((1.00, 1.18), 0.7), 2: ((2.00, 2.36), 0.35)}),
}
OVERRIDES = {
    'require_clean_build': 'true',
    'num_threads': '0',
    'laser_a0': None,  # per run
    'laser_wing_sigma_cutoff': str(CUTOFF),
    'laser_flat_duration': '0.0 cycles_adim',
    'laser_delay': '0.0 cycles_adim',
    'laser_lg_p': '0',
    'laser_lg_l': '2',
    'laser_lg_w0': f'{W0_OVER_LAMBDA:.1f} lambda',
    'detector_type': 'rectangular',
    'rectangular_detector_x_min': '0.0 w0',
    'rectangular_detector_y_min': '0.0 w0',
    'rectangular_detector_y_max': '0.0 w0',
    'rectangular_detector_Ny': '1',
    'dense_frequency_spectrum': 'true',
}


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
    raise KeyError(f'{key} not found in {BASE}')


def main():
    os.makedirs(OUT, exist_ok=True)
    with open(BASE) as f:
        base = f.readlines()
    print(f'{"id":<26}{"D":>8}{"Z [w0]":>10}{"x_max [w0]":>12}{"theta_max [urad]":>18}{"d [lambda]":>12}'
          f'{"omega/omega_1":>16}{"N_omega":>9}')
    for laser, (freq, pz) in LASERS.items():
        D, gamma = doppler(pz)
        for tag, p in CASES.items():
            a0 = p['a0']
            Zs_over_w0 = D * D / (1.0 + a0 * a0 / 2.0) * math.pi * W0_OVER_LAMBDA
            Z = p['Z_over_Zs'] * Zs_over_w0
            x_max = p['theta_max'] * p['Z_over_Zs']
            for N, ((wmin, wmax), d) in p['N'].items():
                run_id = f'{laser}_{tag}_N{N}'
                n_omega = int(round((wmax - wmin) / p['step'])) + 1
                lines = list(base)
                values = dict(OVERRIDES)
                values.update({
                    'output_folder': f'{OUTPUT_ROOT}/{run_id}',
                    'laser_a0': f'{a0}',
                    'laser_wing_sigma': f'{p["N_c"]:.1f} cycles_adim',
                    'laser_frequency': f'{freq}',
                    'average_pz': f'{pz} mc',
                    'rectangular_detector_distance': f'{Z:.1f} w0',
                    'rectangular_detector_x_max': f'{x_max:.3f} w0',
                    'rectangular_detector_Nx': str(p['Nx']),
                    'omega_min': f'{wmin} first_harmonic_frequency',
                    'omega_max': f'{wmax} first_harmonic_frequency',
                    'N_omega': str(n_omega),
                    'beam_lattice_spacing': f'{d} lambda',
                })
                for key, value in values.items():
                    set_key(lines, key, value)
                theta_max = math.atan(x_max / Z)
                header = (f'# GENERATED by reports/2026_thomson_paper/prelim/fig_spectrum/make_configs.py: do not edit.\n'
                          f'# {run_id}: Sec. IV.B spectrum figure (preliminary). LG (0, 2), a0 = {a0}, N_c = {p["N_c"]:g},'
                          f' cutoff {CUTOFF}, harmonic N = {N}.\n'
                          f'# D = {D:.4f} (gamma = {gamma:.4f}); Z/Z_* = {p["Z_over_Zs"]}, theta_max = '
                          f'{p["theta_max"]:g} theta_* = {theta_max * 1e6:.1f} urad.\n')
                with open(os.path.join(OUT, f'{run_id}.cfg'), 'w') as f:
                    f.write(header + ''.join(lines))
                print(f'{run_id:<26}{D:>8.4f}{Z:>10.1f}{x_max:>12.3f}{theta_max * 1e6:>18.1f}{d:>12}'
                      f'{f"{wmin}-{wmax}":>16}{n_omega:>9}')


if __name__ == '__main__':
    main()
