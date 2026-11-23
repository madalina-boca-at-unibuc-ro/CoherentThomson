"""Report Sec. 5: point-by-point comparison of runs A-D (65536 random electrons) with the same runs on the polar
lattice, and the lattice convergence check (run D, spacing halved).

For each pair and harmonic N, with F the total Faraday tensor (six components) on the screen:
  c = <F_lat, F_x> / ||F_lat||^2      the best complex scale; |c| ~ (electrons in x) / (electrons in the lattice run)
  eta = ||F_x - c F_lat||^2 / ||c F_lat||^2   the power of x not explained by the lattice field: for a random run
                                               its sampling speckle (incoherent part), for the fine lattice the
                                               lattice discretization error (squared)
The circular screen is equal-area, so every point has the same weight. Also printed per run: the orbital angular
momentum per photon with more digits than summarize_runs.py.

Run from anywhere: python3 compare_lattice_random.py
"""
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'future_dev' / 'tests' / 'bunch_length'))
from compare_bunch_length import read_field  # noqa: E402

REPORT = Path.home() / 'output/progress_report_20261029'
PAIRS = [('A', 'runs_A-D/20261029_132340', 'runs_A-D_lattice/20261031_093730'),
         ('B', 'runs_A-D/20261029_143739', 'runs_A-D_lattice/20261031_094539'),
         ('C', 'runs_A-D/20261029_140030', 'runs_A-D_lattice/20261031_094054'),
         ('D', 'runs_A-D/20261029_151457', 'runs_A-D_lattice/20261031_094906'),
         ('D fine lattice', 'runs_A-D_lattice/20261031_095357', 'runs_A-D_lattice/20261031_094906')]
L_THEORY = {('A', 1): 0, ('A', 2): -1, ('B', 1): 0, ('B', 2): -1, ('C', 1): -2, ('C', 2): -5, ('D', 1): -2, ('D', 2): -5}


def electron_count(run):
    return int(re.search(r'particle count\s+(\d+)', (run / 'run_log.txt').read_text()).group(1))


def orbital_per_photon(run):
    """omega_N Lambda_zz / P_z per harmonic, from plot_observables.py's section of run_log.txt."""
    text = (run / 'run_log.txt').read_text()
    omega_1 = float(re.search(r'fundamental_frequency\s+(\S+) a\.u\.', text).group(1))
    block = text.split('Screen-integrated observables')[1].split('Flux/density ratios')[0]
    rows = [list(map(float, l.split())) for l in block.splitlines() if re.match(r'\s+\d+\s+[\d.]+\s', l)]
    return {round(r[1]): r[1] * omega_1 * r[6] / r[9] for r in rows}


def main():
    print(f"{'pair':16} {'N':>2} {'|c|':>9} {'n_x/n_lat':>9} {'arg c':>9} {'eta':>9}")
    for label, x, lat in PAIRS:
        x, lat = REPORT / x, REPORT / lat
        om_x, F_x = read_field(x)
        om_l, F_l = read_field(lat)
        assert np.allclose(om_x, om_l)
        for i, N in enumerate(om_l):
            a, b = F_l[i].ravel(), F_x[i].ravel()
            c = np.vdot(a, b) / np.vdot(a, a).real
            eta = np.linalg.norm(b - c * a) ** 2 / (abs(c) ** 2 * np.vdot(a, a).real)
            print(f'{label:16} {N:2.0f} {abs(c):9.4f} {electron_count(x) / electron_count(lat):9.4f} '
                  f'{np.angle(c):9.2e} {eta:9.2e}')
    print('\norbital angular momentum per photon (hbar), random vs lattice, and expected N(m + sigma) - sigma')
    for label, x, lat in PAIRS[:4]:
        lx, ll = orbital_per_photon(REPORT / x), orbital_per_photon(REPORT / lat)
        for N in (1, 2):
            print(f'  {label} N = {N}: random {lx[N]:+.5f}  lattice {ll[N]:+.5f}  expected {L_THEORY[label, N]:+d}')


if __name__ == '__main__':
    main()
