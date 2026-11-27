#!/usr/bin/env python3
"""Compare the 10 um and 800 nm field maps (prelim/fig_maps) point by point in scaled units.

Both screens are the same grid in units of w0 (circular, equal-area, Z/Z_* = 2.85), so screen point i is the same
(theta/theta_*, phi) in both runs. For each harmonic: E_- = (E_x + i E_y)/sqrt(2) (dominant helicity), times
exp(-i k (R - Z)) with each run's own k (spherical carrier removed), normalized, global phase fitted out. Reports
the mismatch 1 - |<a,b>| / (|a||b|) and the relative L2 difference, also inside theta < 4, 8 theta_*.
Usage: python3 maps_twin_compare.py   (reads ~/output/paper_2026/prelim/fig_maps/)
"""
import glob
import os
import re

import numpy as np

RUNS = os.path.expanduser('~/output/paper_2026/prelim/fig_maps')
C = 137.035999084


def latest(run_id):
    return sorted(glob.glob(os.path.join(RUNS, run_id, '*/')))[-1]


def load(run):
    cols = open(run + 'radiation_field.dat').read().split('\n', 4)[3].split()
    names = ['LR_F10', 'LR_F20', 'SR_F10', 'SR_F20', 'BR_F10', 'BR_F20']
    use = [cols.index(n + s) for n in names for s in ('_re', '_im')]
    d = np.loadtxt(run + 'radiation_field.dat', comments='#', skiprows=4, usecols=[1] + use)
    omega = d[0, 0]
    v = d[:, 1:]
    F = lambda k: v[:, 2 * k] + 1j * v[:, 2 * k + 1]
    Ex = C * (F(0) + F(2) + F(4))
    Ey = C * (F(1) + F(3) + F(5))
    log = open(run + 'run_log.txt').read()
    cfg = open(run + 'config.cfg').read()
    g = lambda k: float(re.search(rf'^{k}\s+(\S+)', cfg, re.M).group(1))
    NR, NP = int(g('circular_detector_N_R')), int(g('circular_detector_N_phi'))
    w0 = float(re.search(r'LG w0\s+([0-9.e+-]+) a\.u\.', log).group(1))
    k1 = float(re.search(r'fundamental / laser omega\s+([0-9.e+-]+)', log).group(1)) * \
        float(re.search(r'^laser_frequency\s+(\S+)', cfg, re.M).group(1)) / C
    Z = g('circular_detector_distance') * w0
    Rm = g('circular_detector_R_max') * w0
    r = np.sqrt(np.arange(NR) / (NR - 1)) * Rm
    rho = np.repeat(r, NP)
    k = omega * k1
    Em = (Ex + 1j * Ey) / np.sqrt(2) * np.exp(-1j * k * (np.sqrt(rho ** 2 + Z ** 2) - Z))
    theta_star = 2.0 / (k1 * w0)
    t = np.arctan(rho / Z) / theta_star
    keep = np.tile(np.arange(NP) < NP - 1, NR)          # drop the repeated phi = 2 pi column (equal-area: equal weights)
    return omega, Em[keep], t[keep]


def compare(a, b):
    ov = np.vdot(a, b)
    mis = 1 - abs(ov) / (np.linalg.norm(a) * np.linalg.norm(b))
    ph = np.exp(-1j * np.angle(ov))
    diff = np.linalg.norm(a / np.linalg.norm(a) - ph * b / np.linalg.norm(b))
    return mis, diff


for N in (1, 2):
    wa, a, ta = load(latest(f'10um_a0_1p0_Nc8_N{N}'))
    wb, b, tb = load(latest(f'800nm_a0_1p0_Nc8_N{N}'))
    assert np.allclose(ta, tb, rtol=1e-3), 'scaled grids differ'
    out = []
    for cut in (4, 8, 1e9):
        m = ta <= cut
        mis, diff = compare(a[m], b[m])
        out.append(f'theta<={cut if cut < 1e8 else "all":>3}: 1-overlap {mis:.1e}, rel. L2 {diff:.1e}')
    print(f'N = {N} (omega {wa:.4f} / {wb:.4f}):  ' + ';  '.join(out))
