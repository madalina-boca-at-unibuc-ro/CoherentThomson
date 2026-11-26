#!/usr/bin/env python3
"""
Sec. IV.B figure (spectrum): one PDF panel per case, assembled into one figure in the manuscript
(TeX/sections/04-flat_sheet.tex). Each panel: left, the (theta, omega) map of |E_sigma|^2 (the dominant helicity) of the
10 um run versus theta/theta_* (utils/scaled_parameters.py), log scale, normalized to its maximum; right, sharing the omega axis, the screen-integrated photon spectrum
dN/d(omega/omega_1) of the 10 um run and of its 800 nm scaled twin, each normalized to its own peak (their photon
numbers differ by ~(D_10um/D_800nm)^(2N) at the same number of electrons; quoted in the text, not plotted).
A separate PDF holds the shared colorbar.

Data: the line-spectrum runs of reports/2026_thomson_paper/prelim/fig_spectrum/ (latest run of each id under
~/output/paper_2026/prelim/fig_spectrum/<id>/).

Usage: python3 fig_spectrum.py [case ...]      (default: the main-text cases; see CASES)
Output: ~/Dropbox/pub/ms/2026-Thomson/TeX/figures/fig_spectrum_<case>.pdf and fig_spectrum_colorbar.pdf
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'py_scripts'))
sys.path.insert(0, os.path.join(REPO, 'py_scripts', 'radiation'))
sys.argv, _ARGS = sys.argv[:1], sys.argv[1:]   # the py_scripts modules read sys.argv at import
import plot_spectrum_map as psm  # noqa: E402
from utils.scaled_parameters import scaled_parameters  # noqa: E402

RUNS = os.path.expanduser('~/output/paper_2026/prelim/fig_spectrum')
OUT = os.path.expanduser('~/Dropbox/pub/ms/2026-Thomson/TeX/figures')

# case -> (label, run-id suffix); the 10 um and 800 nm ids are '10um_' + suffix and '800nm_' + suffix. The label is
# not drawn (the caption names the cases); it stays for reference.
CASES = {
    'a0_0p1_N1': (r'$a_0=0.1$, $N_c=30$, $N=1$', 'a0_0p1_N1'),
    'a0_0p1_N2': (r'$a_0=0.1$, $N_c=30$, $N=2$', 'a0_0p1_N2'),
    'a0_1p0_Nc8_N1': (r'$a_0=1$, $N_c=8$, $N=1$', 'a0_1p0_Nc8_N1'),
    'a0_1p0_Nc8_N2': (r'$a_0=1$, $N_c=8$, $N=2$', 'a0_1p0_Nc8_N2'),
    'a0_0p5_Nc30_N1': (r'$a_0=0.5$, $N_c=30$, $N=1$', 'a0_0p5_Nc30_N1'),
    'a0_0p5_Nc30_N2': (r'$a_0=0.5$, $N_c=30$, $N=2$', 'a0_0p5_Nc30_N2'),
    'a0_1p0_N1': (r'$a_0=1$, $N_c=30$, $N=1$', 'a0_1p0_N1'),          # first batch, line too short
    'a0_1p0_N2': (r'$a_0=1$, $N_c=30$, $N=2$', 'a0_1p0_N2'),
}
MAIN = ['a0_0p1_N1', 'a0_1p0_Nc8_N1', 'a0_1p0_Nc8_N2']

PANEL_SIZE = (2.33, 2.1)        # inches: three panels across a 7 in double column
MAP_RANGE = 1e-4                # colour scale of |E_sigma|^2 / max
FLUX_FRACTION = 0.999           # theta axis: up to where this fraction of the flux is reached
SPECTRUM_FLOOR = 1e-3           # omega axis: where the integrated spectrum exceeds this fraction of its peak ...
OMEGA_CROP = {'a0_0p1_N1': 0.05, 'a0_0p1_N2': 0.05}   # ... except here: N omega_1 x (1 -+ crop)
plt.rcParams.update({'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'legend.fontsize': 6, 'mathtext.fontset': 'cm', 'axes.linewidth': 0.6})


def latest_run(run_id):
    runs = sorted(glob.glob(os.path.join(RUNS, run_id, '*/')))
    if not runs:
        raise FileNotFoundError(f'no run of {run_id} under {RUNS}')
    return runs[-1]


def load(run_id):
    """theta (rad), omega/omega_1, |E_sigma|^2 map [omega, theta], helicity label, dN/d(omega/omega_1) (arb.)."""
    run = latest_run(run_id)
    theta, omega, maps, _ = psm.build_maps(os.path.join(run, 'radiation_field.dat'))
    E_plus, E_minus = maps['E_plus'] ** 2, maps['E_minus'] ** 2
    dominant, label = (E_minus, '-') if np.nansum(E_minus) >= np.nansum(E_plus) else (E_plus, '+')
    dN = np.nansum(maps['dPz'] * psm.annulus_areas(theta)[None, :], axis=1) / omega
    flux_theta = np.nansum(maps['dPz'] * psm.annulus_areas(theta)[None, :], axis=0)
    sp = scaled_parameters(run)
    return dict(run=run, theta=theta, theta_star=sp.theta_star, omega=omega, E2=dominant, sign=label, dN=dN,
                flux_theta=flux_theta)


def panel(case, legend=False):
    label, suffix = CASES[case]
    a, b = load('10um_' + suffix), load('800nm_' + suffix)
    theta_s = a['theta'] / a['theta_star']          # theta/theta_* (spec 11)
    cum = np.cumsum(a['flux_theta']) / a['flux_theta'].sum()
    theta_max = theta_s[min(np.searchsorted(cum, FLUX_FRACTION), len(cum) - 1)]
    if case in OMEGA_CROP:
        N = int(case[-1])
        w_lo, w_hi, pad = N * (1 - OMEGA_CROP[case]), N * (1 + OMEGA_CROP[case]), 0.0
    else:
        above = np.nonzero(a['dN'] > SPECTRUM_FLOOR * a['dN'].max())[0]
        w_lo, w_hi = a['omega'][above[0]], a['omega'][above[-1]]
        pad = 0.05 * (w_hi - w_lo)

    fig = plt.figure(figsize=PANEL_SIZE)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.06, left=0.2, right=0.97, bottom=0.17, top=0.97)
    ax_map = fig.add_subplot(gs[0])
    ax_spec = fig.add_subplot(gs[1], sharey=ax_map)

    E2 = a['E2'] / np.nanmax(a['E2'])
    ax_map.pcolormesh(theta_s, a['omega'], np.clip(E2, MAP_RANGE, None), norm=LogNorm(MAP_RANGE, 1.0),
                      cmap='viridis', shading='nearest', rasterized=True)
    ax_map.set_xlim(0, theta_max)
    ax_map.set_ylim(w_lo - pad, w_hi + pad)
    ax_map.set_xlabel(r'$\theta/\theta_*$')
    ax_map.set_ylabel(r'$\omega/\omega_1$')
    ax_map.text(0.96, 0.04, rf'$|E_{a["sign"]}|^2$', transform=ax_map.transAxes, color='w', va='bottom', ha='right',
                bbox=dict(boxstyle='round,pad=0.15', fc='k', ec='none', alpha=0.45))
    ax_map.xaxis.set_major_locator(plt.MaxNLocator(3))

    ax_spec.plot(a['dN'] / a['dN'].max(), a['omega'], color='k', lw=0.8, label=r'10 $\mu$m')
    ax_spec.plot(b['dN'] / b['dN'].max(), b['omega'], color='tab:orange', lw=0.8, ls=(0, (2, 1.5)), label='800 nm')
    ax_spec.set_xlim(0, 1.08)
    ax_spec.set_xticks([0, 1])
    ax_spec.set_xlabel(r'$dN/d\omega$ [norm.]')
    ax_spec.tick_params(labelleft=False)
    if legend:   # first panel only: the spectra fill the others
        ax_spec.legend(loc='upper right', frameon=False, handlelength=1.4, borderaxespad=0.1, labelspacing=0.2)

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f'fig_spectrum_{case}.pdf')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f'{path}\n  10 um: {a["run"]}\n  800 nm: {b["run"]}')


def colorbar():
    fig = plt.figure(figsize=(0.45, PANEL_SIZE[1]))
    ax = fig.add_axes([0.08, 0.17, 0.22, 0.80])
    sm = plt.cm.ScalarMappable(norm=LogNorm(MAP_RANGE, 1.0), cmap='viridis')
    fig.colorbar(sm, cax=ax, label=r'$|E_\sigma|^2/\max$')
    path = os.path.join(OUT, 'fig_spectrum_colorbar.pdf')
    fig.savefig(path)
    plt.close(fig)
    print(path)


if __name__ == '__main__':
    for i, case in enumerate(_ARGS or MAIN):
        panel(case, legend=(i == 0))
    colorbar()
