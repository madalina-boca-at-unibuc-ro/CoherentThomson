#!/usr/bin/env python3
"""
Sec. IV.C figure (emitted field and angular-momentum fluxes): one PDF row per case, assembled into one figure in the
manuscript (TeX/sections/04-flat_sheet.tex). Each row, at the frequency of the spectral peak of that harmonic:
  (left)   Re[E_- e^{-ik(R-Z)}], the dominant helicity component with the spherical carrier about the origin removed,
           normalized to max|E_-|, on the screen in lab (x, y) / (D theta_*);
  (middle) arg[E_- e^{-ik(R-Z)}], masked where |E_-| < PHASE_MASK max|E_-|; the predicted l is written on it;
  (right)  radial profiles per unit theta/theta_* of the normal energy flux P_z (normalized to its screen integral)
           and of the spin and orbital angular-momentum fluxes Sigma_zz, Lambda_zz times omega, normalized by the
           same energy integral, so that their areas are the spin and orbital AM per photon (hbar = 1).
The 800 nm twins agree with these maps to <= 2.6% (L2; checks/maps_twin_compare.py) and are not drawn.

Data: the 10 um runs of reports/2026_thomson_paper/prelim/fig_maps/ (latest run of each id under
~/output/paper_2026/prelim/fig_maps/<id>/).

Usage: python3 fig_maps.py [case ...] [--preview PNG]   (default: both cases; --preview also writes the rows stacked)
Output: ~/Dropbox/pub/ms/2026-Thomson/TeX/figures/fig_maps_<case>.pdf
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'py_scripts'))
sys.path.insert(0, os.path.join(REPO, 'py_scripts', 'radiation'))
sys.argv, _ARGS = sys.argv[:1], sys.argv[1:]   # the py_scripts modules read sys.argv at import
from plot_field import read_config_value, read_length, get_screen_coordinates, DETECTOR_AXES_UNIT  # noqa: E402
from faraday_frame_utils import C_LIGHT, convert_unit_to_number, get_detector_local_rotation  # noqa: E402
from field_components import helicity_components  # noqa: E402
from plot_observables import (compute_observables, get_screen_area_weights_au,  # noqa: E402
                              _read_fundamental_frequency_au)
from lg_decomposition import _get_E_total, screen_lab_xy  # noqa: E402
from utils.scaled_parameters import scaled_parameters  # noqa: E402

RUNS = os.path.expanduser('~/output/paper_2026/prelim/fig_maps')
OUT = os.path.expanduser('~/Dropbox/pub/ms/2026-Thomson/TeX/figures')

# case -> (run id, harmonic N, predicted l = N(m + sigma) - sigma with m = -2, sigma = -1)
CASES = {
    'N1': ('10um_a0_1p0_Nc8_N1', 1, -2),
    'N2': ('10um_a0_1p0_Nc8_N2', 2, -5),
}
ROW_SIZE = (7.0, 2.2)           # inches: one row of three panels across a double column
# Row layout in inches (row_axes): margins, colorbar width and gaps; the maps are square, as high as the row allows.
ROW_BOTTOM, ROW_TOP, LEFT, RIGHT = 0.2, 0.17, 0.27, 0.12
CBAR_PAD, CBAR_W, CBAR_GAP, PROFILE_PAD = 0.05, 0.08, 0.42, 0.12
THETA_MAX = 12.0                # map half-width and profile range, in theta_* (99% of the power within 11.1 / 9.8)
PHASE_MASK = 1e-3               # phase masked below this fraction of max|E_-| (outside the emission only)
LABEL_BOX = dict(boxstyle='round,pad=0.15', fc='w', ec='none', alpha=0.85)
plt.rcParams.update({'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'legend.fontsize': 6, 'mathtext.fontset': 'cm', 'axes.linewidth': 0.6})


def latest_run(run_id):
    runs = sorted(glob.glob(os.path.join(RUNS, run_id, '*/')))
    if not runs:
        raise FileNotFoundError(f'no run of {run_id} under {RUNS}')
    return runs[-1]


def load(run_id):
    """E_- (carrier removed) on the (N_R + 1, N_phi) node grid with the axis prepended, lab (x, y) in theta_*, and
    the radial flux profiles per unit theta/theta_*."""
    run = latest_run(run_id)
    field = os.path.join(run, 'radiation_field.dat')
    cfg = os.path.join(run, 'config.cfg')
    data = pd.read_csv(field, sep=' ', comment='#').sort_values('i_screen')
    if data['i_omega'].nunique() != 1:
        raise ValueError(f'{run}: expected one frequency')
    N_R = int(read_config_value('circular_detector_N_R', cfg)[0])
    N_phi = int(read_config_value('circular_detector_N_phi', cfg)[0])
    s = abs(read_length('circular_detector_distance', cfg, 'a.u.'))
    theta_star = scaled_parameters(run).theta_star

    x_loc, y_loc, _, _ = get_screen_coordinates(field, cfg)
    to_au = convert_unit_to_number(DETECTOR_AXES_UNIT, cfg)
    x, y = screen_lab_xy(x_loc * to_au, y_loc * to_au, get_detector_local_rotation(cfg))
    omega_ratio = data['omega'].iloc[0]
    omega_au = omega_ratio * _read_fundamental_frequency_au(run)
    k = omega_au / C_LIGHT
    Ex, Ey, _ = _get_E_total(data)
    _, E_minus = helicity_components(Ex, Ey)
    E_minus = E_minus * np.exp(-1j * k * (np.sqrt(s ** 2 + x ** 2 + y ** 2) - s))
    shape = (N_R, N_phi)
    X, Y = x.reshape(shape) / (s * theta_star), y.reshape(shape) / (s * theta_star)
    E = E_minus.reshape(shape)
    # The field of a vortex (l != 0) vanishes on the axis: prepend r = 0 so the map has no central hole.
    X = np.vstack([np.zeros((1, N_phi)), X])
    Y = np.vstack([np.zeros((1, N_phi)), Y])
    E = np.vstack([np.zeros((1, N_phi), complex), E])

    obs, _, _ = compute_observables(field)
    obs = obs.sort_values('i_screen')
    w = get_screen_area_weights_au('circular', cfg).reshape(shape)
    ring = {q: (w * obs[q].to_numpy().reshape(shape)).sum(axis=1) for q in ('dPz', 'dSigmazz', 'dLambdazz')}
    r = np.hypot(x, y).reshape(shape)[:, 0]
    w_r2 = 2.0 * w.sum(axis=1) / (2.0 * np.pi)              # trapezoid weights in r^2
    dt = w_r2 / (2.0 * r) / (s * theta_star)                 # ring width in theta/theta_*
    P_tot = ring['dPz'].sum()
    prof = dict(t=r / (s * theta_star),
                P=ring['dPz'] / P_tot / dt,
                S=omega_au * ring['dSigmazz'] / P_tot / dt,
                L=omega_au * ring['dLambdazz'] / P_tot / dt,
                S_int=omega_au * ring['dSigmazz'].sum() / P_tot,
                L_int=omega_au * ring['dLambdazz'].sum() / P_tot)
    return dict(run=run, X=X, Y=Y, E=E / np.abs(E).max(), omega=omega_ratio, prof=prof)


def draw_row(fig, axes, case, cbar=True):
    run_id, N, l_pred = CASES[case]
    d = load(run_id)
    ax_re, ax_ph, ax_pr = axes[:3]
    lim = (-THETA_MAX, THETA_MAX)

    m_re = ax_re.pcolormesh(d['X'], d['Y'], d['E'].real, cmap='RdBu_r', vmin=-1, vmax=1, shading='gouraud',
                            rasterized=True)
    phase = np.ma.masked_where(np.abs(d['E']) < PHASE_MASK, np.angle(d['E']))
    m_ph = ax_ph.pcolormesh(d['X'], d['Y'], phase, cmap='twilight', vmin=-np.pi, vmax=np.pi, shading='nearest',
                            rasterized=True)
    for ax in (ax_re, ax_ph):
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect('equal')
        ax.set_xticks([-10, -5, 0, 5])
        ax.set_yticks([-10, -5, 0, 5, 10])
        end_label_x(ax, r'$\tilde\theta_x$')
    end_label_y(ax_re, r'$\tilde\theta_y$')
    ax_ph.tick_params(labelleft=False)
    ax_re.text(0.03, 0.97, rf'$N={N}$', transform=ax_re.transAxes, va='top', ha='left', bbox=LABEL_BOX)
    ax_re.text(0.97, 0.97, r'Re $E_-$', transform=ax_re.transAxes, va='top', ha='right', bbox=LABEL_BOX)
    ax_ph.text(0.97, 0.97, rf'arg $E_-$, $l={l_pred}$', transform=ax_ph.transAxes, va='top', ha='right',
               bbox=LABEL_BOX)
    if cbar:
        fig.colorbar(m_re, cax=axes[3], ticks=[-1, 0, 1])
        cb = fig.colorbar(m_ph, cax=axes[4], ticks=[-np.pi, 0, np.pi])
        cb.ax.set_yticklabels([r'$-\pi$', '0', r'$\pi$'])

    p = d['prof']
    ax_pr.axhline(0, color='0.6', lw=0.4)
    ax_pr.plot(p['t'], p['P'], color='k', lw=0.8, label=r'$\bar P_z$ (1)')
    ax_pr.plot(p['t'], p['S'], color='tab:blue', lw=0.8, ls=(0, (3, 1.5)),
               label=rf'$\omega\bar\Sigma_{{zz}}$ (${p["S_int"]:.3f}$)')
    ax_pr.plot(p['t'], p['L'], color='tab:red', lw=0.8,
               label=rf'$\omega\bar\Lambda_{{zz}}$ (${p["L_int"]:.3f}$)')
    ax_pr.set_xlim(0, THETA_MAX)
    ax_pr.set_xticks(np.arange(0, THETA_MAX - 1, 2))
    end_label_x(ax_pr, r'$\tilde\theta$')
    end_label_y(ax_pr, r'$d\bar Q/d\tilde\theta$', ha='center')
    if N > 1:   # the curves fill the lower part: make room for the legend above them
        lo, hi = ax_pr.get_ylim()
        ax_pr.set_ylim(lo, hi + 0.45 * (hi - lo))
    ax_pr.legend(loc='upper right' if N > 1 else 'lower right', frameon=False, handlelength=1.6,
                 borderaxespad=0.2, labelspacing=0.25)
    print(f'{case}: {d["run"]}  omega/omega_1 = {d["omega"]}  spin {p["S_int"]:+.6f}  orbital {p["L_int"]:+.6f}'
          f'  (predicted {-1:+d}, {l_pred:+d})')


def end_label_x(ax, text):
    """Axis label at the right end of the x axis, level with the tick labels (the last tick stays unlabelled)."""
    ax.annotate(text, xy=(1, 0), xycoords='axes fraction', xytext=(0, -7), textcoords='offset points',
                ha='right', va='center')


def end_label_y(ax, text, ha='right'):
    """Axis label above the top end of the y axis."""
    ax.annotate(text, xy=(0, 1), xycoords='axes fraction', xytext=(-2 if ha == 'right' else 0, 3),
                textcoords='offset points', ha=ha, va='bottom')


def row_axes(fig, y0):
    """The five axes of one row (Re map, arg map, profiles, two colorbars), laid out in inches from the bottom y0 of
    the row; returned in draw_row's order (maps, profile, colorbars)."""
    W, H = fig.get_size_inches()
    h = ROW_SIZE[1] - ROW_BOTTOM - ROW_TOP              # map side = panel height
    x, boxes = LEFT, []
    for _ in range(2):                                  # map, gap, colorbar, colorbar tick labels + gap
        boxes.append((x, h))
        boxes.append((x + h + CBAR_PAD, CBAR_W))
        x += h + CBAR_PAD + CBAR_W + CBAR_GAP
    x += PROFILE_PAD
    boxes.append((x, W - RIGHT - x))
    ax = [fig.add_axes([bx / W, (y0 + ROW_BOTTOM) / H, bw / W, h / H]) for bx, bw in boxes]
    return [ax[0], ax[2], ax[4], ax[1], ax[3]]


def row_figure(case):
    fig = plt.figure(figsize=ROW_SIZE)
    axes = row_axes(fig, 0.0)
    draw_row(fig, axes, case)
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f'fig_maps_{case}.pdf')
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(path)


def preview(cases, path):
    fig = plt.figure(figsize=(ROW_SIZE[0], ROW_SIZE[1] * len(cases)))
    for i, case in enumerate(cases):
        draw_row(fig, row_axes(fig, ROW_SIZE[1] * (len(cases) - 1 - i)), case)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(path)


if __name__ == '__main__':
    args = list(_ARGS)
    preview_path = None
    if '--preview' in args:
        i = args.index('--preview')
        preview_path = args[i + 1]
        del args[i:i + 2]
    cases = args or list(CASES)
    for case in cases:
        row_figure(case)
    if preview_path:
        preview(cases, preview_path)
