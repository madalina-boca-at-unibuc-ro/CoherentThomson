"""Builds the report's ponderomotive-lens figure (fig8) from the lattice run at a0 = 0.1 with the report's pure
Gaussian pulse, electrons at rest (~/output_tests/ponderomotive_lens/, theory/future_dev/ponderomotive_lens.md, Results).

Top: the dressing factor d_N(r) of one electron (amplitude, and phase relative to the axis). Then one row per
harmonic: the azimuthally averaged |E_-|^2 on the screen and its wavefront phase (along phi = 0, with the phase of a
beam focused at the interaction point, k r^2/(2Z), removed): the simulation, and the coherent-source prediction
without and with the dressing factor, each prediction with its best complex scale (the only free parameter).

Run from this folder: python3 make_lens_figure.py
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'py_scripts'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'py_scripts', 'radiation'))
import lg_decomposition as L
from faraday_frame_utils import C_LIGHT, convert_unit_to_number
from plot_field import get_screen_coordinates, read_length
from plot_observables import _get_E_total, _read_fundamental_frequency_au
import pandas as pd

RUN = os.path.expanduser('~/output_tests/ponderomotive_lens/20261030_111559')
OUT = 'figures/fig8_ponderomotive_lens.jpg'
INK, MUTED = '#2b2b2b', '#8a8a85'
BLUE, ORANGE = '#2a78d6', '#eb6834'      # categorical slots 1 and 2 (validated pair)

cfg = os.path.join(RUN, 'config.cfg')
field_path = os.path.join(RUN, 'radiation_field.dat')
model = L._source_model(RUN, cfg)
lambda_L = convert_unit_to_number('lambda', cfg)
w0 = model['w0']
s = read_length('circular_detector_distance', cfg, 'a.u.')
omega_1 = _read_fundamental_frequency_au(RUN)
x, y, _, _ = get_screen_coordinates(field_path, cfg)
x_au, y_au = (v * convert_unit_to_number(L.DETECTOR_AXES_UNIT, cfg) for v in (x, y))
r, phi = np.hypot(x_au, y_au), np.arctan2(y_au, x_au)
cols = open(field_path).readlines()[3].split()
data = pd.read_csv(field_path, sep=r'\s+', names=cols, header=None, skiprows=4)


def dressing(rr, N, k):
    A_r = 0.5 * model['a0'] ** 2 * L.laser_intensity_profile(rr, model['p'], model['l'], w0)
    return L.dressing_factor(rr, N, k * lambda_L / (2 * np.pi), A_r, model['envelope'])


def radial_mean(power):
    rho, index = np.unique(np.round(r, 6), return_inverse=True)
    return rho, np.bincount(index, power) / np.bincount(index)


fig, axes = plt.subplots(3, 2, figsize=(8.0, 9.6))
rr = np.linspace(0, 2 * w0, 400)
for N, color in ((1, BLUE), (2, ORANGE)):
    subset = data[data['i_omega'] == N - 1].sort_values('i_screen')
    k = subset['omega'].iloc[0] * omega_1 / C_LIGHT
    d = dressing(rr, N, k)
    axes[0, 0].plot(rr / w0, np.abs(d), color=color, lw=2, label=f'N = {N}')
    axes[0, 1].plot(rr / w0, np.unwrap(np.angle(d / d[0])), color=color, lw=2, label=f'N = {N}')

    Ex, Ey, _ = _get_E_total(subset)
    E = L.helicity_components(Ex, Ey)[1]                               # E_-, l = 1 - N
    ax, ax_phase = axes[N]
    ray = np.isclose(phi, 0.0)                                          # one ray, phi = 0, ordered in r
    order = np.argsort(r[ray])
    r_ray = r[ray][order]
    theta_max = 45 if N == 1 else 30
    # Up to where the power is small; without r = 0 for N >= 2, where the l = 1 - N vortex leaves no phase.
    shown = (r_ray / s * 1e3 <= (25 if N == 1 else 28)) & ((r_ray > 0) | (N == 1))

    def wavefront(field):
        f = field[ray][order] * np.exp(-0.5j * k * r_ray ** 2 / s)
        phase = np.unwrap(np.angle(f[shown]))
        return phase - phase[np.argmax(np.abs(f[shown]))]

    rho, measured = radial_mean(np.abs(E) ** 2)
    norm = measured.max()
    theta = rho / s * 1e3                                              # mrad
    ax.plot(theta, measured / norm, color=INK, lw=2.4, label='simulation')
    ax_phase.plot(r_ray[shown] / s * 1e3, wavefront(E), color=INK, lw=2.4, label='simulation')
    for dressed, style, label in ((False, dict(color=ORANGE, lw=2, ls='--'), 'prediction, no dressing'),
                                  (True, dict(color=BLUE, lw=1.6), 'prediction with $d_N(r)$')):
        pred, _ = L.coherent_source_prediction(r, phi, s, k, N, 1 - N, 0, model['p'], w0, model['R_beam'],
                                               dressing=(lambda q, N=N, k=k: dressing(q, N, k)) if dressed else None)
        scale = np.vdot(pred, E) / np.vdot(pred, pred)
        _, predicted = radial_mean(np.abs(scale * pred) ** 2)
        ax.plot(theta, predicted / norm, **style, label=label)
        ax_phase.plot(r_ray[shown] / s * 1e3, wavefront(scale * pred), **style, label=label)
    ax.set_xlim(0, theta_max)
    ax.set_title(f'N = {N}: power on the screen, $|E_-|^2$', fontsize=10, color=INK, loc='left')
    ax.set_ylabel('azimuthal mean (simulation peak = 1)')
    ax_phase.set_xlim(0, theta_max)
    ax_phase.set_title(f'N = {N}: wavefront phase on the screen', fontsize=10, color=INK, loc='left')
    ax_phase.set_ylabel('rad (0 at the power peak)')
    for a in (ax, ax_phase):
        a.set_xlabel('angle from the axis $r/Z$ [mrad]')

axes[0, 0].set_title('dressing factor: amplitude $|d_N(r)|$', fontsize=10, color=INK, loc='left')
axes[0, 1].set_title('dressing factor: phase, $\\arg[d_N(r)/d_N(0)]$', fontsize=10, color=INK, loc='left')
axes[0, 1].set_ylabel('rad')
for ax in axes[0]:
    ax.set_xlabel('electron radius $r/w_0$')
    ax.axvline(np.sqrt(np.log(2) / 2), color=MUTED, lw=1, ls=':')
axes[0, 0].text(np.sqrt(np.log(2) / 2) + 0.04, 0.6, 'resonant\nradius', color=MUTED, fontsize=8)
axes[0, 0].set_ylim(0.5, 1.02)
for ax in axes.flat:
    ax.grid(color='#e6e6e3', lw=0.8)
    ax.tick_params(colors=MUTED, labelsize=8)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(MUTED)
    ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)
    ax.xaxis.label.set_size(9); ax.yaxis.label.set_size(9)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
for ax in axes[1:, 1]:
    ax.get_legend().remove()
axes[2, 0].legend(frameon=False, fontsize=8, labelcolor=INK, loc='upper right')
axes[2, 0].set_ylim(top=1.5)
fig.tight_layout()
fig.savefig(OUT, dpi=200, pil_kwargs={'quality': 85, 'optimize': True, 'progressive': True})
print('wrote', OUT)
