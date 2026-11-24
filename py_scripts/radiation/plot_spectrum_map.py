"""
Spectrum along a screen line (theory/specs/15_spectrum_map.md): maps over (theta, omega/omega_1) of a
dense_frequency_spectrum run on a line of screen points (rectangular Nx or Ny = 1, circular N_phi = 1).

- Map figure: dP_z/domega, |E_+| and |E_-| (total field), log colour from max/R to max (|E_+| and |E_-| share
  one scale). E_+- = (E_x -+ i E_y)/sqrt2,
  spin +-1 about lab +z (lg_decomposition.helicity_components); dP_z/domega as plot_point_spectrum_observables.py.
  These don't depend on phi for an axisymmetric source, so one azimuth stands for the screen; Cartesian F components
  do (both helicities present), which is why they aren't mapped. Solid lines: the computed harmonics N omega_1.
  Dashed: N omega_1 (1 + a0^2/2), the undressed limit (head-on, on axis; the angle dependence is negligible at
  gamma theta << 1).
- Line-integrated spectrum: dW/domega = sum_i dP_z/domega(theta_i) A_i, A_i the annulus of the i-th point (edges
  halfway between points), assuming axisymmetry: the beam's spectrum through the screen's disk.
- Summary per harmonic window omega/omega_1 in [N - 1/2, N + 1/2), printed and written to run_log.txt under its own
  section: the integrated spectrum's peak and its offset from N, FWHM, the central 90% band, the energy fraction
  within +-1% of N, and theta of the dP_z/domega maximum at N omega_1 and at the peak.

Spectra at a fixed theta: single-point runs and plot_point_spectrum*.py.

Usage: python3 plot_spectrum_map.py [run_folder] [--range R] [--omega-min W] [--omega-max W]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from plot_point_spectrum import OMEGA_LABEL, read_radiation_data, resolve_radiation_file
from plot_point_spectrum_observables import compute_point_observables, get_observation_directions
from plot_field import read_config_value, read_length
from plot_observables import _get_E_total, _write_run_log_section
from lg_decomposition import helicity_components
from faraday_frame_utils import get_detector_local_rotation

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))
_RUN_LOG_SECTION_HEADER = "Spectrum map (plot_spectrum_map.py)"


def require_line_run(config_path, n_screen):
    """Exits unless the run is a dense spectrum on a line of at least two points."""
    dense, _ = read_config_value('dense_frequency_spectrum', config_path)
    detector_type, _ = read_config_value('detector_type', config_path)
    counts = {'rectangular': ('rectangular_detector_Nx', 'rectangular_detector_Ny'),
              'circular': ('circular_detector_N_R', 'circular_detector_N_phi')}
    if dense != 'true' or detector_type not in counts:
        print("Error: needs a dense_frequency_spectrum run on a rectangular or circular line; for other runs use "
              "plot_field.py / plot_point_spectrum*.py.")
        sys.exit(1)
    n1, n2 = (int(read_config_value(key, config_path)[0]) for key in counts[detector_type])
    is_line = (n1 == 1) != (n2 == 1) if detector_type == 'rectangular' else (n2 == 1 and n1 > 1)
    if not is_line or n_screen < 2:
        print(f"Error: {n_screen} screen point(s) on a {n1} x {n2} {detector_type} grid is not a line; for a single "
              "point use plot_point_spectrum*.py.")
        sys.exit(1)


def annulus_areas(rho):
    """Area of the annulus each (sorted) radius stands for, edges halfway between neighbours."""
    if len(rho) == 1:
        return np.array([np.pi * rho[0] ** 2])
    mid = 0.5 * (rho[1:] + rho[:-1])
    edges = np.concatenate([[max(0.0, rho[0] - (mid[0] - rho[0]))], mid, [rho[-1] + (rho[-1] - mid[-1])]])
    return np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)


def build_maps(radiation_filepath):
    """theta (rad, sorted), omega/omega_1, and [omega, theta] arrays of dP_z/domega, |E_+|, |E_-|."""
    run_dir = os.path.dirname(radiation_filepath)
    config_path = os.path.join(run_dir, 'config.cfg')
    data = read_radiation_data(radiation_filepath)
    n_screen = data['i_screen'].nunique()
    require_line_run(config_path, n_screen)

    observables, _ = compute_point_observables(radiation_filepath, data)
    data = data[data['omega'] != 0.0]
    Ex, Ey, _ = _get_E_total(data)
    E_plus, E_minus = helicity_components(Ex, Ey)

    directions = get_observation_directions(radiation_filepath, config_path, n_screen)
    normal = get_detector_local_rotation(config_path)[:, 2]
    theta = np.arccos(np.clip(directions @ normal, -1.0, 1.0))
    order = np.argsort(theta)

    omega = np.sort(data['omega'].unique())
    i_omega = np.searchsorted(omega, data['omega'].to_numpy())
    i_screen = data['i_screen'].to_numpy()

    def as_map(values):
        grid = np.full((len(omega), n_screen), np.nan)
        grid[i_omega, i_screen] = values
        return grid[:, order]

    maps = {'dPz': as_map(observables['dPz'].to_numpy()), 'E_plus': as_map(np.abs(E_plus)),
            'E_minus': as_map(np.abs(E_minus))}
    return theta[order], omega, maps, config_path


def angle_scale(theta_max):
    if theta_max < 2e-3:
        return 1e6, r'$\mu$rad'
    if theta_max < 2.0:
        return 1e3, 'mrad'
    return 1.0, 'rad'


def harmonic_lines(ax, omega, a0, vertical=False):
    """Solid at N omega_1 (computed), dashed at N omega_1 (1 + a0^2/2) (undressed, head-on), inside the range."""
    undressed = 1.0 + a0 ** 2 / 2.0
    draw = ax.axvline if vertical else ax.axhline
    for N in range(1, int(np.floor(omega.max())) + 1):
        if omega.min() <= N <= omega.max():
            draw(N, color='white' if not vertical else 'tab:red', linewidth=0.8, alpha=0.8)
        if omega.min() <= N * undressed <= omega.max():
            draw(N * undressed, color='white' if not vertical else 'tab:red', linewidth=0.8, alpha=0.8,
                 linestyle='--')


def draw_maps(theta, omega, maps, a0, dynamic_range, png_dir, title):
    scale, unit = angle_scale(theta.max())
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), sharey=True)
    panels = (('dPz', r'$dP_z/d\omega$'), ('E_plus', r'$|E_+|$'), ('E_minus', r'$|E_-|$'))
    # |E_+| and |E_-| share one scale, so the panels show how strong each helicity is relative to the other.
    helicity_max = np.nanmax([np.nanmax(np.abs(maps['E_plus'])), np.nanmax(np.abs(maps['E_minus']))])
    for ax, (key, label) in zip(axes, panels):
        values = np.abs(maps[key])
        vmax = np.nanmax(values) if key == 'dPz' else helicity_max
        if not np.isfinite(vmax) or vmax <= 0.0:
            ax.set_title(f"{label}: zero")
            continue
        vmin = vmax / dynamic_range
        image = ax.pcolormesh(theta * scale, omega, np.clip(values, vmin, None), shading='nearest',
                              norm=LogNorm(vmin=vmin, vmax=vmax), cmap='viridis')
        harmonic_lines(ax, omega, a0)
        ax.set_title(label + (' (|value|)' if key == 'dPz' and np.nanmin(maps[key]) < 0 else ''))
        ax.set_xlabel(rf'$\theta$ [{unit}]')
        fig.colorbar(image, ax=ax, pad=0.02)
    axes[0].set_ylabel(OMEGA_LABEL)
    fig.suptitle(f"{title}\nsolid: $N\\omega_1$; dashed: $N\\omega_1(1 + a_0^2/2)$, undressed limit (head-on, on axis)",
                 fontsize=10)
    fig.tight_layout()
    path = os.path.join(png_dir, 'spectrum_map.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Successfully saved plot to {path}")


def draw_integrated(omega, dW, a0, png_dir, title):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
    for ax, log in zip(axes, (False, True)):
        ax.plot(omega, np.abs(dW) if log else dW, color='tab:blue', linewidth=1.0)
        if log:
            ax.set_yscale('log')
            positive = np.abs(dW)[np.abs(dW) > 0]
            if positive.size:
                ax.set_ylim(positive.max() * 1e-12, positive.max() * 2)
        harmonic_lines(ax, omega, a0, vertical=True)
        ax.set_xlabel(OMEGA_LABEL)
        ax.set_ylabel(r'$dW/d\omega$ through the line' + ' (|value|)' * log + ' [a.u.]')
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: line-integrated spectrum (axisymmetric), $\\sum_i dP_z/d\\omega\\, A_i$", fontsize=10)
    fig.tight_layout()
    path = os.path.join(png_dir, 'spectrum_map_integrated.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Successfully saved plot to {path}")


def summary_lines(theta, omega, dPz, dW):
    scale, unit = angle_scale(theta.max())
    unit = unit.replace(r'$\mu$', 'u')
    lines = [f"  integrated over {len(theta)} points, theta {theta[0] * scale:.4g} to {theta[-1] * scale:.4g} {unit}; "
             f"omega/omega_1 {omega[0]:.6g} to {omega[-1]:.6g}, step {omega[1] - omega[0]:.3g}"]
    for N in range(1, int(np.floor(omega.max() + 0.5)) + 1):
        window = (omega >= N - 0.5) & (omega < N + 0.5)
        if window.sum() < 3:
            continue
        w_omega, w_dW = omega[window], dW[window]
        partial = omega[0] > N - 0.5 or omega[-1] < N + 0.5 - (omega[1] - omega[0])
        energy = w_dW.sum()
        i_peak = int(np.argmax(w_dW))
        above = np.nonzero(w_dW >= 0.5 * w_dW[i_peak])[0]
        cumulative = np.cumsum(w_dW) / energy if energy > 0 else np.zeros_like(w_dW)
        lo, hi = w_omega[np.searchsorted(cumulative, 0.05)], w_omega[min(np.searchsorted(cumulative, 0.95), len(w_omega) - 1)]
        near = np.abs(w_omega - N) <= 0.01 * N
        at_N = int(np.argmin(np.abs(omega - N)))
        at_peak = int(np.nonzero(window)[0][i_peak])
        lines += [f"  N = {N}{' (window only partly scanned)' if partial else ''}:",
                  f"    integrated peak at omega/omega_1 = {w_omega[i_peak]:.6g} (offset from N {w_omega[i_peak] - N:+.4g})",
                  f"    FWHM span {w_omega[above[0]]:.6g} to {w_omega[above[-1]]:.6g} (width {w_omega[above[-1]] - w_omega[above[0]]:.4g})",
                  f"    central 90% of the window's energy: {lo:.6g} to {hi:.6g}",
                  f"    energy within +-1% of N omega_1: {w_dW[near].sum() / energy if energy > 0 else float('nan'):.4g} of the window",
                  f"    theta of max dP_z/domega: {theta[np.nanargmax(dPz[at_N])] * scale:.4g} {unit} at omega = N omega_1, "
                  f"{theta[np.nanargmax(dPz[at_peak])] * scale:.4g} {unit} at the peak"]
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('run_folder', nargs='?', default=None)
    parser.add_argument('--range', type=float, default=1e6, help='log colour range max/min (default 1e6)')
    parser.add_argument('--omega-min', type=float, default=None, help='crop the plots (omega/omega_1)')
    parser.add_argument('--omega-max', type=float, default=None, help='crop the plots (omega/omega_1)')
    args = parser.parse_args()

    radiation_filepath = resolve_radiation_file(args.run_folder)
    theta, omega, maps, config_path = build_maps(radiation_filepath)
    a0 = float(read_config_value('laser_a0', config_path)[0])
    detector_type, _ = read_config_value('detector_type', config_path)
    distance = read_length(f'{detector_type}_detector_distance', config_path)
    rho = distance * np.tan(theta)
    dW = np.nansum(maps['dPz'] * annulus_areas(rho)[None, :], axis=1)

    keep = np.ones(len(omega), dtype=bool)
    if args.omega_min is not None:
        keep &= omega >= args.omega_min
    if args.omega_max is not None:
        keep &= omega <= args.omega_max
    run_dir = os.path.dirname(radiation_filepath)
    png_dir = os.path.join(run_dir, 'png_folder', MODULE_NAME, 'spectrum_map')
    os.makedirs(png_dir, exist_ok=True)
    title = f"{os.path.basename(run_dir)}: a0 = {a0:g}, {len(theta)} points"
    draw_maps(theta, omega[keep], {k: v[keep] for k, v in maps.items()}, a0, args.range, png_dir, title)
    draw_integrated(omega[keep], dW[keep], a0, png_dir, title)

    lines = summary_lines(theta, omega, maps['dPz'], dW)
    print('\n'.join(lines))
    _write_run_log_section(run_dir, _RUN_LOG_SECTION_HEADER, lines)


if __name__ == "__main__":
    main()
