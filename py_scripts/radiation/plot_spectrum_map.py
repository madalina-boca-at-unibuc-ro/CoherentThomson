"""
Spectrum along a screen line (theory/specs/15_spectrum_map.md): maps over (theta, omega/omega_1) of a
dense_frequency_spectrum run on a line of screen points (rectangular Nx or Ny = 1, circular N_phi = 1).

- Map figure: dP_z/domega, |E_+| and |E_-| (total field), log colour from max/R to max (|E_+| and |E_-| share
  one scale). E_+- = (E_x -+ i E_y)/sqrt2,
  spin +-1 about lab +z (lg_decomposition.helicity_components); dP_z/domega as plot_point_spectrum_observables.py.
  These don't depend on phi for an axisymmetric source, so one azimuth stands for the screen; Cartesian F components
  do (both helicities present), which is why they aren't mapped. Solid lines: the computed harmonics N omega_1.
  Dashed: N omega_1 (1 + a0^2/2), the undressed limit (head-on, on axis; the angle dependence is negligible at
  gamma theta << 1). Dotted: the red edge of the laser's actual peak <a^2> (see the summary below).
- Line-integrated spectrum: dW/domega = sum_i dP_z/domega(theta_i) A_i, A_i the annulus of the i-th point (edges
  halfway between points), assuming axisymmetry: the beam's spectrum through the screen's disk.
- Summary per harmonic window, printed and written to run_log.txt under its own section: the integrated spectrum's
  peak and its offset from N, FWHM, the central 90% band, the energy fractions within +-1% of N omega_1 and in the
  reachable band, theta of the dP_z/domega maximum at N omega_1 and at the peak, the energy-weighted centroid, the
  local maxima above 10% of the window maximum, and a ready omega_list line of those maxima (spec 16; they are peaks
  of the integrated spectrum, not of a given island of the map: pick others from the map). Harmonic N can reach
  N x [(1 + a0^2/2)/(1 + g a0^2/2), 1 + a0^2/2] omega_1 (g = the laser's peak |a|^2/a0^2, 0.541 for the code's LG
  (0, 2); omega_1 itself is dressed with g = 1); windows split halfway between neighbouring bands.

Spectra at a fixed theta: single-point runs and plot_point_spectrum*.py.

Usage: python3 plot_spectrum_map.py [run_folder] [--range R] [--omega-min W] [--omega-max W] [--harmonic N] [--linear]
--linear (colour from 0 to max) needs a plotted range holding one harmonic: --harmonic N, or --omega-min/--omega-max.
PNGs: spectrum_map{,_integrated}[_N<N>][_linear].png.
"""
import argparse
import math
import os
import sys

import numpy as np
import scipy.special
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


def peak_a2_fraction(config_path):
    """Peak local |a|^2 / a0^2 of the laser's transverse profile at the waist: 1 for a plane wave; for LG (p, l) the
    maximum over rho of |Npn (sqrt2 rho/w0)^|l| 1F1(-p, |l|+1, 2 rho^2/w0^2) exp(-rho^2/w0^2)|^2, the C++
    LaguerreGaussLaser amplitude with its custom Npn (0.541 for (0, 2), 2 on axis for (0, 0))."""
    laser_type, _ = read_config_value('laser_type', config_path)
    if laser_type != 'laguerre_gauss':
        return 1.0
    p = int(read_config_value('laser_lg_p', config_path)[0])
    n = abs(int(read_config_value('laser_lg_l', config_path)[0]))
    Npn = np.sqrt(2.0) / math.factorial(n) * np.sqrt(math.factorial(p + n) / math.factorial(p))
    u = np.linspace(0.0, 10.0 + 4.0 * (p + n), 200001)  # rho^2 / w0^2
    amplitude = Npn * (2.0 * u) ** (n / 2.0) * scipy.special.hyp1f1(-p, n + 1, 2.0 * u) * np.exp(-u)
    return float(np.max(amplitude ** 2))


def harmonic_bands(a0, g):
    """(lower, upper) in omega/omega_1 per unit N: the emission a peak local <a^2> = g a0^2/2 can reach, from its
    red edge (1 + a0^2/2)/(1 + g a0^2/2) (head-on, on axis) up to the undressed 1 + a0^2/2."""
    undressed = 1.0 + a0 ** 2 / 2.0
    return undressed / (1.0 + g * a0 ** 2 / 2.0), undressed


def window_boundary(N, lower, upper):
    """omega/omega_1 between harmonic band N and band N + 1 (halfway between N upper and (N + 1) lower); 0 -> below N = 1."""
    return 0.5 * lower if N == 0 else 0.5 * (N * upper + (N + 1) * lower)


def harmonics_in(omega_min, omega_max, lower, upper):
    """The harmonics whose window overlaps [omega_min, omega_max]."""
    return [N for N in range(1, int(np.floor(omega_max / lower)) + 2)
            if window_boundary(N - 1, lower, upper) < omega_max and window_boundary(N, lower, upper) > omega_min]


def harmonic_lines(ax, omega, a0, g, vertical=False):
    """Solid at N omega_1 (computed); dashed at N omega_1 (1 + a0^2/2) (undressed); dotted at the red edge of the
    laser's actual peak <a^2> (harmonic_bands), when it differs from N omega_1. Inside the plotted range."""
    lower, upper = harmonic_bands(a0, g)
    draw = ax.axvline if vertical else ax.axhline
    color = 'tab:red' if vertical else 'white'
    for N in range(1, int(np.floor(omega.max() / lower)) + 1):
        for value, style in ((N, '-'), (N * upper, '--'), (N * lower, ':')):
            if style == ':' and abs(lower - 1.0) < 1e-6:
                continue
            if omega.min() <= value <= omega.max():
                draw(value, color=color, linewidth=0.8 if style != ':' else 1.2, alpha=0.8, linestyle=style)


def draw_maps(theta, omega, maps, a0, g, dynamic_range, png_dir, title, linear=False, suffix=''):
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
        if linear:
            image = ax.pcolormesh(theta * scale, omega, values, shading='nearest', vmin=0.0, vmax=vmax, cmap='viridis')
        else:
            vmin = vmax / dynamic_range
            image = ax.pcolormesh(theta * scale, omega, np.clip(values, vmin, None), shading='nearest',
                                  norm=LogNorm(vmin=vmin, vmax=vmax), cmap='viridis')
        harmonic_lines(ax, omega, a0, g)
        ax.set_title(label + (' (|value|)' if key == 'dPz' and np.nanmin(maps[key]) < 0 else ''))
        ax.set_xlabel(rf'$\theta$ [{unit}]')
        fig.colorbar(image, ax=ax, pad=0.02)
    axes[0].set_ylabel(OMEGA_LABEL)
    fig.suptitle(f"{title}\nsolid: $N\\omega_1$; dashed: $N\\omega_1(1 + a_0^2/2)$, undressed; dotted: red edge of the "
                 f"peak $\\langle a^2\\rangle = {g:.3g}\\,a_0^2/2$ (head-on, on axis)", fontsize=10)
    fig.tight_layout()
    path = os.path.join(png_dir, f'spectrum_map{suffix}.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Successfully saved plot to {path}")


def draw_integrated(omega, dW, a0, g, png_dir, title, suffix=''):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
    for ax, log in zip(axes, (False, True)):
        ax.plot(omega, np.abs(dW) if log else dW, color='tab:blue', linewidth=1.0)
        if log:
            ax.set_yscale('log')
            positive = np.abs(dW)[np.abs(dW) > 0]
            if positive.size:
                ax.set_ylim(positive.max() * 1e-12, positive.max() * 2)
        harmonic_lines(ax, omega, a0, g, vertical=True)
        ax.set_xlabel(OMEGA_LABEL)
        ax.set_ylabel(r'$dW/d\omega$ through the line' + ' (|value|)' * log + ' [a.u.]')
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title}: line-integrated spectrum (axisymmetric), $\\sum_i dP_z/d\\omega\\, A_i$", fontsize=10)
    fig.tight_layout()
    path = os.path.join(png_dir, f'spectrum_map_integrated{suffix}.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Successfully saved plot to {path}")


def local_peaks(omega, values, fraction):
    """(omega, height / max) of the local maxima of `values` above `fraction` of its maximum, highest first: a sample
    larger than its neighbours on the scan's own grid, no smoothing; a flat top counts once (its first sample)."""
    top = np.nanmax(values) if len(values) else 0.0
    if not np.isfinite(top) or top <= 0.0:
        return []
    peaks = []
    for i in range(1, len(values) - 1):
        j = i
        while j + 1 < len(values) - 1 and values[j + 1] == values[i]:
            j += 1
        if values[i] > values[i - 1] and values[i] > values[j + 1] and values[i] >= fraction * top:
            peaks.append((omega[i], values[i] / top))
    return sorted(peaks, key=lambda peak: -peak[1])


def summary_lines(theta, omega, dPz, dW, a0, g):
    scale, unit = angle_scale(theta.max())
    unit = unit.replace(r'$\mu$', 'u')
    lines = [f"  integrated over {len(theta)} points, theta {theta[0] * scale:.4g} to {theta[-1] * scale:.4g} {unit}; "
             f"omega/omega_1 {omega[0]:.6g} to {omega[-1]:.6g}, step {omega[1] - omega[0]:.3g}"]
    lower, upper = harmonic_bands(a0, g)
    all_peaks = []
    lines.append(f"  peak local <a^2> = {g:.4g} a0^2/2 (laser profile): harmonic N reachable over N x [{lower:.6g}, "
                 f"{upper:.6g}] omega_1 (head-on, on axis); windows split halfway between neighbouring bands")

    def boundary(N):  # between band N and band N + 1
        return window_boundary(N, lower, upper)

    for N in range(1, int(np.floor(omega.max() / lower)) + 1):
        window = (omega >= boundary(N - 1)) & (omega < boundary(N))
        if window.sum() < 3:
            continue
        w_omega, w_dW = omega[window], dW[window]
        partial = omega[0] > boundary(N - 1) or omega[-1] < boundary(N) - (omega[1] - omega[0])
        band = (w_omega >= N * lower) & (w_omega <= N * upper)
        energy = w_dW.sum()
        i_peak = int(np.argmax(w_dW))
        above = np.nonzero(w_dW >= 0.5 * w_dW[i_peak])[0]
        cumulative = np.cumsum(w_dW) / energy if energy > 0 else np.zeros_like(w_dW)
        lo, hi = w_omega[np.searchsorted(cumulative, 0.05)], w_omega[min(np.searchsorted(cumulative, 0.95), len(w_omega) - 1)]
        near = np.abs(w_omega - N) <= 0.01 * N
        at_N = int(np.argmin(np.abs(omega - N)))
        at_peak = int(np.nonzero(window)[0][i_peak])
        peaks = local_peaks(w_omega, w_dW, 0.1)
        all_peaks += [(N, w) for w, _ in peaks]
        lines += [f"  N = {N}{' (window only partly scanned)' if partial else ''}:",
                  f"    integrated peak at omega/omega_1 = {w_omega[i_peak]:.6g} (offset from N {w_omega[i_peak] - N:+.4g})",
                  f"    FWHM span {w_omega[above[0]]:.6g} to {w_omega[above[-1]]:.6g} (width {w_omega[above[-1]] - w_omega[above[0]]:.4g})",
                  f"    central 90% of the window's energy: {lo:.6g} to {hi:.6g}",
                  f"    energy within +-1% of N omega_1: {w_dW[near].sum() / energy if energy > 0 else float('nan'):.4g} of the window",
                  f"    energy in the reachable band {N * lower:.6g} to {N * upper:.6g}: "
                  f"{w_dW[band].sum() / energy if energy > 0 else float('nan'):.4g} of the window",
                  f"    theta of max dP_z/domega: {theta[np.nanargmax(dPz[at_N])] * scale:.4g} {unit} at omega = N omega_1, "
                  f"{theta[np.nanargmax(dPz[at_peak])] * scale:.4g} {unit} at the peak",
                  f"    centroid (energy-weighted) at omega/omega_1 = {np.sum(w_omega * w_dW) / energy if energy > 0 else float('nan'):.6g}",
                  "    local maxima above 10% of the window maximum (omega/omega_1: height): "
                  + ", ".join(f"{w:.5g}: {h:.3g}" for w, h in peaks)]
    if all_peaks:
        ordered = sorted(all_peaks, key=lambda peak: (peak[0], peak[1]))
        lines.append("  omega_list line for a screen run at these peaks (theory/specs/16_omega_list.md; a peak of the "
                     "integrated spectrum, not of a given island of the map):")
        lines.append("    omega_list    " + " ".join(f"{N}:{w:.6g}" for N, w in ordered) + "   first_harmonic_frequency")
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('run_folder', nargs='?', default=None)
    parser.add_argument('--range', type=float, default=1e6, help='log colour range max/min (default 1e6)')
    parser.add_argument('--omega-min', type=float, default=None, help='crop the plots (omega/omega_1)')
    parser.add_argument('--omega-max', type=float, default=None, help='crop the plots (omega/omega_1)')
    parser.add_argument('--harmonic', type=int, default=None,
                        help='crop the plots to harmonic N: its window, narrowed to where its spectrum is above 1e-3 '
                             'of the peak (unless --omega-min/--omega-max are given)')
    parser.add_argument('--linear', action='store_true',
                        help='linear colour scale from 0; needs a range holding one harmonic (--harmonic or --omega-*)')
    args = parser.parse_args()

    radiation_filepath = resolve_radiation_file(args.run_folder)
    theta, omega, maps, config_path = build_maps(radiation_filepath)
    a0 = float(read_config_value('laser_a0', config_path)[0])
    g = peak_a2_fraction(config_path)
    detector_type, _ = read_config_value('detector_type', config_path)
    distance = read_length(f'{detector_type}_detector_distance', config_path)
    rho = distance * np.tan(theta)
    dW = np.nansum(maps['dPz'] * annulus_areas(rho)[None, :], axis=1)

    lower, upper = harmonic_bands(a0, g)
    keep = np.ones(len(omega), dtype=bool)
    if args.harmonic is not None:
        # The harmonic's window, then (unless --omega-* says otherwise) cropped to where its integrated spectrum is
        # above 1e-3 of its peak, padded by 25% of that span: at small a0 the line is far narrower than the window.
        keep &= (omega >= window_boundary(args.harmonic - 1, lower, upper)) & (
            omega < window_boundary(args.harmonic, lower, upper))
        if args.omega_min is None and args.omega_max is None and keep.any():
            strong = keep & (dW >= 1e-3 * np.nanmax(np.where(keep, dW, np.nan)))
            lo, hi = omega[strong].min(), omega[strong].max()
            pad = 0.25 * max(hi - lo, omega[1] - omega[0])
            keep &= (omega >= lo - pad) & (omega <= hi + pad)
    if args.omega_min is not None:
        keep &= omega >= args.omega_min
    if args.omega_max is not None:
        keep &= omega <= args.omega_max
    if not keep.any():
        print("Error: no frequency left after cropping.")
        sys.exit(1)
    harmonics = harmonics_in(omega[keep].min(), omega[keep].max(), lower, upper)
    if args.linear and len(harmonics) > 1:
        print(f"Error: --linear needs one harmonic in the plotted range, which holds N = {harmonics}: harmonics differ by "
              "orders of magnitude. Use --harmonic N or --omega-min/--omega-max.")
        sys.exit(1)
    suffix = (f'_N{args.harmonic}' if args.harmonic is not None else '') + ('_linear' if args.linear else '')
    run_dir = os.path.dirname(radiation_filepath)
    png_dir = os.path.join(run_dir, 'png_folder', MODULE_NAME, 'spectrum_map')
    os.makedirs(png_dir, exist_ok=True)
    title = f"{os.path.basename(run_dir)}: a0 = {a0:g}, {len(theta)} points"
    draw_maps(theta, omega[keep], {k: v[keep] for k, v in maps.items()}, a0, g, args.range, png_dir, title,
              linear=args.linear, suffix=suffix)
    draw_integrated(omega[keep], dW[keep], a0, g, png_dir, title, suffix=suffix)

    lines = summary_lines(theta, omega, maps['dPz'], dW, a0, g)
    print('\n'.join(lines))
    _write_run_log_section(run_dir, _RUN_LOG_SECTION_HEADER, lines)


if __name__ == "__main__":
    main()
