"""
Local radiation observables vs omega/omega_1 at the screen point(s) of a dense_frequency_spectrum run, from the
total (LR+SR+BR) field in radiation_field.dat. Same formulas and prefactors as plot_observables.py (see its doc
comment and theory/numerical_calculation_of_angular_momentum.md), restricted to the ones that need no spatial
derivative:

- dS_z/domega        spin angular momentum density, (4 eps0/omega) Im[Ex* Ey]
- dSigma_zz/domega   spin angular momentum flux through the z-plane
- du/domega          energy density
- dP_z/domega        energy flux (Poynting vector, z component)

The orbital quantities L_z/Lambda_zz need x d/dy - y d/dx (or d/dphi) across the screen, so they can't be
computed at a single point; use plot_observables.py on a gridded run for them.

Derived, per frequency:
- omega Sigma_zz/P_z and omega S_z/u: spin per photon about lab +z (hbar = 1). Both equal h n_z for radiation of
  helicity h moving along n.
- V n_z, with V = 2 Im[E1* E2]/(|E1|^2 + |E2|^2) the degree of circular polarization about the observation
  direction n = (screen point)/|screen point| (E1, E2 the components along R x, R y, R =
  rotation_matrix_from_direction(n), the convention the laser's sigma uses). Far from the source this matches the
  two ratios above; a mismatch means the screen is not in the far field or the field there is not transverse.
- P_z/(c u), expected to be n_z for the same reason.
Ratios are masked where P_z is below 1e-3 of its maximum within the same harmonic window round(omega/omega_1)
(between peaks the field is noise; per window so weak harmonics aren't masked by the strongest).
Two figures: point_spectrum_observables.png (linear) and point_spectrum_observables_log.png (|value| of the four
spectral panels on a log axis, for scans over several harmonics orders of magnitude apart).
The summary (per harmonic window: values at its P_z peak, spin per photon over the window) is printed and
written to the run's run_log.txt under its own section, replaced on a re-run.

The dotted reference line is the laser spin sigma = 2 Im(zeta_1* zeta_2) times sign(n_laser,z), the emitted
photon's spin that plot_observables.py's theory column assumes ("spin stays sigma").

All lab-frame components, so any detector direction works. The observation direction comes from
detector_scatter.dat when the run wrote it (plot_detector_scatter=true), otherwise it is rebuilt from the config
(rectangular/circular detectors only).

Usage: python3 plot_point_spectrum_observables.py [path_to_run_folder]
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plot_point_spectrum import (
    OMEGA_LABEL,
    get_point_spectrum_png_dir,
    mark_harmonics,
    read_radiation_data,
    resolve_radiation_file,
)
from plot_field import read_config_value, read_length, get_screen_coordinates
from plot_observables import (
    EPSILON_0,
    _get_B_total,
    _get_E_total,
    _read_fundamental_frequency_au,
    _read_laser_direction,
    _read_laser_spin,
    _write_run_log_section,
)
from faraday_frame_utils import C_LIGHT, get_detector_local_rotation, rotation_matrix_from_direction

_RATIO_MASK_FRACTION = 1e-3
_LOG_DYNAMIC_RANGE = 1e30
_RUN_LOG_SECTION_HEADER = "Point-spectrum observables (plot_point_spectrum_observables.py)"


def get_observation_directions(radiation_filepath, config_path, n_points):
    """Unit vectors from the origin to each screen point (i_screen order), lab frame."""
    scatter_path = os.path.join(os.path.dirname(radiation_filepath), 'detector_scatter.dat')
    if os.path.exists(scatter_path):
        points = pd.read_csv(scatter_path, sep=' ', comment='#')[['x', 'y', 'z']].to_numpy()
    else:
        detector_type, _ = read_config_value('detector_type', config_path)
        if detector_type not in ('rectangular', 'circular'):
            raise ValueError(f"detector_type '{detector_type}': no detector_scatter.dat in the run folder to read "
                             "the screen points from; re-run with plot_detector_scatter true.")
        x, y, _, _ = get_screen_coordinates(radiation_filepath, config_path)
        distance = read_length(f'{detector_type}_detector_distance', config_path)
        local = np.stack([x, y, np.full_like(x, distance)], axis=1)
        points = local @ get_detector_local_rotation(config_path).T
    if len(points) != n_points:
        raise ValueError(f"{len(points)} detector points but {n_points} screen points in radiation_field.dat")
    return points / np.linalg.norm(points, axis=1, keepdims=True)


def compute_point_observables(radiation_filepath, data):
    """One row per radiation_field.dat row with omega != 0: the local observables and derived ratios."""
    run_dir = os.path.dirname(radiation_filepath)
    config_path = os.path.join(run_dir, 'config.cfg')
    if not os.path.exists(config_path):
        print(f"Error: '{config_path}' not found -- this run predates per-run config snapshots; re-run the solver.")
        sys.exit(1)

    data = data[data['omega'] != 0.0]
    if data.empty:
        print("Error: no frequency with omega/omega_1 != 0 -- nothing to compute.")
        sys.exit(1)

    omega_au = data['omega'].to_numpy() * _read_fundamental_frequency_au(run_dir)
    Ex, Ey, Ez = _get_E_total(data)
    Bx, By, Bz = _get_B_total(data)

    dSz = (4.0 * EPSILON_0 / omega_au) * np.imag(np.conj(Ex) * Ey)
    dSigmazz = (2.0 * EPSILON_0 * C_LIGHT ** 2 / omega_au) * np.imag(
        np.conj(Bx) * Ex + np.conj(By) * Ey - np.conj(Bz) * Ez)
    du = (EPSILON_0 * (np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2)
          + EPSILON_0 * C_LIGHT ** 2 * (np.abs(Bx) ** 2 + np.abs(By) ** 2 + np.abs(Bz) ** 2))
    dPz = 2.0 * EPSILON_0 * C_LIGHT ** 2 * np.real(Ex * np.conj(By) - Ey * np.conj(Bx))

    n_screen = data['i_screen'].nunique()
    directions = get_observation_directions(radiation_filepath, config_path, n_screen)
    n_obs = directions[data['i_screen'].to_numpy()]
    E = np.stack([Ex, Ey, Ez], axis=1)
    E1 = np.empty(len(data), dtype=complex)
    E2 = np.empty(len(data), dtype=complex)
    for i_screen, n in enumerate(directions):
        rows = data['i_screen'].to_numpy() == i_screen
        R = rotation_matrix_from_direction(*n)
        E1[rows] = E[rows] @ R[:, 0]
        E2[rows] = E[rows] @ R[:, 1]
    V = 2.0 * np.imag(np.conj(E1) * E2) / (np.abs(E1) ** 2 + np.abs(E2) ** 2)

    with np.errstate(divide='ignore', invalid='ignore'):
        spin_flux_ratio = omega_au * dSigmazz / dPz
        spin_density_ratio = omega_au * dSz / du
        flux_density_ratio = dPz / (C_LIGHT * du)

    result = pd.DataFrame({
        'omega': data['omega'].to_numpy(), 'i_screen': data['i_screen'].to_numpy(), 'omega_au': omega_au,
        'dSz': dSz, 'dSigmazz': dSigmazz, 'du': du, 'dPz': dPz,
        'spin_flux_ratio': spin_flux_ratio, 'spin_density_ratio': spin_density_ratio,
        'V_nz': V * n_obs[:, 2], 'flux_density_ratio': flux_density_ratio, 'n_z': n_obs[:, 2],
    })

    laser_dir, _ = _read_laser_direction(config_path)
    sigma_ref = _read_laser_spin(config_path) * np.sign(laser_dir[2])
    return result, sigma_ref


def _masked(subset, column):
    """
    `column` with NaN where |P_z| is below _RATIO_MASK_FRACTION of its peak within the same harmonic window
    (round(omega/omega_1) = N), so a scan over several harmonics orders of magnitude apart keeps the weak ones.
    """
    values = subset[column].to_numpy().copy()
    abs_pz = np.abs(subset['dPz'])
    window_peak = abs_pz.groupby(np.round(subset['omega'])).transform('max').to_numpy()
    values[abs_pz.to_numpy() < _RATIO_MASK_FRACTION * window_peak] = np.nan
    return values


def summary_lines(result, sigma_ref):
    """
    Per screen point and harmonic window (round(omega/omega_1) = N): the values at that window's P_z peak and
    the photon-averaged spin over the window.
    """
    lines = [f"  Reference: laser sigma * sign(n_laser,z) = {sigma_ref:+.6g} (spin per photon about +z)."]
    for i_screen, subset in result.groupby('i_screen'):
        subset = subset.sort_values('omega')
        lines.append(f"  screen point {i_screen}: n_z = {subset['n_z'].iloc[0]:.12f}, scan omega/omega_1 = "
                     f"{subset['omega'].min():.6g} .. {subset['omega'].max():.6g} ({len(subset)} points)")
        for N, window in subset.groupby(np.round(subset['omega']).astype(int)):
            peak = window.loc[window['dPz'].abs().idxmax()]
            omega_au = window['omega_au'].to_numpy()
            # Spin per photon over the window: int Sigma_zz domega / int (P_z/omega) domega.
            if len(window) > 1:
                photons = np.trapezoid(window['dPz'].to_numpy() / omega_au, omega_au)
                mean_spin = np.trapezoid(window['dSigmazz'].to_numpy(), omega_au) / photons
            else:
                mean_spin = peak['spin_flux_ratio']
            lines += [
                f"    N = {N} (omega/omega_1 = {window['omega'].min():.6g} .. {window['omega'].max():.6g}, "
                f"{len(window)} points): P_z peak at omega/omega_1 = {peak['omega']:.6g}",
                f"      dP_z/domega = {peak['dPz']:.6e} a.u., du/domega = {peak['du']:.6e} a.u.",
                f"      omega Sigma_zz/P_z = {peak['spin_flux_ratio']:+.10f}",
                f"      omega S_z/u        = {peak['spin_density_ratio']:+.10f}",
                f"      V n_z              = {peak['V_nz']:+.10f}",
                f"      P_z/(c u)          = {peak['flux_density_ratio']:+.10f}",
                f"      spin per photon over the window (int Sigma_zz / int P_z/omega) = {mean_spin:+.10f}",
            ]
    return lines


def draw_observables_figure(result, sigma_ref, png_dir, log_scale):
    """
    The 2x3 observables figure. With `log_scale`, the four spectral panels show |value| on a log axis
    (S_z, Sigma_zz are signed; P_z is negative on a backward screen), clipped _LOG_DYNAMIC_RANGE below the
    largest value; the two ratio panels stay linear.
    """
    n_screen = result['i_screen'].nunique()

    fig, axes = plt.subplots(2, 3, figsize=(17, 9), sharex=True)
    panels = (
        (axes[0, 0], 'dSz', r'd\mathcal{S}_z/d\omega', 'Spin AM density'),
        (axes[0, 1], 'dSigmazz', r'd\Sigma_{zz}/d\omega', 'Spin AM flux'),
        (axes[1, 0], 'du', r'du/d\omega', 'Energy density'),
        (axes[1, 1], 'dPz', r'dP_z/d\omega', 'Energy flux'),
    )
    ax_spin, ax_flux = axes[0, 2], axes[1, 2]
    multi = n_screen > 1

    for i_screen, subset in result.groupby('i_screen'):
        subset = subset.sort_values('omega')
        omega = subset['omega'].to_numpy()
        suffix = f' (point {i_screen})' if multi else ''
        for ax, column, _, _ in panels:
            values = subset[column].to_numpy()
            ax.plot(omega, np.abs(values) if log_scale else values, marker='.', markersize=3,
                    label=f'point {i_screen}')
        ax_spin.plot(omega, _masked(subset, 'spin_flux_ratio'), marker='.', markersize=3,
                     label=r'$\omega\,\Sigma_{zz}/P_z$' + suffix)
        ax_spin.plot(omega, _masked(subset, 'spin_density_ratio'), linestyle='--',
                     label=r'$\omega\,\mathcal{S}_z/u$' + suffix)
        ax_spin.plot(omega, _masked(subset, 'V_nz'), linestyle=':', label=r'$V\,n_z$' + suffix)
        ax_flux.plot(omega, _masked(subset, 'flux_density_ratio'), marker='.', markersize=3,
                     label=r'$P_z/(c\,u)$' + suffix)
        ax_flux.axhline(subset['n_z'].iloc[0], color='gray', linestyle=':', label=r'$n_z$' + suffix)

    ax_spin.axhline(sigma_ref, color='black', linestyle=':', linewidth=1,
                    label=rf'laser $\sigma\,\mathrm{{sign}}(n_{{L,z}}) = {sigma_ref:+.3g}$')

    for ax, column, symbol, name in panels:
        ax.set_title(f"{name} $|{symbol}|$" if log_scale else f"{name} ${symbol}$")
        if log_scale:
            ax.set_yscale('log')
            top = np.abs(result[column]).max()
            if top > 0:
                ax.set_ylim(bottom=max(ax.get_ylim()[0], top / _LOG_DYNAMIC_RANGE))
        if multi:
            ax.legend(fontsize=8)
    # All four ratios lie in [-1, 1]; fixed limits keep matplotlib from switching to an offset axis when they
    # agree to ~1e-8 (the deviations are in the printed summary).
    ax_spin.set_title('Spin per photon about $+z$')
    ax_spin.set_ylim(-1.1, 1.1)
    ax_flux.set_title('Energy flux / density')
    ax_flux.set_ylim(-1.1, 1.1)
    ax_spin.legend(fontsize=8)
    ax_flux.legend(fontsize=8)

    for ax in axes.flat:
        mark_harmonics(ax, result['omega'].to_numpy())
        ax.grid(True, which='major')
    for ax in axes[-1]:
        ax.set_xlabel(OMEGA_LABEL)

    fig.suptitle(f"Local observables vs. {OMEGA_LABEL}{' (log scale)' if log_scale else ''} ({n_screen} screen "
                 f"point{'s' if multi else ''}; ratios masked where $P_z < {_RATIO_MASK_FRACTION:g}\\,P_z^{{max}}$ "
                 "of each harmonic)", fontsize=13, fontweight='bold')
    plt.tight_layout()

    output_img = os.path.join(png_dir, f"point_spectrum_observables{'_log' if log_scale else ''}.png")
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved plot to {output_img}")
    plt.close(fig)


def plot_point_spectrum_observables(radiation_filepath):
    data = read_radiation_data(radiation_filepath)
    result, sigma_ref = compute_point_observables(radiation_filepath, data)
    png_dir = get_point_spectrum_png_dir(radiation_filepath)
    for log_scale in (False, True):
        draw_observables_figure(result, sigma_ref, png_dir, log_scale)

    lines = summary_lines(result, sigma_ref)
    print('\n'.join(lines))
    _write_run_log_section(os.path.dirname(radiation_filepath), _RUN_LOG_SECTION_HEADER, lines)


if __name__ == "__main__":
    input_file = resolve_radiation_file(sys.argv[1] if len(sys.argv) >= 2 else None)
    plot_point_spectrum_observables(input_file)
