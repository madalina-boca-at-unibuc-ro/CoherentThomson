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
Ratios are masked where P_z is below 1e-3 of its maximum (between harmonics the field is noise).

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
)
from faraday_frame_utils import C_LIGHT, get_detector_local_rotation, rotation_matrix_from_direction

_RATIO_MASK_FRACTION = 1e-3


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
    values = subset[column].to_numpy().copy()
    values[np.abs(subset['dPz'].to_numpy()) < _RATIO_MASK_FRACTION * np.abs(subset['dPz']).max()] = np.nan
    return values


def print_summary(result):
    """Per screen point: the values at the P_z peak and the photon-averaged spin over the whole scan."""
    for i_screen, subset in result.groupby('i_screen'):
        subset = subset.sort_values('omega')
        peak = subset.loc[subset['dPz'].abs().idxmax()]
        omega_au = subset['omega_au'].to_numpy()
        # Spin per photon over the scan: int Sigma_zz domega / int (P_z/omega) domega.
        photons = np.trapezoid(subset['dPz'].to_numpy() / omega_au, omega_au)
        mean_spin = np.trapezoid(subset['dSigmazz'].to_numpy(), omega_au) / photons
        print(f"screen point {i_screen} (n_z = {peak['n_z']:.9f}):")
        print(f"  P_z peak at omega/omega_1 = {peak['omega']:.6g}: omega Sigma_zz/P_z = {peak['spin_flux_ratio']:+.6g}, "
              f"omega S_z/u = {peak['spin_density_ratio']:+.6g}, V n_z = {peak['V_nz']:+.6g}, "
              f"P_z/(c u) = {peak['flux_density_ratio']:.6g}")
        print(f"  spin per photon over the scan (int Sigma_zz / int P_z/omega) = {mean_spin:+.6g}")


def plot_point_spectrum_observables(radiation_filepath):
    data = read_radiation_data(radiation_filepath)
    result, sigma_ref = compute_point_observables(radiation_filepath, data)
    png_dir = get_point_spectrum_png_dir(radiation_filepath)
    n_screen = result['i_screen'].nunique()

    fig, axes = plt.subplots(2, 3, figsize=(17, 9), sharex=True)
    panels = (
        (axes[0, 0], 'dSz', r'Spin AM density $d\mathcal{S}_z/d\omega$'),
        (axes[0, 1], 'dSigmazz', r'Spin AM flux $d\Sigma_{zz}/d\omega$'),
        (axes[1, 0], 'du', r'Energy density $du/d\omega$'),
        (axes[1, 1], 'dPz', r'Energy flux $dP_z/d\omega$'),
    )
    ax_spin, ax_flux = axes[0, 2], axes[1, 2]
    multi = n_screen > 1

    for i_screen, subset in result.groupby('i_screen'):
        subset = subset.sort_values('omega')
        omega = subset['omega'].to_numpy()
        suffix = f' (point {i_screen})' if multi else ''
        for ax, column, _ in panels:
            ax.plot(omega, subset[column].to_numpy(), marker='.', markersize=3, label=f'point {i_screen}')
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

    for ax, _, title in panels:
        ax.set_title(title)
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

    fig.suptitle(f"Local observables vs. {OMEGA_LABEL} ({n_screen} screen point{'s' if multi else ''}; "
                 f"ratios masked where $P_z < {_RATIO_MASK_FRACTION:g}\\,P_z^{{max}}$)",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()

    output_img = os.path.join(png_dir, "point_spectrum_observables.png")
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved plot to {output_img}")
    plt.close(fig)

    print_summary(result)


if __name__ == "__main__":
    input_file = resolve_radiation_file(sys.argv[1] if len(sys.argv) >= 2 else None)
    plot_point_spectrum_observables(input_file)
