"""
Field components of a run's radiation_field.dat (theory/specs/13_helicity_plots.md §3): the total (LR + SR + BR)
E and B, the helicity components E_+- (lab or local transverse basis), and the laser's spin.

The only implementation of each: plot_observables.py and lg_decomposition.py re-export the old names
(_get_total_field, _get_E_total, _get_B_total, _read_laser_spin; helicity_components).

Imports: plot_field.py imports this module only inside functions, because faraday_frame_utils (imported here)
imports plot_field at module level.
"""
import math

import numpy as np

from plot_field import read_config_value
from faraday_frame_utils import C_LIGHT, rotation_matrix_from_direction


def _complex_column(data, prefix, mu, nu):
    return (data[f'{prefix}_F{mu}{nu}_re'] + 1j * data[f'{prefix}_F{mu}{nu}_im']).to_numpy()


def get_total_field(data, index_pairs, scale=1.0):
    """
    Generic total-(LR+SR+BR) field extractor: one complex numpy array per (mu, nu) pair in
    `index_pairs`, each scaled by `scale`. Shared by get_E_total (E_i=c*F^{i0}, scale=C_LIGHT) and
    get_B_total (B_i=F^{jk} directly, scale=1). Only the total is physical: the long/short/boundary split
    differs between radiation formulas. BR columns are treated as zero if the file predates the boundary term.
    """
    def component(mu, nu):
        value = _complex_column(data, 'LR', mu, nu) + _complex_column(data, 'SR', mu, nu)
        if {f'BR_F{mu}{nu}_re', f'BR_F{mu}{nu}_im'}.issubset(data.columns):
            value = value + _complex_column(data, 'BR', mu, nu)
        return scale * value

    return tuple(component(mu, nu) for mu, nu in index_pairs)


def get_E_total(data):
    """Returns (Ex, Ey, Ez) from the total Faraday tensor. E_i = c*F^{i0}."""
    return get_total_field(data, ((1, 0), (2, 0), (3, 0)), scale=C_LIGHT)


def get_B_total(data):
    """Returns (Bx, By, Bz) from the total Faraday tensor. Bx=F^{32}, By=F^{13}, Bz=F^{21}."""
    return get_total_field(data, ((3, 2), (1, 3), (2, 1)), scale=1.0)


def helicity_components(Ex, Ey):
    """(E_+, E_-): amplitudes on e_+- = (x + i y)/sqrt2 (lab x, y), i.e. spin +-1 about lab +z for
    exp(-i omega t) -- S_z ~ Im(Ex* Ey) = |E_+|^2 - |E_-|^2."""
    return (Ex - 1j * Ey) / math.sqrt(2.0), (Ex + 1j * Ey) / math.sqrt(2.0)


def local_helicity_components(Ex, Ey, Ez, directions):
    """
    (E_+, E_-) in the transverse basis of each point's own direction n (spec 13 §2, --local-basis):
    E_1 = E . R x, E_2 = E . R y with R = rotation_matrix_from_direction(n), then E_+- = (E_1 -+ i E_2)/sqrt2.
    The projection plot_point_spectrum_observables.py uses. `directions` (one unit vector per element of Ex)
    comes from plot_point_spectrum_observables.get_observation_directions, indexed by i_screen.
    Agrees with helicity_components to O(theta^2).
    """
    E = np.stack([Ex, Ey, Ez], axis=1)
    E1 = np.empty(len(E), dtype=complex)
    E2 = np.empty(len(E), dtype=complex)
    unique, inverse = np.unique(np.asarray(directions), axis=0, return_inverse=True)
    for k, n in enumerate(unique):
        rows = inverse.ravel() == k
        R = rotation_matrix_from_direction(*n)
        E1[rows] = E[rows] @ R[:, 0]
        E2[rows] = E[rows] @ R[:, 1]
    return helicity_components(E1, E2)


def read_laser_spin(config_path):
    """
    The laser's spin per photon along its own propagation direction, sigma = 2 Im(conj(zeta_1) zeta_2) with
    (zeta_1, zeta_2) normalized as create_laser does: +1 for the default circular zeta_1=(1,0), zeta_2=(0,1)
    (checked against the run's own S_z/P_z), 0 for linear, in between for elliptical.
    """
    zeta = [complex(float(read_config_value(f'laser_zeta_{i}_re', config_path)[0]),
                    float(read_config_value(f'laser_zeta_{i}_im', config_path)[0])) for i in (1, 2)]
    norm_sq = abs(zeta[0]) ** 2 + abs(zeta[1]) ** 2
    return 2.0 * (np.conj(zeta[0]) * zeta[1]).imag / norm_sq
