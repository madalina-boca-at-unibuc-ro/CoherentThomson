"""
Shared frame/rotation helpers for turning radiation_field.dat's stored Faraday tensor into
Cartesian E/B components projected onto a detector's own local direction. Used by
radiation/plot_spherical_components.py to project the field onto spherical components.
"""
import numpy as np

from plot_field import read_config_value

# Core::PhysUtils::AtomicUnits (phys_utils.hpp): c in the solver's own atomic units -- c is 1/alpha
# (2018 CODATA). Mirrors phys_utils.hpp's own value independently (no shared constants module
# between C++ and Python) -- keep in sync if that one changes.
C_LIGHT = 137.035999084


def convert_unit_to_number(unit, config_path):
    """Python mirror of Core::IoUtils::convert_unit_to_number (io_utils.hpp) -- only the branches
    detector-geometry/direction keys actually use ('lambda', 'pi', 'a.u.', bare)."""
    unit = unit.lower()
    if unit == 'pi':
        return np.pi
    if unit == 'lambda':
        omega = float(read_config_value('laser_frequency', config_path)[0])
        return 2.0 * np.pi / omega * C_LIGHT
    if unit == 'mc':
        return C_LIGHT
    if unit == 'cycles_adim':
        return 2.0 * np.pi
    if unit == 'omega_laser':
        return float(read_config_value('laser_frequency', config_path)[0])
    if unit in ('a.u.', ''):
        return 1.0
    print(f"Warning: Unknown unit '{unit}'. Assuming 1.0.")
    return 1.0


def rotation_matrix_from_direction(nx, ny, nz):
    """
    Python mirror of Core::MathUtils::rotation_matrix_from_direction (math_utils.hpp): the 3x3
    rotation mapping local +z onto the given direction, including its two hardcoded near-degenerate
    cases (exactly +z / exactly -z).
    """
    norm = np.sqrt(nx * nx + ny * ny + nz * nz)
    if norm < 1e-12:
        return np.eye(3)
    nx, ny, nz = nx / norm, ny / norm, nz / norm
    if nz > 0.999999:
        return np.eye(3)
    if nz < -0.999999:
        return np.diag([1.0, -1.0, -1.0])
    k = 1.0 / (1.0 + nz)
    return np.array([
        [1.0 - nx * nx * k, -nx * ny * k, nx],
        [-nx * ny * k, 1.0 - ny * ny * k, ny],
        [-nx, -ny, nz],
    ])


def get_detector_local_rotation(config_path):
    """
    Returns the 3x3 rotation R such that v_canonical = R @ v_local -- a Python mirror of
    Core::Detector::Detector_2D's own local_rotation, built from detector_direction_theta/phi.
    """
    theta_str, theta_unit = read_config_value('detector_direction_theta', config_path)
    phi_str, phi_unit = read_config_value('detector_direction_phi', config_path)
    theta = float(theta_str) * convert_unit_to_number(theta_unit, config_path)
    phi = float(phi_str) * convert_unit_to_number(phi_unit, config_path)
    dir_x, dir_y, dir_z = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)
    return rotation_matrix_from_direction(dir_x, dir_y, dir_z)


def extract_rotated_faraday_fields(data):
    """
    Returns a dict of six complex numpy arrays per range (Ex/Ey/Ez/Bx/By/Bz, suffixed '_l' for
    long-range, '_s' for short-range, '_b' for boundary), read directly out of radiation_field.dat's
    F^{mu nu} columns (the laser always propagates along canonical Oz, so no rotation is needed).
    '_b' is identically zero when the run used radiation_formula="direct" -- see
    theory/FT_Faraday_tensor-direct_and_simplified_forms.md's "Form 2's boundary term F_b" section.

    Column mapping mirrors Core::Laser::LaserField::get_faraday_tensor's F^{mu nu} <-> E/B sign
    convention (F^{i0}=E_i, F^{jk}=-eps_jkl B_l): Ex=F10, Ey=F20, Ez=F30, Bx=F32, By=F13, Bz=F21.
    """
    fields = {}
    for prefix, suffix in (('LR', 'l'), ('SR', 's'), ('BR', 'b')):
        def col(mu, nu):
            return (data[f'{prefix}_F{mu}{nu}_re'] + 1j * data[f'{prefix}_F{mu}{nu}_im']).to_numpy()

        Ex, Ey, Ez = col(1, 0), col(2, 0), col(3, 0)
        Bx, By, Bz = col(3, 2), col(1, 3), col(2, 1)
        fields[f'Ex_{suffix}'], fields[f'Ey_{suffix}'], fields[f'Ez_{suffix}'] = Ex, Ey, Ez
        fields[f'Bx_{suffix}'], fields[f'By_{suffix}'], fields[f'Bz_{suffix}'] = Bx, By, Bz
    return fields
