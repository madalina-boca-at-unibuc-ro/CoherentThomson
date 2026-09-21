"""
Shared helpers for locating a run's output files.

Every plotting script needs to find "the most recent run directory" that
Core::IoUtils::make_run_output_directory created under the .cfg file's
'output_folder' key (a '<output_folder>/YYYYMMDD_HHMMSS/' subfolder). This
module centralizes that lookup so it isn't duplicated per script.
"""

import math
import os
import re

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DEFAULT_CONFIG_PATH = os.path.join(PROJECT_ROOT, "config", "config.cfg")

RUN_DIR_PATTERN = re.compile(r"^\d{8}_\d{6}$")

def get_output_folder(config_path=DEFAULT_CONFIG_PATH):
    """
    Reads the 'output_folder' key out of the .cfg file (the same key
    Core::IoUtils::make_run_output_directory reads) and resolves it relative to
    the home directory, mirroring Core::IoUtils::make_run_output_directory:
    the config value is always a subfolder of '~', not an absolute path, so the
    same .cfg still works after sync_to_remote.sh copies the repo to a
    different machine/user.
    """
    with open(config_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            parts = stripped.split()
            if parts[0] == 'output_folder' and len(parts) >= 2:
                return os.path.join(os.path.expanduser('~'), parts[1])
    raise ValueError(f"No 'output_folder' key found in '{config_path}'")

def find_latest_run_dir(output_folder):
    """
    Finds the most recent '<output_folder>/YYYYMMDD_HHMMSS/' run directory
    (created by IoUtils::make_run_output_directory) and returns its path.
    """
    if not os.path.isdir(output_folder):
        raise FileNotFoundError(f"Output folder '{output_folder}' does not exist")

    run_dirs = [d for d in os.listdir(output_folder) if RUN_DIR_PATTERN.match(d)]
    if not run_dirs:
        raise FileNotFoundError(f"No run directories found under '{output_folder}'")

    return os.path.join(output_folder, max(run_dirs))

def find_latest_output_file(filename, output_folder=None):
    """
    Returns the path to 'filename' inside the most recent run directory under
    output_folder (read from the .cfg file if not given explicitly).
    """
    if output_folder is None:
        output_folder = get_output_folder()

    latest_run_dir = find_latest_run_dir(output_folder)
    filepath = os.path.join(latest_run_dir, filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"'{filepath}' does not exist")
    return filepath

def read_config_value(config_path, key):
    """
    Reads a single 'key value [unit]' entry out of a .cfg file (the raw value token only, unit
    suffix if any left unparsed) -- mirrors Core::IoUtils::get_required's key lookup. Meant to be
    called with a specific run's own 'config.cfg' (written by
    Core::IoUtils::copy_config_to_run_directory into every run directory), so a plot reads back the
    config that actually produced that run's data rather than whatever config/coherent_thomson.cfg
    currently contains.
    """
    with open(config_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            parts = stripped.split()
            if parts[0] == key and len(parts) >= 2:
                return parts[1]
    raise ValueError(f"No '{key}' key found in '{config_path}'")

def get_laser_period(run_dir):
    """
    Returns the laser period T = 2*pi/omega, reading 'laser_frequency' (bare omega, no unit
    suffix -- see Core::IoUtils::get_laser_frequency) from run_dir's own 'config.cfg'.
    """
    config_path = os.path.join(run_dir, "config.cfg")
    omega = float(read_config_value(config_path, "laser_frequency"))
    return 2 * math.pi / omega

# Core::PhysUtils::AtomicUnits::c (phys_utils.hpp): c in the solver's own atomic units (1/alpha,
# 2018 CODATA); m_0 = 1.0 there too, so mc is just this constant. Mirrored independently (no shared
# constants module between C++ and Python) -- keep in sync if that value changes.
_C_LIGHT = 137.035999084

def _read_config_value_and_unit(config_path, key):
    """Like read_config_value, but also returns the optional unit suffix token ('' if none)."""
    with open(config_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            parts = stripped.split()
            if parts[0] == key and len(parts) >= 2:
                return parts[1], (parts[2] if len(parts) >= 3 else '')
    raise ValueError(f"No '{key}' key found in '{config_path}'")

def _momentum_component_au(config_path, key):
    """
    Reads one average_p{x,y,z} component and converts it to atomic units. Only the 'mc' unit
    (Core::IoUtils::convert_unit_to_number's 'mc' branch, = m_0*c = _C_LIGHT since m_0=1) and bare/
    'a.u.' values are handled -- the only units this repo's own configs ever give these keys in.
    """
    value_str, unit = _read_config_value_and_unit(config_path, key)
    value = float(value_str)
    if unit in ('', 'a.u.'):
        return value
    if unit.lower() == 'mc':
        return value * _C_LIGHT
    raise ValueError(f"Unsupported unit '{unit}' for '{key}' in '{config_path}'")

def get_doppler_shifted_period(run_dir):
    """
    Returns the proper-time period T_D = T/D corresponding to one laser cycle as actually seen by an
    electron with the beam's average initial four-momentum, where D = (p^0-p^3)/mc is exactly the
    same Doppler factor Core::Simulation::init_simulation_parameters (simulation.cpp) already divides
    the RK4 trajectory's time step d_tau_traj by, to keep NT_traj points per cycle resolving the field
    as seen in the electron's own proper time (not the bare lab-frame T) whenever average_px/py/pz
    isn't zero. Reads average_px/py/pz back out of run_dir's own 'config.cfg' snapshot, the same way
    get_laser_period reads back laser_frequency.
    """
    config_path = os.path.join(run_dir, "config.cfg")
    T = get_laser_period(run_dir)
    px = _momentum_component_au(config_path, "average_px")
    py = _momentum_component_au(config_path, "average_py")
    pz = _momentum_component_au(config_path, "average_pz")
    mc = _C_LIGHT
    p0 = math.sqrt(px * px + py * py + pz * pz + mc * mc)
    doppler_factor = (p0 - pz) / mc
    return T / doppler_factor
