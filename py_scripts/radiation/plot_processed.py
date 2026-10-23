"""
Post-processing of the emitted field (theory/post_processing.md): plots a processed version of
radiation_field.dat's F^{mu nu} in plot_field.py's format, into png_folder/radiation/processed/<processing>/.

Processings:
- radial_phase: F^{mu nu} * exp(-i k_N (R - Z)), removing the spherical-wave phase of a source at the
  origin. k_N = omega_N/c is each frequency row's own wavevector, Z the screen distance and
  R = sqrt(X^2 + Y^2 + Z^2) with X, Y the screen coordinates (each screen point sits at local (X, Y, Z),
  rotated with the detector, so this holds for a tilted detector too). The minus sign matches the
  solver's exp(+i omega/c (x^0 + R)) convention (radiation.cpp's radiation_phase_argument), in which
  the emitted field goes as exp(+i k R). Rectangular and circular detectors only.
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file
from plot_field import plot_radiation_component, read_config_value, read_length, get_screen_coordinates
from plot_observables import _read_fundamental_frequency_au
from plot_all_components import RANGE_TYPES, COMPONENT_PAIRS

PROCESSINGS = ('radial_phase',)


def make_radial_phase_transform(radiation_filepath, config_path):
    """Returns field_transform(field, omega_ratio) multiplying by exp(-i k_N (R - Z)). Lengths are in
    lambda, so k_N = 2*pi*omega_N/omega_laser per lambda, with omega_N = omega_ratio * omega_1."""
    detector_type, _ = read_config_value('detector_type', config_path)
    if detector_type not in ('rectangular', 'circular'):
        raise ValueError(f"radial_phase needs a flat (rectangular/circular) detector, not '{detector_type}': "
                         "a spherical screen has constant R and no single screen distance Z")

    x, y, _, _ = get_screen_coordinates(radiation_filepath, config_path)
    Z = read_length(f'{detector_type}_detector_distance', config_path, target_unit='lambda')
    rho_sq = x**2 + y**2
    # R - Z = rho^2/(R + Z): direct subtraction loses digits on a far screen, where R ~ Z.
    R_minus_Z = rho_sq / (np.sqrt(rho_sq + Z**2) + Z)

    omega_1_over_omega_laser = (_read_fundamental_frequency_au(os.path.dirname(radiation_filepath))
                                / float(read_config_value('laser_frequency', config_path)[0]))

    def transform(field, omega_ratio):
        k_N = 2.0 * np.pi * omega_ratio * omega_1_over_omega_laser
        return field * np.exp(-1j * k_N * R_minus_Z)

    return transform


TRANSFORM_FACTORIES = {'radial_phase': make_radial_phase_transform}
TITLE_NOTES = {'radial_phase': r'\, e^{-ik_N(R-Z)}'}


if __name__ == "__main__":
    usage = (f"Usage: python3 {sys.argv[0]} <{'|'.join(PROCESSINGS)}> <long|short|boundary|total> <mu> <nu> "
             f"[run_folder]\n       python3 {sys.argv[0]} <{'|'.join(PROCESSINGS)}> all [run_folder]")
    args = sys.argv[1:]
    if len(args) < 2 or args[0] not in PROCESSINGS:
        print(usage)
        sys.exit(1)
    processing = args[0]

    if args[1] == 'all':
        combinations = [(range_type, mu, nu) for range_type in RANGE_TYPES for mu, nu in COMPONENT_PAIRS]
        rest = args[2:]
    else:
        if len(args) < 4 or args[1] not in RANGE_TYPES:
            print(usage)
            sys.exit(1)
        try:
            mu, nu = int(args[2]), int(args[3])
        except ValueError:
            print("Error: mu and nu must be integers")
            sys.exit(1)
        if not (0 <= mu <= 3 and 0 <= nu <= 3) or mu == nu:
            print("Error: mu and nu must each be in [0, 3] and differ")
            sys.exit(1)
        combinations = [(args[1], mu, nu)]
        rest = args[4:]

    if rest:
        input_file = os.path.join(rest[0], "radiation_field.dat")
        if not os.path.exists(input_file):
            print(f"Error: '{input_file}' not found")
            sys.exit(1)
    else:
        try:
            input_file = find_latest_output_file("radiation_field.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            sys.exit(1)
        print(f"No folder given; using latest run: {input_file}")

    config_path = os.path.join(os.path.dirname(input_file), 'config.cfg')
    try:
        transform = TRANSFORM_FACTORIES[processing](input_file, config_path)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    for range_type, mu, nu in combinations:
        print(f"Plotting {processing} {range_type} F^{{{mu}{nu}}}...")
        plot_radiation_component(range_type, mu, nu, input_file, field_transform=transform,
                                 png_subdir=os.path.join("processed", processing), title_note=TITLE_NOTES[processing])
