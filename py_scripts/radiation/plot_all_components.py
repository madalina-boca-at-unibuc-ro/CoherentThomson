"""
Convenience driver over plot_field.py's plot_radiation_component: calls it for every
(range_type, mu, nu) combination instead of requiring one invocation per combination -- 4 range
types ('long', 'short', 'boundary', 'total') x the 6 independent upper-triangle (mu < nu) Faraday
tensor components (F^{mu nu} = -F^{nu mu}, diagonal identically zero), 24 combinations total, each
producing one PNG per configured frequency in png_folder/radiation/ (same output plot_field.py's own
CLI produces one combination at a time).

Angular-momentum flux (radiation/plot_angular_momentum_flux.py) is a separate, detector-restricted
post-processing step and is intentionally not looped over here.
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file
from plot_field import plot_radiation_component

RANGE_TYPES = ('long', 'short', 'boundary', 'total')
COMPONENT_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

if __name__ == "__main__":
    if len(sys.argv) >= 2:
        input_file = os.path.join(sys.argv[1], "radiation_field.dat")
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

    for range_type in RANGE_TYPES:
        for mu, nu in COMPONENT_PAIRS:
            print(f"Plotting {range_type} F^{{{mu}{nu}}}...")
            plot_radiation_component(range_type, mu, nu, input_file)
