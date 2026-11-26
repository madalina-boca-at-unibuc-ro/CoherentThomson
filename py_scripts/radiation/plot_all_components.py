"""
Convenience driver over plot_field.py's plot_radiation_component: calls it for every
(range_type, mu, nu) combination instead of requiring one invocation per combination -- 4 range
types ('long', 'short', 'boundary', 'total') x the 6 independent upper-triangle (mu < nu) Faraday
tensor components (F^{mu nu} = -F^{nu mu}, diagonal identically zero), 24 combinations total, each
producing one PNG per configured frequency in png_folder/radiation/emitted/ (same output plot_field.py's
own CLI produces one combination at a time).

Accepts --incident: plots incident_field.dat (Core::Radiation::export_incident_field_fourier's
known-analytic incident-beam reference field) instead of radiation_field.dat.
plot_radiation_component itself routes those PNGs into png_folder/radiation/incident/ instead (see
its own doc comment), so they never collide with the real run's own plots even when looping over
every combination like this.

After the 24 F plots, the two helicity components E_+- (plot_field.plot_helicity_component, spec 13), raw and, for
the emitted field, with the spherical carrier removed (plot_field_carrier_removed.py plus|minus).
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file
from plot_field import plot_radiation_component, plot_helicity_component

RANGE_TYPES = ('long', 'short', 'boundary', 'total')
COMPONENT_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

if __name__ == "__main__":
    args = sys.argv[1:]
    use_incident = "--incident" in args
    if use_incident:
        args.remove("--incident")
    field_filename = "incident_field.dat" if use_incident else "radiation_field.dat"

    if len(args) >= 1:
        input_file = os.path.join(args[0], field_filename)
        if not os.path.exists(input_file):
            print(f"Error: '{input_file}' not found")
            sys.exit(1)
    else:
        try:
            input_file = find_latest_output_file(field_filename)
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            sys.exit(1)
        print(f"No folder given; using latest run: {input_file}")

    for range_type in RANGE_TYPES:
        for mu, nu in COMPONENT_PAIRS:
            print(f"Plotting {range_type} F^{{{mu}{nu}}}...")
            plot_radiation_component(range_type, mu, nu, input_file)

    for sign in ('plus', 'minus'):
        print(f"Plotting E_{sign}...")
        plot_helicity_component(sign, input_file)
    if not use_incident:
        # Imported here: plot_field_carrier_removed imports this module.
        from plot_field_carrier_removed import make_radial_phase_transform, PNG_SUBDIR, TITLE_NOTE
        transform = make_radial_phase_transform(input_file, os.path.join(os.path.dirname(input_file), 'config.cfg'))
        for sign in ('plus', 'minus'):
            print(f"Plotting E_{sign} with the carrier removed...")
            plot_helicity_component(sign, input_file, field_transform=transform, png_subdir=PNG_SUBDIR,
                                    title_note=TITLE_NOTE)
