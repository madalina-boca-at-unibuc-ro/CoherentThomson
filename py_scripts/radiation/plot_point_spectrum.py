import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

def plot_point_spectrum(range_type, mu, nu, radiation_filepath):
    """
    Plots the real part, imaginary part, magnitude, and phase (2x2 grid) of
    F^{mu nu} of the requested (long_range/short_range/boundary) Faraday tensor as a
    function of omega, read from radiation_field.dat -- the "fine spectrum at
    a point" counterpart of radiation/plot_field.py's plot_radiation_component,
    which instead renders one 2D field-map PNG per frequency (the wrong shape
    of plot when frequency, not screen position, is the interesting axis; see
    the dense_frequency_spectrum config key).

    'boundary' is identically zero when the run used radiation_formula="direct" -- see
    theory/FT_Faraday_tensor-direct_and_simplified_forms.md's "Form 2's boundary term F_b" section.

    Draws one line per distinct i_screen (normally just one, if the detector
    that produced radiation_field.dat was collapsed to a single point per the
    dense_frequency_spectrum convention -- degrades gracefully to an overlaid
    multi-line plot otherwise).
    """
    prefix = {'long': 'LR', 'short': 'SR', 'boundary': 'BR'}[range_type]
    re_col = f'{prefix}_F{mu}{nu}_re'
    im_col = f'{prefix}_F{mu}{nu}_im'

    try:
        data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    except Exception as e:
        print(f"Error reading file '{radiation_filepath}': {e}")
        sys.exit(1)

    if not {re_col, im_col}.issubset(data.columns):
        print(f"Error: File must contain columns '{re_col}' and '{im_col}'")
        sys.exit(1)

    # Per-module subfolder of the run directory's 'png_folder' (mirroring py_scripts/'s own
    # laser/detector/particle/radiation/debug layout).
    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder", MODULE_NAME)
    os.makedirs(png_dir, exist_ok=True)

    component_label = f"F^{{{mu}{nu}}}_{{\\mathrm{{{range_type}}}}}"

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ax_re, ax_im, ax_abs, ax_phase = axes.flat

    for i_screen, subset in data.groupby('i_screen'):
        subset = subset.sort_values('omega')
        omega = subset['omega'].to_numpy()
        re_values = subset[re_col].to_numpy()
        im_values = subset[im_col].to_numpy()
        abs_values = np.hypot(re_values, im_values)
        phase_values = np.arctan2(im_values, re_values)

        label = f'screen point {i_screen}'
        ax_re.plot(omega, re_values, marker='.', label=label)
        ax_im.plot(omega, im_values, marker='.', label=label)
        ax_abs.plot(omega, abs_values, marker='.', label=label)
        ax_phase.plot(omega, phase_values, marker='.', label=label)

    n_screen = data['i_screen'].nunique()

    ax_re.set_title(f"$\\mathrm{{Re}}({component_label})$")
    ax_im.set_title(f"$\\mathrm{{Im}}({component_label})$")
    ax_abs.set_title(f"$|{component_label}|$")
    ax_phase.set_title(f"$\\arg({component_label})$")
    ax_phase.set_ylim(-np.pi, np.pi)

    for ax in axes.flat:
        ax.set_xlabel("$\\omega / \\omega_1$ (units of the fundamental)")
        ax.grid(True, which='major')
        if n_screen > 1:
            ax.legend(fontsize=8)

    fig.suptitle(f"Radiation spectrum vs. $\\omega/\\omega_1$ ({n_screen} screen point"
                 f"{'s' if n_screen != 1 else ''})", fontsize=13, fontweight='bold')
    plt.tight_layout()

    output_name = f"point_spectrum_{range_type}_F{mu}{nu}.png"
    output_img = os.path.join(png_dir, output_name)
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved plot to {output_img}")
    plt.close(fig)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: python3 {sys.argv[0]} <long|short|boundary> <mu> <nu> [path_to_run_folder]")
        sys.exit(1)

    range_type = sys.argv[1]
    if range_type not in ('long', 'short', 'boundary'):
        print("Error: first argument must be 'long', 'short', or 'boundary'")
        sys.exit(1)

    try:
        mu = int(sys.argv[2])
        nu = int(sys.argv[3])
    except ValueError:
        print("Error: mu and nu must be integers")
        sys.exit(1)

    if not (0 <= mu <= 3 and 0 <= nu <= 3):
        print("Error: mu and nu must each be in [0, 3]")
        sys.exit(1)

    if len(sys.argv) >= 5:
        input_file = os.path.join(sys.argv[4], "radiation_field.dat")
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

    plot_point_spectrum(range_type, mu, nu, input_file)
