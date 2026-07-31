import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file

def plot_point_spectrum(range_type, mu, nu, radiation_filepath):
    """
    Plots the real part, imaginary part, magnitude, and phase (2x2 grid) of
    F^{mu nu} of the requested (long_range/short_range) Faraday tensor as a
    function of omega, read from radiation_field.dat -- the "fine spectrum at
    a point" counterpart of plot_radiation_field.py's plot_radiation_component,
    which instead renders one 2D field-map PNG per frequency (the wrong shape
    of plot when frequency, not screen position, is the interesting axis; see
    the dense_frequency_spectrum config key).

    Draws one line per distinct i_screen (normally just one, if the detector
    that produced radiation_field.dat was collapsed to a single point per the
    dense_frequency_spectrum convention -- degrades gracefully to an overlaid
    multi-line plot otherwise).
    """
    prefix = {'long': 'LR', 'short': 'SR'}[range_type]
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

    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder")
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
        ax.set_xlabel("$\\omega$")
        if n_screen > 1:
            ax.legend(fontsize=8)

    fig.suptitle(f"Radiation spectrum vs. $\\omega$ ({n_screen} screen point"
                 f"{'s' if n_screen != 1 else ''})", fontsize=13, fontweight='bold')
    plt.tight_layout()

    output_name = f"point_spectrum_{range_type}_F{mu}{nu}.png"
    output_img = os.path.join(png_dir, output_name)
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved plot to {output_img}")
    plt.close(fig)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: python3 {sys.argv[0]} <long|short> <mu> <nu> [path_to_radiation_field.dat]")
        sys.exit(1)

    range_type = sys.argv[1]
    if range_type not in ('long', 'short'):
        print("Error: first argument must be 'long' or 'short'")
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
        input_file = sys.argv[4]
    else:
        try:
            input_file = find_latest_output_file("radiation_field.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            sys.exit(1)
        print(f"No file given; using latest run's radiation field: {input_file}")

    plot_point_spectrum(range_type, mu, nu, input_file)
