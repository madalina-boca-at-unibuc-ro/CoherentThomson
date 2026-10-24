"""
Fine spectrum at a screen point: Re/Im/Abs/Phase of one F^{mu nu} vs omega/omega_1, read from
radiation_field.dat. The counterpart of plot_field.py's plot_radiation_component, which renders one 2D field
map per frequency instead (the wrong plot when frequency, not screen position, is the interesting axis; see the
dense_frequency_spectrum config key, which collapses the run to a single screen point).

Shared with plot_point_spectrum_all.py and plot_point_spectrum_observables.py: the data reader, the per-range
component extraction (get_point_component), the output folder and the harmonic markers.
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

RANGE_TYPES = ('long', 'short', 'boundary', 'total')
_RANGE_PREFIXES = {'long': ('LR',), 'short': ('SR',), 'boundary': ('BR',), 'total': ('LR', 'SR', 'BR')}
OMEGA_LABEL = "$\\omega / \\omega_1$"


def read_radiation_data(radiation_filepath):
    try:
        return pd.read_csv(radiation_filepath, sep=' ', comment='#')
    except Exception as e:
        print(f"Error reading file '{radiation_filepath}': {e}")
        sys.exit(1)


def get_point_component(data, range_type, mu, nu):
    """
    Complex F^{mu nu} of `range_type` (long/short/boundary/total), one value per row of `data`.

    'total' is the physical field LR+SR+BR. BR is identically zero for radiation_formula="direct" (see
    theory/FT_Faraday_tensor-direct_and_simplified_forms-v2.md, Form 2's boundary term), and SR for the
    long_distance_* formulas; in the LD modes only the total is meaningful. BR columns missing from an old run
    count as zero in 'total'.
    """
    value = np.zeros(len(data), dtype=complex)
    for prefix in _RANGE_PREFIXES[range_type]:
        re_col, im_col = f'{prefix}_F{mu}{nu}_re', f'{prefix}_F{mu}{nu}_im'
        if not {re_col, im_col}.issubset(data.columns):
            if range_type == 'total' and prefix == 'BR':
                continue
            print(f"Error: File must contain columns '{re_col}' and '{im_col}'")
            sys.exit(1)
        value += data[re_col].to_numpy() + 1j * data[im_col].to_numpy()
    return value


def get_point_spectrum_png_dir(radiation_filepath):
    """png_folder/radiation/point_spectrum/ of the run, created if needed."""
    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder", MODULE_NAME, "point_spectrum")
    os.makedirs(png_dir, exist_ok=True)
    return png_dir


def mark_harmonics(ax, omega):
    """Dashed vertical lines at the integer harmonics N inside the plotted omega/omega_1 range."""
    for N in range(int(np.ceil(omega.min())), int(np.floor(omega.max())) + 1):
        if N >= 1:
            ax.axvline(N, color='tab:red', linestyle='--', linewidth=0.8, alpha=0.6, zorder=0)


def plot_point_spectrum(range_type, mu, nu, radiation_filepath, data=None):
    """
    Plots Re, Im, |.| and arg (2x2 grid) of F^{mu nu} of `range_type` vs omega/omega_1, one line per distinct
    i_screen (normally one: dense_frequency_spectrum keeps a single screen point). `data` lets a caller that
    plots many components read radiation_field.dat once.
    """
    if data is None:
        data = read_radiation_data(radiation_filepath)
    field = get_point_component(data, range_type, mu, nu)
    png_dir = get_point_spectrum_png_dir(radiation_filepath)

    component_label = f"F^{{{mu}{nu}}}_{{\\mathrm{{{range_type}}}}}"

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ax_re, ax_im, ax_abs, ax_phase = axes.flat

    for i_screen, subset in data.groupby('i_screen'):
        subset = subset.sort_values('omega')
        omega = subset['omega'].to_numpy()
        values = field[subset.index.to_numpy()]

        label = f'screen point {i_screen}'
        ax_re.plot(omega, values.real, marker='.', label=label)
        ax_im.plot(omega, values.imag, marker='.', label=label)
        ax_abs.plot(omega, np.abs(values), marker='.', label=label)
        ax_phase.plot(omega, np.angle(values), marker='.', label=label)

    n_screen = data['i_screen'].nunique()

    ax_re.set_title(f"$\\mathrm{{Re}}({component_label})$")
    ax_im.set_title(f"$\\mathrm{{Im}}({component_label})$")
    ax_abs.set_title(f"$|{component_label}|$")
    ax_phase.set_title(f"$\\arg({component_label})$")
    ax_phase.set_ylim(-np.pi, np.pi)

    for ax in axes.flat:
        mark_harmonics(ax, data['omega'].to_numpy())
        ax.set_xlabel(OMEGA_LABEL)
        ax.grid(True, which='major')
        if n_screen > 1:
            ax.legend(fontsize=8)

    fig.suptitle(f"Radiation spectrum vs. {OMEGA_LABEL} ({n_screen} screen point"
                 f"{'s' if n_screen != 1 else ''})", fontsize=13, fontweight='bold')
    plt.tight_layout()

    output_img = os.path.join(png_dir, f"point_spectrum_{range_type}_F{mu}{nu}.png")
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved plot to {output_img}")
    plt.close(fig)


def resolve_radiation_file(run_folder):
    """radiation_field.dat of `run_folder`, or of the latest run when `run_folder` is None."""
    if run_folder is not None:
        input_file = os.path.join(run_folder, "radiation_field.dat")
        if not os.path.exists(input_file):
            print(f"Error: '{input_file}' not found")
            sys.exit(1)
        return input_file
    try:
        input_file = find_latest_output_file("radiation_field.dat")
    except (ValueError, FileNotFoundError) as e:
        print(f"Usage error: {e}")
        sys.exit(1)
    print(f"No folder given; using latest run: {input_file}")
    return input_file


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: python3 {sys.argv[0]} <long|short|boundary|total> <mu> <nu> [path_to_run_folder]")
        sys.exit(1)

    range_type = sys.argv[1]
    if range_type not in RANGE_TYPES:
        print("Error: first argument must be 'long', 'short', 'boundary', or 'total'")
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

    input_file = resolve_radiation_file(sys.argv[4] if len(sys.argv) >= 5 else None)
    plot_point_spectrum(range_type, mu, nu, input_file)
