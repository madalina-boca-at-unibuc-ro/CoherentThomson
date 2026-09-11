import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file

def plot_debug_integrand(range_type, mu, nu, debug_filepath):
    """
    Plots the real part, imaginary part, magnitude, and phase (2x2 grid) of the requested
    (long_range/short_range/boundary) Faraday bivector component's per-tau integrand, read from
    debug_integrand.dat (Core::Debug::export_radiation_integrand, only written when the config's
    'debug' key is true) -- the raw per-trajectory-point terms Radiation::compute_radiation sums
    over tau to build the coherent field, for the single electron/screen point/wavenumber k (= omega/c)
    that debug mode requires. The "point spectrum" counterpart (plot_point_spectrum.py) instead plots
    the already-summed field vs. omega; this plots the unsummed tau-integrand at one k.

    'boundary' is nonzero only at the first/last tau (everywhere else it's exactly zero by
    construction -- it isn't itself a sum over tau, unlike long/short -- see
    theory/FT_Faraday_tensor-direct_and_simplified_forms.md's "Form 2's boundary term F_b" section),
    so this plot shows two spikes rather than a smooth curve.

    debug_integrand.dat only stores the 6 independent upper-triangle (mu < nu) bivector elements
    (F^{mu nu} = -F^{nu mu}, diagonal zero), so mu == nu is rejected and mu > nu is resolved by
    reading the (nu, mu) column and negating it.
    """
    if mu == nu:
        print(f"Error: mu and nu must differ (F^{{{mu}{nu}}} is identically zero)")
        sys.exit(1)
    sign = 1.0 if mu < nu else -1.0
    a, b = (mu, nu) if mu < nu else (nu, mu)

    prefix = {'long': 'LR', 'short': 'SR', 'boundary': 'BR'}[range_type]
    re_col = f'{prefix}_F{a}{b}_re'
    im_col = f'{prefix}_F{a}{b}_im'

    try:
        data = pd.read_csv(debug_filepath, sep=' ', comment='#')
    except Exception as e:
        print(f"Error reading file '{debug_filepath}': {e}")
        sys.exit(1)

    if not {re_col, im_col}.issubset(data.columns):
        print(f"Error: File must contain columns '{re_col}' and '{im_col}'")
        sys.exit(1)

    png_dir = os.path.join(os.path.dirname(debug_filepath), "png_folder")
    os.makedirs(png_dir, exist_ok=True)

    data = data.sort_values('tau')
    tau = data['tau'].to_numpy()
    re_values = sign * data[re_col].to_numpy()
    im_values = sign * data[im_col].to_numpy()
    abs_values = np.hypot(re_values, im_values)
    phase_values = np.arctan2(im_values, re_values)

    component_label = f"F^{{{mu}{nu}}}_{{\\mathrm{{{range_type}}}}}"

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ax_re, ax_im, ax_abs, ax_phase = axes.flat

    ax_re.plot(tau, re_values, marker='.')
    ax_im.plot(tau, im_values, marker='.')
    ax_abs.plot(tau, abs_values, marker='.')
    ax_phase.plot(tau, phase_values, marker='.')

    ax_re.set_title(f"$\\mathrm{{Re}}({component_label})$")
    ax_im.set_title(f"$\\mathrm{{Im}}({component_label})$")
    ax_abs.set_title(f"$|{component_label}|$")
    ax_phase.set_title(f"$\\arg({component_label})$")
    ax_phase.set_ylim(-np.pi, np.pi)

    for ax in axes.flat:
        ax.set_xlabel(r"$\tau$")

    fig.suptitle(f"Radiation integrand vs. $\\tau$ ({len(tau)} trajectory points)", fontsize=13, fontweight='bold')
    plt.tight_layout()

    output_name = f"debug_integrand_{range_type}_F{mu}{nu}.png"
    output_img = os.path.join(png_dir, output_name)
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved plot to {output_img}")
    plt.show()
    plt.close(fig)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: python3 {sys.argv[0]} <long|short|boundary> <mu> <nu> [path_to_debug_integrand.dat]")
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
        input_file = sys.argv[4]
    else:
        try:
            input_file = find_latest_output_file("debug_integrand.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            sys.exit(1)
        print(f"No file given; using latest run's debug integrand: {input_file}")

    plot_debug_integrand(range_type, mu, nu, input_file)
