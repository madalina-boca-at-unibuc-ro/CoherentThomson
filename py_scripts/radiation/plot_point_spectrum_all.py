"""
Driver over plot_point_spectrum.py: for each range type ('long', 'short', 'boundary', 'total') draws one
overview figure with |F^{mu nu}| of the 6 independent components (mu < nu; F^{nu mu} = -F^{mu nu}, diagonal
zero) vs omega/omega_1, then the per-component Re/Im/Abs/Phase figure for all 24 (range, mu, nu) combinations.
Everything goes to png_folder/radiation/point_spectrum/.

Usage: python3 plot_point_spectrum_all.py [path_to_run_folder]
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

from plot_point_spectrum import (
    OMEGA_LABEL,
    RANGE_TYPES,
    get_point_component,
    get_point_spectrum_png_dir,
    mark_harmonics,
    plot_point_spectrum,
    read_radiation_data,
    resolve_radiation_file,
)

COMPONENT_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


def plot_point_spectrum_overview(range_type, radiation_filepath, data):
    """2x3 grid of |F^{mu nu}| vs omega/omega_1 for the 6 independent components of `range_type`."""
    png_dir = get_point_spectrum_png_dir(radiation_filepath)
    n_screen = data['i_screen'].nunique()

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True)
    for ax, (mu, nu) in zip(axes.flat, COMPONENT_PAIRS):
        field = get_point_component(data, range_type, mu, nu)
        for i_screen, subset in data.groupby('i_screen'):
            subset = subset.sort_values('omega')
            ax.plot(subset['omega'].to_numpy(), np.abs(field[subset.index.to_numpy()]), marker='.',
                    markersize=3, label=f'screen point {i_screen}')
        mark_harmonics(ax, data['omega'].to_numpy())
        ax.set_title(f"$|F^{{{mu}{nu}}}_{{\\mathrm{{{range_type}}}}}|$")
        ax.grid(True, which='major')
        if n_screen > 1:
            ax.legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel(OMEGA_LABEL)

    fig.suptitle(f"|F| spectrum, {range_type} ({n_screen} screen point{'s' if n_screen != 1 else ''})",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()

    output_img = os.path.join(png_dir, f"point_spectrum_overview_{range_type}.png")
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved plot to {output_img}")
    plt.close(fig)


if __name__ == "__main__":
    input_file = resolve_radiation_file(sys.argv[1] if len(sys.argv) >= 2 else None)
    data = read_radiation_data(input_file)

    for range_type in RANGE_TYPES:
        plot_point_spectrum_overview(range_type, input_file, data)
        for mu, nu in COMPONENT_PAIRS:
            print(f"Plotting {range_type} F^{{{mu}{nu}}}...")
            plot_point_spectrum(range_type, mu, nu, input_file, data=data)
