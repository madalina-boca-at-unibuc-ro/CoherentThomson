import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file, get_laser_period

# Core::PhysUtils::AtomicUnits (phys_utils.hpp): c in the solver's own atomic units (1/alpha, 2018
# CODATA), and m_0 = 1.0 there -- so mc (the momentum unit) is just C_LIGHT itself, no separate mass
# factor needed. Mirrors phys_utils.hpp's own value independently (no shared constants module between
# C++ and Python) -- keep in sync if that one changes.
C_LIGHT = 137.035999084

def _plot_four_vector_figure(data, electron_ids, electron_colors, cols, panel_titles, suptitle, output_img,
                             scale, ylabel):
    """
    Renders one 2x2 figure, one panel per four-vector component (cols[i] -> panel_titles[i]),
    each with its own independently auto-scaled y-axis -- unlike an earlier version of this script,
    which overlaid x1/x2/x3 (or p1/p2/p3) on one shared panel: components with very different
    absolute magnitudes (e.g. a beam's large-but-nearly-constant transverse starting offset in
    x1/x2 vs. x3's genuine large excursion from 0) would then visually flatten the smaller-looking
    one even though it's actually varying substantially, purely because the shared axis had to
    stretch to fit the other components too. Splitting into one panel per component sidesteps that
    entirely. Color still encodes electron identity, shared across both this figure and its
    position/momentum counterpart via the same electron_colors mapping.

    `scale` converts each plotted component from atomic units into the desired display unit
    (divides the raw column value) -- lambda for position (all four components are lengths, x0=ct
    included), mc for momentum (all four are momentum-like, p0=E/c included) -- and `ylabel` names
    that unit; see plot_electron_trajectory's own doc comment for how each scale is computed.
    """
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8), sharex=True)
    axes_flat = axes.flatten()

    for ax, col, title in zip(axes_flat, cols, panel_titles):
        for eid in electron_ids:
            electron_data = data[data['electron_id'] == eid]
            ax.plot(electron_data['tau_over_T'], electron_data[col] / scale, color=electron_colors[eid],
                    linewidth=1.5)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_ylabel(ylabel)
        ax.grid(True, linestyle='--', alpha=0.6)

    # Electron-identity legend (color), shown once since it applies to every panel.
    electron_handles = [plt.Line2D([0], [0], color=electron_colors[eid], linewidth=2, label=f'electron {eid}')
                        for eid in electron_ids]
    fig.legend(handles=electron_handles, loc='upper center', ncol=min(len(electron_ids), 10),
              bbox_to_anchor=(0.5, 0.90), fontsize=9, frameon=False)

    xlabel = r'Proper time $\tau / T$  ($T = 2\pi/\omega$, laser period)'
    axes[1, 0].set_xlabel(xlabel, fontsize=12)
    axes[1, 1].set_xlabel(xlabel, fontsize=12)

    plt.suptitle(suptitle, fontsize=14, fontweight='bold', y=0.99)
    plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.85))

    plt.savefig(output_img, dpi=200)
    print(f"Successfully rendered and saved layout array to {output_img}")

def plot_electron_trajectory(filepath):
    """
    Reads trajectory data (one or more electrons, tagged by 'electron_id') from a
    space-delimited text file and plots the position and momentum four-vectors as a function of
    proper time (in units of the laser period T = 2*pi/omega, read from the run's own config.cfg)
    as two separate figures -- one for position (x0..x3, in units of lambda), one for momentum
    (p0..p3, in units of mc) -- each a 2x2 grid with one panel per component. Color encodes which
    electron a line belongs to.

    lambda (the laser wavelength, in atomic units) is derived from the same T this script already
    reads for the x-axis: lambda = c*T (T = 2*pi/omega is a time, lambda = 2*pi*c/omega a length --
    they differ by exactly the c factor, mirroring Core::IoUtils::convert_unit_to_number's own
    'lambda' unit branch, which is 2*pi/omega*c). mc is just C_LIGHT (m_0 = 1.0 in these atomic
    units, see Core::PhysUtils::AtomicUnits).
    """
    # 1. Read data cleanly using pandas (electron.dat starts with a '#'-prefixed
    # metadata line before the header row; blank lines separate each electron's
    # block and are skipped automatically)
    try:
        data = pd.read_csv(filepath, sep=" ", comment='#')
    except Exception as e:
        print(f"Error reading file '{filepath}': {e}")
        sys.exit(1)

    # Verify that required columns exist
    required_cols = {'electron_id', 'tau', 'x0', 'x1', 'x2', 'x3', 'p0', 'p1', 'p2', 'p3'}
    if not required_cols.issubset(data.columns):
        print(f"Error: File must contain headers: {list(required_cols)}")
        sys.exit(1)

    T = get_laser_period(os.path.dirname(filepath))
    lambda_au = T * C_LIGHT
    data = data.copy()
    data['tau_over_T'] = data['tau'] / T

    electron_ids = sorted(data['electron_id'].unique())
    # tab10 is a fixed 10-color qualitative palette; electron.dat never carries more
    # than 10 electrons (see main.cpp), so identity never wraps/cycles.
    cmap = plt.get_cmap('tab10')
    electron_colors = {eid: cmap(i) for i, eid in enumerate(electron_ids)}

    png_dir = os.path.join(os.path.dirname(filepath), "png_folder")
    os.makedirs(png_dir, exist_ok=True)
    basename = os.path.basename(filepath).rsplit('.', 1)[0]

    _plot_four_vector_figure(
        data, electron_ids, electron_colors,
        cols=['x0', 'x1', 'x2', 'x3'],
        panel_titles=['$x^0\\, (c t)$', '$x^1$', '$x^2$', '$x^3$'],
        suptitle=f"Electron Trajectory (Position) for {os.path.basename(filepath)}",
        output_img=os.path.join(png_dir, f"{basename}_position_profile.png"),
        scale=lambda_au, ylabel=r'Position [$\lambda$]',
    )

    _plot_four_vector_figure(
        data, electron_ids, electron_colors,
        cols=['p0', 'p1', 'p2', 'p3'],
        panel_titles=['$p^0\\,(E/c)$', '$p^1$', '$p^2$', '$p^3$'],
        suptitle=f"Electron Trajectory (Momentum) for {os.path.basename(filepath)}",
        output_img=os.path.join(png_dir, f"{basename}_momentum_profile.png"),
        scale=C_LIGHT, ylabel=r'Momentum [$mc$]',
    )

    # Display both interactive windows.
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        # No explicit path given: locate electron.dat under the most recent
        # run directory inside the configured output_folder.
        try:
            input_file = find_latest_output_file("electron.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            print(f"Run command like: python3 {sys.argv[0]} <path_to_electron.dat>")
            sys.exit(1)
        print(f"No file given; using latest run's trajectory: {input_file}")
    else:
        input_file = sys.argv[1]

    plot_electron_trajectory(input_file)
