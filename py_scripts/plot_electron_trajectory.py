import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file

def plot_electron_trajectory(filepath):
    """
    Reads trajectory data from a space-delimited text file and plots
    the trajectory components (position, momentum) as a function of proper time in a 2x2 matrix array. To the left the temporal component is plotted, to the right the three spatial parts
    """
    # 1. Read data cleanly using pandas (electron.dat starts with '#'-prefixed
    # metadata lines, e.g. '# store_trajectory value ...', before the header row)
    try:
        data = pd.read_csv(filepath, sep=" ", comment='#')
    except Exception as e:
        print(f"Error reading file '{filepath}': {e}")
        sys.exit(1)

    # Verify that required columns exist
    required_cols = {'tau', 'x0', 'x1', 'x2', 'x3', 'p0', 'p1', 'p2', 'p3'}
    if not required_cols.issubset(data.columns):
        print(f"Error: File must contain headers: {list(required_cols)}")
        sys.exit(1)

    # 2. Initialize a 2x2 subplot grid layout
    # sharex=True couples the X-axis so zooming into one panel zooms all of them seamlessly
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 7), sharex=True)

    # Flatten the 2D axes array to a 1D list for easy sequential looping
    axes_flat = axes.flatten()

    # Define the mapping profiles for the 6 panels. Timelike components (left
    # column) are a single line; spacelike components (right column) are the
    # three x/y/z lines overlaid on the same panel.
    components = [
        {'cols': ['x0'], 'title': '$x^0\\, (c t)$', 'colors': ['crimson']},
        {'cols': ['x1', 'x2', 'x3'], 'title': '${\\bf x}$ $(x^1,\\, x^2,\\, x^3)$',
         'colors': ['navy', 'steelblue', 'skyblue']},
        {'cols': ['p0'], 'title': '$p^0\\,(E/c)$', 'colors': ['darkgreen']},
        {'cols': ['p1', 'p2', 'p3'], 'title': '${\\bf p}$ $(p^1,\\, p^2,\\,p^3)$',
         'colors': ['darkorange', 'goldenrod', 'khaki']},
    ]

    # 3. Loop and populate each subplot panel dynamically
    for i, comp in enumerate(components):
        ax = axes_flat[i]
        for col, color in zip(comp['cols'], comp['colors']):
            ax.plot(data['tau'], data[col], color=color, linewidth=1.5, label=col)
        ax.set_title(comp['title'], fontsize=11, fontweight='bold')
        ax.set_ylabel('Amplitude')
        ax.grid(True, linestyle='--', alpha=0.6)
        if len(comp['cols']) > 1:
            ax.legend(fontsize=9)

    # 4. Label the bottom-most row X-axes explicitly
    axes[1, 0].set_xlabel(r'Proper time $\tau$', fontsize=12)
    axes[1, 1].set_xlabel(r'Proper time $\tau$', fontsize=12)

    # Adjust layout padding so titles and axes labels don't overlap
    plt.suptitle(f"Electron Trajectory for {os.path.basename(filepath)}", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()

    # Save a high-res diagnostic image into a 'png_folder' subfolder of the
    # run directory, alongside the .dat files.
    png_dir = os.path.join(os.path.dirname(filepath), "png_folder")
    os.makedirs(png_dir, exist_ok=True)
    output_name = os.path.basename(filepath).rsplit('.', 1)[0] + "_profile.png"
    output_img = os.path.join(png_dir, output_name)
    plt.savefig(output_img, dpi=200)
    print(f"Successfully rendered and saved layout array to {output_img}")

    # Display the interactive window
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
