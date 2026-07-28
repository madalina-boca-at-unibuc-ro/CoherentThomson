import sys
import os
import typing
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from run_output_utils import find_latest_output_file

def read_length_units(filepath):
    """
    Peeks at the file's leading '# length units: <name>' comment line (written
    by Core::Particle::plot_electron_beam_scatter) and returns the unit name,
    or None if the file has no such line.
    """
    with open(filepath) as f:
        first_line = f.readline().strip()
    if first_line.startswith('#'):
        parts = first_line.lstrip('#').split(':', 1)
        if len(parts) == 2 and parts[0].strip() == 'length units':
            return parts[1].strip()
    return None

def plot_scatter(filepath):
    """
    Reads each electron's initial x/y/z position from a space-delimited text
    file and makes a 3D scatter plot, to check the generated beam's position
    and spread.
    """
    length_units = read_length_units(filepath)

    try:
        data = pd.read_csv(filepath, sep=" ", comment='#')
    except Exception as e:
        print(f"Error reading file '{filepath}': {e}")
        sys.exit(1)

    required_cols = {'x', 'y', 'z'}
    if not required_cols.issubset(data.columns):
        print(f"Error: File must contain headers: {list(required_cols)}")
        sys.exit(1)

    unit_suffix = f" [{length_units}]" if length_units else ""

    fig = plt.figure(figsize=(8, 8))
    ax = typing.cast(Axes3D, fig.add_subplot(projection='3d'))
    ax.scatter(data['x'].to_numpy(), data['y'].to_numpy(), data['z'].to_numpy(), c='royalblue', s=10, alpha=0.7,
               edgecolors='none')
    ax.scatter([0], [0], [0], c='red', marker='x', s=60, label='origin')
    ax.set_title("Electron Beam Initial Positions (lab frame)", fontsize=12, fontweight='bold')
    ax.set_xlabel(f"x{unit_suffix}", fontsize=11)
    ax.set_ylabel(f"y{unit_suffix}", fontsize=11)
    ax.set_zlabel(f"z{unit_suffix}", fontsize=11)
    ax.legend()

    # Save output plot into a 'png_folder' subfolder of the run directory,
    # alongside the .dat files.
    png_dir = os.path.join(os.path.dirname(filepath), "png_folder")
    os.makedirs(png_dir, exist_ok=True)
    output_name = os.path.basename(filepath).rsplit('.', 1)[0] + "_plot.png"
    output_img = os.path.join(png_dir, output_name)
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved scatter plot to {output_img}")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        # No explicit path given: locate electron_beam_scatter.dat under the
        # most recent run directory inside the configured output_folder.
        try:
            input_file = find_latest_output_file("electron_beam_scatter.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            print(f"Run command like: python3 {sys.argv[0]} <path_to_data_file>")
            sys.exit(1)
        print(f"No file given; using latest run's electron beam scatter positions: {input_file}")
    else:
        input_file = sys.argv[1]

    plot_scatter(input_file)
