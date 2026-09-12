import sys
import os
import typing
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

def read_length_units(filepath):
    """
    Peeks at the file's leading '# length units: <name>' comment line (written
    by Core::Detector::plot_detector_scatter) and returns the unit name, or
    None if the file has no such line.
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
    Reads real x/y/z detector point positions from a space-delimited text
    file and makes a 3D scatter plot, to check the detector's position and
    orientation.
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
    ax.set_title("Detector Point Positions (lab frame)", fontsize=12, fontweight='bold')
    ax.set_xlabel(f"x{unit_suffix}", fontsize=11)
    ax.set_ylabel(f"y{unit_suffix}", fontsize=11)
    ax.set_zlabel(f"z{unit_suffix}", fontsize=11)
    ax.legend()

    # Save output plot into a per-module subfolder of the run directory's 'png_folder'
    # (mirroring py_scripts/'s own laser/detector/particle/radiation/debug layout).
    png_dir = os.path.join(os.path.dirname(filepath), "png_folder", MODULE_NAME)
    os.makedirs(png_dir, exist_ok=True)
    output_name = os.path.basename(filepath).rsplit('.', 1)[0] + "_plot.png"
    output_img = os.path.join(png_dir, output_name)
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved scatter plot to {output_img}")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        # No explicit run folder given: locate detector_scatter.dat under the
        # most recent run directory inside the configured output_folder.
        try:
            input_file = find_latest_output_file("detector_scatter.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            print(f"Run command like: python3 {sys.argv[0]} <path_to_run_folder>")
            sys.exit(1)
        print(f"No folder given; using latest run: {input_file}")
    else:
        input_file = os.path.join(sys.argv[1], "detector_scatter.dat")
        if not os.path.exists(input_file):
            print(f"Error: '{input_file}' not found")
            sys.exit(1)

    plot_scatter(input_file)
