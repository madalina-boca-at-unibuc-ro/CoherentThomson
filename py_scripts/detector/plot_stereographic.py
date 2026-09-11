import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file

def read_axes_unit(filepath):
    """
    Peeks at the file's leading '# axes_unit <name>' comment line (written by
    Core::Detector::plot_detector) and returns the unit name, or None if the
    file has no such line.
    """
    with open(filepath) as f:
        first_line = f.readline().strip()
    if first_line.startswith('#'):
        parts = first_line.lstrip('#').split()
        if len(parts) == 2 and parts[0] == 'axes_unit':
            return parts[1]
    return None

def plot_stereographic(filepath):
    """
    Reads stereographic coordinates from a space-delimited text file
    and makes a scatter plot.
    """
    axes_unit = read_axes_unit(filepath)

    try:
        data = pd.read_csv(filepath, sep=" ", comment='#')
    except Exception as e:
        print(f"Error reading file '{filepath}': {e}")
        sys.exit(1)

    required_cols = {'x_proj', 'y_proj'}
    if not required_cols.issubset(data.columns):
        print(f"Error: File must contain headers: {list(required_cols)}")
        sys.exit(1)

    x_label = "$x_{proj}$" + (f" [{axes_unit}]" if axes_unit else "")
    y_label = "$y_{proj}$" + (f" [{axes_unit}]" if axes_unit else "")

    plt.figure(figsize=(8, 8))
    plt.scatter(data['x_proj'], data['y_proj'], c='royalblue', s=10, alpha=0.7, edgecolors='none')
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Stereographic Projection of Spherical Detector", fontsize=12, fontweight='bold')
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.6)

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
        # No explicit path given: locate detector_stereographic.dat under the
        # most recent run directory inside the configured output_folder.
        try:
            input_file = find_latest_output_file("detector_stereographic.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            print(f"Run command like: python3 {sys.argv[0]} <path_to_data_file>")
            sys.exit(1)
        print(f"No file given; using latest run's stereographic projection: {input_file}")
    else:
        input_file = sys.argv[1]

    plot_stereographic(input_file)
