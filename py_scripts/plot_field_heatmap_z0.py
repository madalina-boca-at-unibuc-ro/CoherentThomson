import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file

def read_header_comments(filepath):
    """
    Peeks at the file's leading '# ...' comment lines (written by
    Core::Laser::export_field_heatmap_z0): the 'axes_unit <name>' line and the
    'canonical-frame z=0 snapshot at ...' description line.
    """
    axes_unit = None
    snapshot_description = None
    with open(filepath) as f:
        for line in f:
            stripped = line.strip()
            if not stripped.startswith('#'):
                break
            comment = stripped.lstrip('#').strip()
            parts = comment.split()
            if len(parts) == 2 and parts[0] == 'axes_unit':
                axes_unit = parts[1]
            elif comment.startswith('canonical-frame'):
                snapshot_description = comment
    return axes_unit, snapshot_description

def plot_field_heatmap(filepath):
    """
    Reads the x/y/intensity grid from a space-delimited text file and plots
    it as a heat map, marking the snapshot time (the start of the pulse's
    flat-top plateau, per laser_delay) in the title.
    """
    axes_unit, snapshot_description = read_header_comments(filepath)

    try:
        data = pd.read_csv(filepath, sep=" ", comment='#')
    except Exception as e:
        print(f"Error reading file '{filepath}': {e}")
        sys.exit(1)

    required_cols = {'x', 'y', 'intensity'}
    if not required_cols.issubset(data.columns):
        print(f"Error: File must contain headers: {list(required_cols)}")
        sys.exit(1)

    # Reshape the flat (x, y, intensity) rows into a 2D grid for pcolormesh.
    grid = data.pivot(index='y', columns='x', values='intensity')

    axes_label = f" [{axes_unit}]" if axes_unit else ""

    plt.figure(figsize=(8, 7))
    mesh = plt.pcolormesh(grid.columns, grid.index, grid.values, shading='nearest', cmap='inferno')
    plt.colorbar(mesh, label=r"Intensity $E_x^2+E_y^2+E_z^2$")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel(f"$x${axes_label}", fontsize=12)
    plt.ylabel(f"$y${axes_label}", fontsize=12)

    plt.title("Field Intensity in Canonical z=0 Plane", fontsize=12, fontweight='bold')
    if snapshot_description:
        plt.suptitle(f"{snapshot_description} -- start of flat-top plateau", fontsize=9, y=0.98)

    # Save output plot into a 'png_folder' subfolder of the run directory,
    # alongside the .dat files.
    png_dir = os.path.join(os.path.dirname(filepath), "png_folder")
    os.makedirs(png_dir, exist_ok=True)
    output_name = os.path.basename(filepath).rsplit('.', 1)[0] + "_plot.png"
    output_img = os.path.join(png_dir, output_name)
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved field heatmap to {output_img}")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        # No explicit path given: locate laser_field_heatmap_z0.dat under the
        # most recent run directory inside the configured output_folder.
        try:
            input_file = find_latest_output_file("laser_field_heatmap_z0.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            print(f"Run command like: python3 {sys.argv[0]} <path_to_data_file>")
            sys.exit(1)
        print(f"No file given; using latest run's field heatmap: {input_file}")
    else:
        input_file = sys.argv[1]

    plot_field_heatmap(input_file)
