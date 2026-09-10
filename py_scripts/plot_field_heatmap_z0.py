import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file
from plot_radiation_field import read_config_value

def get_laser_lg_w0_in_axes_units(filepath, axes_unit):
    """
    Returns laser_lg_w0 (the Laguerre-Gauss beam waist) expressed in the same units as the
    heatmap's own x/y columns (axes_unit, read from the .dat file's header -- currently always
    'lambda', since main.cpp hardcodes heatmap_axes_unit to that), or None if this run's
    laser_type isn't 'laguerre_gauss' (no w0 to label with) or laser_lg_w0's own config unit
    doesn't match axes_unit (defensive: laser_lg_w0 is always given in 'lambda' units by this
    project's own config convention, matching axes_unit -- but degrades to None instead of
    silently mislabeling the axis if that ever changes, rather than doing a full unit conversion
    for a case that shouldn't occur in practice).

    Reads the run's own config.cfg snapshot (Core::IoUtils::copy_config_to_run_directory), not
    the live repo config, per the convention every other script here already follows.
    """
    config_path = os.path.join(os.path.dirname(filepath), 'config.cfg')
    if not os.path.exists(config_path):
        return None
    try:
        laser_type, _ = read_config_value('laser_type', config_path)
        if laser_type != 'laguerre_gauss':
            return None
        w0_str, w0_unit = read_config_value('laser_lg_w0', config_path)
    except ValueError:
        return None
    if w0_unit != axes_unit:
        print(f"Warning: laser_lg_w0's unit ('{w0_unit}') does not match the heatmap's axes_unit "
              f"('{axes_unit}') -- skipping the w0-unit axis labels.")
        return None
    return float(w0_str)

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
    w0 = get_laser_lg_w0_in_axes_units(filepath, axes_unit)

    fig, ax = plt.subplots(figsize=(8, 7))
    mesh = ax.pcolormesh(grid.columns, grid.index, grid.values, shading='nearest', cmap='inferno')
    fig.colorbar(mesh, ax=ax, label=r"Intensity $E_x^2+E_y^2+E_z^2$")
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(f"$x${axes_label}", fontsize=12)
    ax.set_ylabel(f"$y${axes_label}", fontsize=12)

    # Secondary top/right axes in units of the beam waist w0 (laguerre_gauss runs only -- w0 is
    # None otherwise, see get_laser_lg_w0_in_axes_units), alongside the primary bottom/left axes
    # already in axes_unit -- both shown at once rather than replacing one with the other, so
    # either physical position or w0-multiples can be read directly off the plot.
    if w0 is not None and w0 > 0:
        secax_x = ax.secondary_xaxis('top', functions=(lambda x: x / w0, lambda x: x * w0))
        secax_x.set_xlabel("$x / w_0$", fontsize=12)
        secax_y = ax.secondary_yaxis('right', functions=(lambda y: y / w0, lambda y: y * w0))
        secax_y.set_ylabel("$y / w_0$", fontsize=12)

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
