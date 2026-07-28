import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from run_output_utils import find_latest_output_file

def plot_laser_fields(filepath):
    """
    Reads field data from a space-delimited text file and plots 
    the 6 electromagnetic components in a 3x2 matrix array.
    """
    # 1. Read data cleanly using pandas
    try:
        data = pd.read_csv(filepath, sep=" ")
    except Exception as e:
        print(f"Error reading file '{filepath}': {e}")
        sys.exit(1)

    # Verify that required columns exist
    required_cols = {'phi', 'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz'}
    if not required_cols.issubset(data.columns):
        print(f"Error: File must contain headers: {list(required_cols)}")
        sys.exit(1)

    # 2. Initialize a 3x2 subplot grid layout
    # sharex=True couples the X-axis so zooming into one panel zooms all of them seamlessly
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(14, 10), sharex=True)
    
    # Flatten the 2D axes array to a 1D list for easy sequential looping
    # Index mapping: 0=Ex, 1=Bx, 2=Ey, 3=By, 4=Ez, 5=Bz
    axes_flat = axes.flatten()

    # Define the mapping profiles for the 6 components
    components = [
        {'col': 'Ex', 'title': '$E_x$ (Electric Field X)', 'color': 'crimson'},
        {'col': 'Bx', 'title': '$B_x$ (Magnetic Field X)', 'color': 'navy'},
        {'col': 'Ey', 'title': '$E_y$ (Electric Field Y)', 'color': 'darkgreen'},
        {'col': 'By', 'title': '$B_y$ (Magnetic Field Y)', 'color': 'darkorange'},
        {'col': 'Ez', 'title': '$E_z$ (Electric Field Z)', 'color': 'purple'},
        {'col': 'Bz', 'title': '$B_z$ (Magnetic Field Z)', 'color': 'teal'}
    ]

    # 3. Loop and populate each subplot panel dynamically
    for i, comp in enumerate(components):
        ax = axes_flat[i]
        ax.plot(data['phi'] / (2 * np.pi), data[comp['col']], color=comp['color'], linewidth=1.5)
        ax.set_title(comp['title'], fontsize=11, fontweight='bold')
        ax.set_ylabel('Amplitude')
        ax.grid(True, linestyle='--', alpha=0.6)

    # 4. Label the bottom-most row X-axes explicitly
    axes[2, 0].set_xlabel('Phase $\\phi$ (cycles)', fontsize=12)
    axes[2, 1].set_xlabel('Phase $\\phi$ (cycles)', fontsize=12)

    # Adjust layout padding so titles and axes labels don't overlap
    plt.suptitle(f"Spacetime Field Profiles for {os.path.basename(filepath)}", fontsize=14, fontweight='bold', y=0.98)
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
        # No explicit path given: locate laser_field.dat under the most
        # recent run directory inside the configured output_folder.
        try:
            input_file = find_latest_output_file("laser_field.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            print(f"Run command like: python3 {sys.argv[0]} <path_to_laser_field.dat>")
            sys.exit(1)
        print(f"No file given; using latest run's laser field: {input_file}")
    else:
        input_file = sys.argv[1]

    plot_laser_fields(input_file)
    