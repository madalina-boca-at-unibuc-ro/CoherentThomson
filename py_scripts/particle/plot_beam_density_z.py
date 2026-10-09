import sys
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import erf

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'radiation'))
from utils.run_output_utils import get_output_folder, find_latest_run_dir
from plot_field import read_config_value, read_length

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

def read_beam_geometry(config_path):
    """
    Reads the longitudinal beam parameters (in lambda) from a run's own config.cfg snapshot:
    the cylinder height H, the Gaussian edge smoothing sigma_z, and the center z_c. Runs made
    before beam_sigma_z existed are read as sigma_z = 0 (the sharp-edged cylinder they used).
    """
    H = read_length("beam_cylinder_height", config_path)
    z_c = read_length("beam_center_z", config_path)
    try:
        sigma_z = read_length("beam_sigma_z", config_path)
    except ValueError:
        print("No 'beam_sigma_z' in this run's config (older run); using sigma_z = 0")
        sigma_z = 0.0
    return H, sigma_z, z_c

def analytic_density(z, H, sigma_z, z_c):
    """
    Normalized longitudinal density of Core::Particle::generate_electron's z = z_c + U(-H/2, H/2) +
    N(0, sigma_z^2) (theory/gaussian_smoothed_electron_cylinder.md): the top-hat convolved with the
    Gaussian, reducing to the plain top-hat at sigma_z = 0 and to the Gaussian alone at H = 0.
    """
    s = z - z_c
    if sigma_z == 0.0:
        return np.where(np.abs(s) <= H / 2, 1.0 / H, 0.0)
    if H == 0.0:
        return np.exp(-s**2 / (2 * sigma_z**2)) / (math.sqrt(2 * math.pi) * sigma_z)
    w = math.sqrt(2) * sigma_z
    return (erf((s + H / 2) / w) - erf((s - H / 2) / w)) / (2 * H)

def read_axial_positions(scatter_path, config_path):
    """
    Returns each electron's position along the beam axis, in lambda, from electron_beam_scatter.dat.
    The stored positions are after the rotation that aligns the cylinder axis with the mean momentum
    (MathUtils::rotation_matrix_from_direction), so they are projected back onto that axis (Oz when the
    mean momentum is zero); the projection is exactly the sampled z (plus beam_center_z).
    """
    with open(scatter_path) as f:
        first_line = f.readline().strip()
    if first_line != "# length units: lambda":
        raise ValueError(f"Expected '# length units: lambda' in '{scatter_path}', got '{first_line}'")
    data = pd.read_csv(scatter_path, sep=" ", comment='#')

    # Only the direction of the mean momentum matters, so each component's own unit is irrelevant as
    # long as they all share one (the repo's configs always give them in 'mc').
    p_mean = np.array([float(read_config_value(k, config_path)[0]) for k in ("average_px", "average_py", "average_pz")])
    norm = np.linalg.norm(p_mean)
    axis = p_mean / norm if norm >= 1e-12 else np.array([0.0, 0.0, 1.0])
    return data[['x', 'y', 'z']].to_numpy() @ axis

def plot_density(run_dir):
    config_path = os.path.join(run_dir, "config.cfg")
    H, sigma_z, z_c = read_beam_geometry(config_path)
    if H == 0.0 and sigma_z == 0.0:
        print("beam_cylinder_height = 0 and beam_sigma_z = 0: every electron sits at z = beam_center_z, "
              "no density to plot")
        sys.exit(0)

    half_extent = H / 2 + 5 * sigma_z
    z_grid = np.linspace(z_c - 1.2 * half_extent, z_c + 1.2 * half_extent, 2001)

    fig, ax = plt.subplots(figsize=(8, 5))

    scatter_path = os.path.join(run_dir, "electron_beam_scatter.dat")
    if os.path.exists(scatter_path):
        z = read_axial_positions(scatter_path, config_path)
        n_bins = int(np.clip(np.sqrt(len(z)), 20, 120))  # ~sqrt(N) bins, so small test beams aren't all noise
        bins = np.linspace(z_grid[0], z_grid[-1], n_bins + 1)
        ax.hist(z, bins=bins, density=True, histtype='stepfilled', color='lightsteelblue',
                label=f"sampled electrons (N = {len(z)})")
        print(f"N = {len(z)}")
        print(f"mean z:   sampled {z.mean():+.5f}   expected {z_c:+.5f}  lambda")
        print(f"var z:    sampled {z.var():.5f}    expected {H**2 / 12 + sigma_z**2:.5f}  lambda^2")
        print(f"fraction outside [z_c - H/2, z_c + H/2]: {np.mean(np.abs(z - z_c) > H / 2):.4f}")
    else:
        print(f"No '{scatter_path}' (needs plot_beam_scatter=true); plotting the analytic density only")

    ax.plot(z_grid, analytic_density(z_grid, H, sigma_z, z_c), color='crimson', lw=1.8, label="analytic density")
    if H > 0.0:
        for edge in (z_c - H / 2, z_c + H / 2):
            ax.axvline(edge, color='gray', ls='--', lw=1)
    ax.set_title(f"Electron density along the beam axis   (H = {H:g} $\\lambda$, $\\sigma_z$ = {sigma_z:g} $\\lambda$)",
                 fontsize=12, fontweight='bold')
    ax.set_xlabel(r"z along beam axis [$\lambda$]", fontsize=11)
    ax.set_ylabel(r"normalized density [1/$\lambda$]", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.legend()

    png_dir = os.path.join(run_dir, "png_folder", MODULE_NAME)
    os.makedirs(png_dir, exist_ok=True)
    output_img = os.path.join(png_dir, "beam_density_z_plot.png")
    plt.savefig(output_img, dpi=200, bbox_inches='tight')
    print(f"Successfully saved density plot to {output_img}")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        try:
            run_dir = find_latest_run_dir(get_output_folder())
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            print(f"Run command like: python3 {sys.argv[0]} <path_to_run_folder>")
            sys.exit(1)
        print(f"No folder given; using latest run: {run_dir}")
    else:
        run_dir = sys.argv[1]
        if not os.path.exists(os.path.join(run_dir, "config.cfg")):
            print(f"Error: '{run_dir}' has no config.cfg")
            sys.exit(1)

    plot_density(run_dir)
