import sys
import os
import math
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import erf

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'radiation'))
from utils.run_output_utils import get_output_folder, find_latest_run_dir
from plot_field import read_config_value, read_length

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

def read_log_length(run_dir, label):
    """A length from run_log.txt's "label  <a.u.> a.u.  <lambda> lambda" line, in lambda; None if absent."""
    with open(os.path.join(run_dir, "run_log.txt")) as log_file:
        match = re.search(r"^  " + re.escape(label) + r"\s+\S+ a\.u\.\s+(\S+) lambda", log_file.read(), re.M)
    return float(match.group(1)) if match else None

def read_beam_geometry(run_dir, config_path):
    """
    Reads the longitudinal beam parameters (in lambda): the cylinder height H (a lattice stack's width), the Gaussian
    edge smoothing sigma_z, and the center z_c. H and sigma_z come from run_log.txt, which has them resolved (they can
    be given in lambda_1, which only the solver knows); runs without run_log.txt fall back to config.cfg. Runs made
    before beam_sigma_z existed are read as sigma_z = 0 (the sharp-edged cylinder they used).
    """
    z_c = read_length("beam_center_z", config_path)
    has_log = os.path.exists(os.path.join(run_dir, "run_log.txt"))
    H = read_log_length(run_dir, "cylinder height") if has_log else read_length("beam_cylinder_height", config_path)
    sigma_z = read_log_length(run_dir, "sigma_z") if has_log else None
    if sigma_z is None:
        try:
            sigma_z = read_length("beam_sigma_z", config_path)
        except ValueError:
            print("No 'beam_sigma_z' in this run's config (older run); using sigma_z = 0")
            sigma_z = 0.0
    # Runs from before the 2026-11-09 rework: a stack of M sheets, spacing d, drawn over H = M d.
    with open(config_path) as config_file:
        sheets = [line.split()[1] for line in config_file if line.split()[:1] == ["beam_stack_sheets"]]
    if sheets and int(sheets[0]) > 1 and read_config_value_or_none("beam_stack_spacing", config_path) is not None:
        H = int(sheets[0]) * read_length("beam_stack_spacing", config_path)
    return H, sigma_z, z_c

def read_config_value_or_none(key, config_path):
    try:
        return read_config_value(key, config_path)[0]
    except ValueError:
        return None

def read_modulation(run_dir, config_path):
    """
    The modulation's map, from the run's log: {zeta, phase, Lambda_mod, H, sigma_z, lattice}, in lambda and rad, or
    None for an unmodulated beam.
    """
    switch = read_config_value_or_none("beam_modulation", config_path)
    legacy = read_config_value_or_none("beam_stack_modulation", config_path)  # before the 2026-11-09 rework
    if switch is not None:
        if switch != "true":
            return None
    elif legacy is not None and legacy != "true":
        return None
    # A few runs from before beam_stack_modulation existed have the beam_modulation_* keys alone, and are bunched.
    zeta = read_config_value_or_none("beam_modulation_zeta", config_path)
    if zeta is None:
        return None
    with open(os.path.join(run_dir, "run_log.txt")) as log_file:
        log = log_file.read()
    phase = math.pi * float(re.search(r"phase\s+(\S+) pi", log).group(1))
    H, sigma_z, _ = read_beam_geometry(run_dir, config_path)
    if switch is None:
        H = read_log_length(run_dir, "stack length H = M d")
    lattice = read_config_value_or_none("beam_distribution", config_path) == "polar_lattice"
    return {"zeta": float(zeta), "phase": phase, "Lambda_mod": read_log_length(run_dir, "Lambda_mod"), "H": H,
            "sigma_z": sigma_z, "lattice": lattice}

def modulated_density(bins, modulation, z_c):
    """
    The modulated density averaged over each histogram bin, from a fine sample of the drift map on the uniform
    density over H (plus the Gaussian wings of a random beam).
    """
    K = 2 * math.pi / modulation["Lambda_mod"]
    H = modulation["H"]
    uniform = (np.arange(4_000_000) + 0.5) / 4_000_000 * H - H / 2
    if modulation["sigma_z"] > 0.0:
        uniform = uniform + np.random.default_rng(0).normal(0.0, modulation["sigma_z"], uniform.size)
    moved = z_c + uniform + modulation["zeta"] / K * np.sin(K * uniform + modulation["phase"])
    return np.histogram(moved, bins=bins, density=True)[0]

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
    Returns each electron's position along the beam axis, in lambda, from electron_beam_scatter.dat, and its charge
    weight (the file's w column, written only when the weights differ: a Gaussian lattice stack), or None.
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
    weights = data['w'].to_numpy() if 'w' in data.columns else None
    return data[['x', 'y', 'z']].to_numpy() @ axis, weights

def lattice_bins(z):
    """
    Bins for a lattice stack's sheets, about four sheets each, since finer bins hold 0, 1 or 2 discrete sheets and show
    gaps that aren't in the density. Sheets are counted from the distinct positions, so this works at H = 0 too.
    Returns (bins, sheet count).
    """
    scale = max(np.ptp(z), 1e-300)
    sheets = len(np.unique(np.round((z - z.min()) / scale, 10)))
    half_sheet = scale / (2 * max(sheets - 1, 1))
    return np.linspace(z.min() - half_sheet, z.max() + half_sheet, max(5, sheets // 4) + 1), sheets

def read_tail_cut(run_dir):
    """A Gaussian lattice stack's beam_stack_sigma_z_cut, from run_log.txt; None for any other beam."""
    with open(os.path.join(run_dir, "run_log.txt")) as log_file:
        match = re.search(r"^  tail cut\s+(\S+) sigma_z", log_file.read(), re.M)
    return float(match.group(1)) if match else None

def plot_modulated_density(run_dir, config_path, modulation):
    """
    A modulated beam: the charge per bin against the drift map's density, with the histogram from
    electron_beam_scatter.dat when the run has it. For a lattice stack the bins cover only the sheets and average about
    four of them (lattice_bins), weighted by the sheets' charges; for a random beam there are ~sqrt(N) bins.
    """
    z_c = read_length("beam_center_z", config_path)
    H = modulation["H"]
    cut = read_tail_cut(run_dir) or 5.0
    spread = H / 2 + cut * modulation["sigma_z"] + modulation["zeta"] * modulation["Lambda_mod"] / (2 * math.pi)
    fig, ax = plt.subplots(figsize=(8, 5))
    scatter_path = os.path.join(run_dir, "electron_beam_scatter.dat")
    weights = None
    if os.path.exists(scatter_path):
        z, weights = read_axial_positions(scatter_path, config_path)
        if modulation["lattice"]:
            bins, sheets = lattice_bins(z)
            label = f"sheets ({sheets}, {len(z)} electrons{', weighted' if weights is not None else ''})"
        else:
            n_bins = int(np.clip(np.sqrt(len(z)), 20, 120))
            bins = np.linspace(z.min(), z.max(), n_bins + 1)
            label = f"sampled electrons (N = {len(z)})"
    else:
        print(f"No '{scatter_path}' (needs plot_beam_scatter=true); plotting the analytic density only")
        z = None
        bins = np.linspace(z_c - 1.05 * spread, z_c + 1.05 * spread, 201)
    if z is not None:
        ax.hist(z, bins=bins, weights=weights, density=True, histtype='stepfilled', color='lightsteelblue',
                label=label)
    ax.stairs(modulated_density(bins, modulation, z_c), bins, color='crimson', lw=1.5, label="drift map")
    for edge in (z_c - H / 2, z_c + H / 2):
        ax.axvline(edge, color='gray', ls='--', lw=1)
    ax.set_title(f"{'Bunched stack' if modulation['lattice'] else 'Modulated beam'} along the beam axis   ($\\zeta$ = {modulation['zeta']:g}, "
                 f"$H/\\Lambda_{{mod}}$ = {H / modulation['Lambda_mod']:g})", fontsize=12, fontweight='bold')
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

def plot_density(run_dir):
    config_path = os.path.join(run_dir, "config.cfg")
    modulation = read_modulation(run_dir, config_path)
    if modulation is not None:
        plot_modulated_density(run_dir, config_path, modulation)
        return

    H, sigma_z, z_c = read_beam_geometry(run_dir, config_path)
    if H == 0.0 and sigma_z == 0.0:
        print("beam_cylinder_height = 0 and beam_sigma_z = 0: every electron sits at z = beam_center_z, "
              "no density to plot")
        sys.exit(0)

    lattice = read_config_value_or_none("beam_distribution", config_path) == "polar_lattice"
    cut = read_tail_cut(run_dir)
    half_extent = H / 2 + (cut or 5.0) * sigma_z
    z_grid = np.linspace(z_c - 1.2 * half_extent, z_c + 1.2 * half_extent, 2001)

    fig, ax = plt.subplots(figsize=(8, 5))

    scatter_path = os.path.join(run_dir, "electron_beam_scatter.dat")
    if os.path.exists(scatter_path):
        z, weights = read_axial_positions(scatter_path, config_path)
        if lattice:
            # A lattice stack's sheets (charges from the Gaussian profile when sigma_z > 0): about four per bin.
            bins, sheets = lattice_bins(z)
            label = f"sheets ({sheets}, {len(z)} electrons{', weighted' if weights is not None else ''})"
        else:
            n_bins = int(np.clip(np.sqrt(len(z)), 20, 120))  # ~sqrt(N) bins, so small test beams aren't all noise
            bins = np.linspace(z_grid[0], z_grid[-1], n_bins + 1)
            label = f"sampled electrons (N = {len(z)})"
        ax.hist(z, bins=bins, weights=weights, density=True, histtype='stepfilled', color='lightsteelblue',
                label=label)
        w = np.ones_like(z) if weights is None else weights
        mean = np.average(z, weights=w)
        var = np.average((z - mean) ** 2, weights=w)
        # A cut Gaussian (a lattice stack's tails stop at cut sigma_z) has a slightly smaller variance than H^2/12 +
        # sigma_z^2; at the default cut of 6 the difference is ~1e-7 relative.
        print(f"N = {len(z)}" + (f" in {sheets} sheets" if lattice else ""))
        print(f"mean z:   {'weighted' if weights is not None else 'sampled'} {mean:+.5e}   expected {z_c:+.5e}  lambda")
        print(f"var z:    {'weighted' if weights is not None else 'sampled'} {var:.5e}    expected "
              f"{H**2 / 12 + sigma_z**2:.5e}  lambda^2")
        print(f"fraction outside [z_c - H/2, z_c + H/2]: {np.sum(w[np.abs(z - z_c) > H / 2]) / np.sum(w):.4f}")
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
