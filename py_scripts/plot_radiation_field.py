import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file, DEFAULT_CONFIG_PATH

def read_config_value(key, config_path=DEFAULT_CONFIG_PATH):
    """
    Reads a single 'key value [unit]' line out of the .cfg file and returns
    (value, unit) as strings (unit is '' if the key has no unit suffix).
    """
    with open(config_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            parts = stripped.split()
            if parts[0] == key and len(parts) >= 2:
                return parts[1], (parts[2] if len(parts) >= 3 else '')
    raise ValueError(f"No '{key}' key found in '{config_path}'")

def _angle_radians(value_str, unit_str):
    """Converts a 'value [pi]' config pair to radians (bare numbers are assumed radians)."""
    value = float(value_str)
    if unit_str == 'pi':
        return value * np.pi
    if unit_str == '':
        return value
    raise ValueError(f"Unsupported angle unit '{unit_str}' (expected 'pi' or none)")

def get_spherical_cell_edges(config_path=DEFAULT_CONFIG_PATH):
    """
    Returns (x_edges, y_edges), each shape (N_theta+1, N_phi+1): the stereographic-projection
    coordinates of the *cell corners* of the spherical detector's (N_theta, N_phi) grid, i.e. the
    same MathUtils formula Core::Detector::SphericalDetector::get_stereographic_projection uses
    (rho = R*sin(theta)/(1+cos(theta)), x=rho*cos(phi), y=rho*sin(phi)), evaluated half a grid
    step beyond each cell center instead of at the centers themselves.

    Needed because pcolormesh's shading='auto'/'nearest' infers cell edges from centers assuming
    they're roughly monotonic -- true for the rectangular detector's grid, but false here:
    x_proj/y_proj individually rise and fall as phi sweeps a full circle, so that inference is
    unreliable exactly at the phi=0/2*pi seam a helical/vortex field pattern needs to be
    continuous across. Passing explicit corners (shading='flat') sidesteps the inference entirely.
    """
    radius, _ = read_config_value('spherical_detector_radius', config_path)  # assumed given in 'lambda',
    radius = float(radius)                                                  # matching detector_axes_unit

    N_theta = int(read_config_value('spherical_detector_N_theta', config_path)[0])
    N_phi = int(read_config_value('spherical_detector_N_phi', config_path)[0])
    theta_min = _angle_radians(*read_config_value('spherical_detector_theta_min', config_path))
    theta_max = _angle_radians(*read_config_value('spherical_detector_theta_max', config_path))
    phi_min = _angle_radians(*read_config_value('spherical_detector_phi_min', config_path))
    phi_max = _angle_radians(*read_config_value('spherical_detector_phi_max', config_path))

    cos_theta_min, cos_theta_max = np.cos(theta_min), np.cos(theta_max)
    d_cos_theta = (cos_theta_max - cos_theta_min) / (N_theta - 1) if N_theta > 1 else 0.0
    d_phi = (phi_max - phi_min) / (N_phi - 1) if N_phi > 1 else 0.0

    # Cell centers sit at cos_theta_min + i*d_cos_theta (i=0..N_theta-1) / phi_min + j*d_phi
    # (j=0..N_phi-1) -- see get_row_coordinate/get_col_coordinate -- so the N_theta+1/N_phi+1
    # cell edges are those same centers shifted back by half a step.
    cos_theta_edges = cos_theta_min + (np.arange(N_theta + 1) - 0.5) * d_cos_theta
    cos_theta_edges = np.clip(cos_theta_edges, -1.0, 1.0)
    phi_edges = phi_min + (np.arange(N_phi + 1) - 0.5) * d_phi

    sin_theta_edges = np.sqrt(1.0 - cos_theta_edges**2)
    rho_edges = radius * sin_theta_edges / (1.0 + cos_theta_edges)
    x_edges = rho_edges[:, None] * np.cos(phi_edges)[None, :]
    y_edges = rho_edges[:, None] * np.sin(phi_edges)[None, :]
    return x_edges, y_edges

def get_circular_cell_edges(config_path=DEFAULT_CONFIG_PATH):
    """
    Returns (x_edges, y_edges), each shape (N_R+1, N_phi+1): the lab-plane coordinates of the *cell
    corners* of the circular detector's (N_R, N_phi) annulus grid, mirroring
    Core::Detector::CircularDetector's own 'r*cos(phi)'/'r*sin(phi)' construction, evaluated half a
    grid step beyond each cell center instead of at the centers themselves -- the same
    half-step-back trick get_spherical_cell_edges uses. r is equal-area spaced (r_i^2 linear in i,
    matching CircularDetector's dR_sq step) rather than equal-distance, so the half-step-back is
    taken in r^2, not r.

    Needed for the same reason as get_spherical_cell_edges: pcolormesh's shading='auto'/'nearest'
    infers cell edges from centers assuming they're roughly monotonic in Cartesian x/y, which breaks
    down for any polar (r, phi) grid mapped onto x/y -- most obviously right at the phi=0/2*pi seam.
    """
    N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
    N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
    R_min, _ = read_config_value('circular_detector_R_min', config_path)
    R_max, _ = read_config_value('circular_detector_R_max', config_path)
    R_min, R_max = float(R_min), float(R_max)

    d_R_sq = (R_max**2 - R_min**2) / (N_R - 1) if N_R > 1 else 0.0
    d_phi = (2.0 * np.pi) / (N_phi - 1) if N_phi > 1 else 0.0

    # Cell centers sit at sqrt(R_min^2 + i*d_R_sq) (i=0..N_R-1) / j*d_phi (j=0..N_phi-1) -- see
    # CircularDetector::get_row_coordinate/get_col_coordinate -- so the N_R+1/N_phi+1 cell edges are
    # those same centers shifted back by half a step, taken in r^2 space to match the equal-area
    # spacing (then clipped at 0 before the sqrt, since r^2 can't go negative).
    r_sq_edges = R_min**2 + (np.arange(N_R + 1) - 0.5) * d_R_sq
    r_edges = np.sqrt(np.clip(r_sq_edges, 0.0, None))
    phi_edges = (np.arange(N_phi + 1) - 0.5) * d_phi

    x_edges = r_edges[:, None] * np.cos(phi_edges)[None, :]
    y_edges = r_edges[:, None] * np.sin(phi_edges)[None, :]
    return x_edges, y_edges

def get_screen_coordinates(radiation_filepath, config_path=DEFAULT_CONFIG_PATH):
    """
    Returns (x, y, x_label, y_label): one screen-plane coordinate pair per
    i_screen row of radiation_field.dat, in the same row order.

    - spherical detector: reuses the stereographic projection already exported
      to detector_stereographic.dat next to radiation_field.dat. Row order
      matches i_screen exactly, since both files are written by iterating the
      same underlying Detector_2D::points in order.
    - rectangular detector: reconstructed directly from the config's grid
      bounds, mirroring Core::Detector::RectangularDetector's own
      'x_min + i * dx' construction (a plain linspace).
    - circular detector: reconstructed from the config's (R_min, R_max, N_R,
      N_phi) annulus parameters, mirroring Core::Detector::CircularDetector's
      own 'r*cos(phi)'/'r*sin(phi)' construction (phi spans [0, 2*pi) as a
      plain linspace, matching the C++ side's inclusive-both-ends convention).
    """
    detector_type, _ = read_config_value('detector_type', config_path)

    if detector_type == 'spherical':
        stereo_path = os.path.join(os.path.dirname(radiation_filepath), 'detector_stereographic.dat')
        stereo = pd.read_csv(stereo_path, sep=' ', comment='#')
        return stereo['x_proj'].to_numpy(), stereo['y_proj'].to_numpy(), '$x_{proj}$', '$y_{proj}$'

    elif detector_type == 'rectangular':
        Nx, _ = read_config_value('rectangular_detector_Nx', config_path)
        Ny, _ = read_config_value('rectangular_detector_Ny', config_path)
        x_min, x_unit = read_config_value('rectangular_detector_x_min', config_path)
        x_max, _ = read_config_value('rectangular_detector_x_max', config_path)
        y_min, y_unit = read_config_value('rectangular_detector_y_min', config_path)
        y_max, _ = read_config_value('rectangular_detector_y_max', config_path)

        x_vals = np.linspace(float(x_min), float(x_max), int(Nx))
        y_vals = np.linspace(float(y_min), float(y_max), int(Ny))
        # Row-major (i outer, j inner), matching RectangularDetector's own point order.
        x_grid, y_grid = np.meshgrid(x_vals, y_vals, indexing='ij')
        x_label = '$x$' + (f' [{x_unit}]' if x_unit else '')
        y_label = '$y$' + (f' [{y_unit}]' if y_unit else '')
        return x_grid.ravel(), y_grid.ravel(), x_label, y_label

    elif detector_type == 'circular':
        N_R, _ = read_config_value('circular_detector_N_R', config_path)
        N_phi, _ = read_config_value('circular_detector_N_phi', config_path)
        R_min, r_unit = read_config_value('circular_detector_R_min', config_path)
        R_max, _ = read_config_value('circular_detector_R_max', config_path)

        N_R = int(N_R)
        N_phi = int(N_phi)
        # Equal-area spacing (r^2 linear in i), matching CircularDetector's own construction.
        r_vals = np.sqrt(np.linspace(float(R_min)**2, float(R_max)**2, N_R))
        phi_vals = np.linspace(0.0, 2.0 * np.pi, N_phi)
        # Row-major (i outer over R, j inner over phi), matching CircularDetector's own point order.
        r_grid, phi_grid = np.meshgrid(r_vals, phi_vals, indexing='ij')
        x_grid = r_grid * np.cos(phi_grid)
        y_grid = r_grid * np.sin(phi_grid)
        x_label = '$x$' + (f' [{r_unit}]' if r_unit else '')
        y_label = '$y$' + (f' [{r_unit}]' if r_unit else '')
        return x_grid.ravel(), y_grid.ravel(), x_label, y_label

    else:
        raise ValueError(f"Unknown detector_type '{detector_type}'")

def get_detector_geometry_label(detector_type, config_path):
    """
    Returns a short 'distance=<value> <unit>' (rectangular/circular) or 'radius=<value> <unit>'
    (spherical) string describing how far the detector screen sits from the beam, read straight
    from the per-run config snapshot -- unit is whatever's in the config (conventionally 'lambda'
    for all three keys, per config/coherent_thomson.cfg).
    """
    if detector_type == 'rectangular':
        value, unit = read_config_value('rectangular_detector_distance', config_path)
        label = 'distance'
    elif detector_type == 'circular':
        value, unit = read_config_value('circular_detector_distance', config_path)
        label = 'distance'
    elif detector_type == 'spherical':
        value, unit = read_config_value('spherical_detector_radius', config_path)
        label = 'radius'
    else:
        raise ValueError(f"Unknown detector_type '{detector_type}'")

    return f"{label}={value}" + (f" {unit}" if unit else "")

def plot_radiation_component(range_type, mu, nu, radiation_filepath):
    """
    Plots the real part, imaginary part, magnitude, and phase (2x2 grid) of
    F^{mu nu} of the requested (long_range/short_range) Faraday tensor over
    the detector screen, one figure per configured frequency.

    Rectangular detector: rendered as a colored scatter over its grid.
    Spherical/circular detector: rendered as a true heatmap (pcolormesh over
    the detector's native (N_theta, N_phi)/(N_R, N_phi) grid) rather than a
    scatter -- a scatter over these points has no notion of which points are
    neighbors, which breaks the azimuthal continuity (phi=0 and phi=2*pi are
    the same physical direction) that any helical/vortex structure in the
    field needs to actually look like a spiral around the vertex.
    """
    prefix = {'long': 'LR', 'short': 'SR'}[range_type]
    re_col = f'{prefix}_F{mu}{nu}_re'
    im_col = f'{prefix}_F{mu}{nu}_im'

    try:
        data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    except Exception as e:
        print(f"Error reading file '{radiation_filepath}': {e}")
        sys.exit(1)

    if not {re_col, im_col}.issubset(data.columns):
        print(f"Error: File must contain columns '{re_col}' and '{im_col}'")
        sys.exit(1)

    # Read back the config Core::IoUtils::copy_config_to_run_directory saved alongside this
    # specific run's .dat files, rather than config/coherent_thomson.cfg -- the project's live
    # config may have been edited since this run (e.g. for a different detector geometry), and
    # every key read below is used to interpret THIS file's data, not to launch a new run.
    config_path = os.path.join(os.path.dirname(radiation_filepath), 'config.cfg')
    if not os.path.exists(config_path):
        print(f"Error: '{config_path}' not found -- this run predates per-run config snapshots "
              "(Core::IoUtils::copy_config_to_run_directory); re-run the solver to regenerate it.")
        sys.exit(1)

    detector_type, _ = read_config_value('detector_type', config_path)
    x, y, x_label, y_label = get_screen_coordinates(radiation_filepath, config_path)
    detector_geometry_label = get_detector_geometry_label(detector_type, config_path)

    # Spherical/circular points are laid out row-major over their native (N_theta, N_phi)/(N_R, N_phi)
    # grid -- see Core::Detector_2D::get_grid_indices -- so the field values below reshape cleanly
    # into that grid; x/y are replaced with the cell-corner coordinates pcolormesh needs (see
    # get_spherical_cell_edges/get_circular_cell_edges) instead of the per-point centers
    # get_screen_coordinates returns. Rectangular detectors keep the scatter rendering below: their
    # grid is already Cartesian-monotonic, so pcolormesh's default edge-inference works fine and a
    # scatter is simpler to keep consistent with get_screen_coordinates' plain linspace.
    grid_shape = None
    if detector_type == 'spherical':
        N_theta = int(read_config_value('spherical_detector_N_theta', config_path)[0])
        N_phi = int(read_config_value('spherical_detector_N_phi', config_path)[0])
        grid_shape = (N_theta, N_phi)
        x, y = get_spherical_cell_edges(config_path)
    elif detector_type == 'circular':
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        grid_shape = (N_R, N_phi)
        x, y = get_circular_cell_edges(config_path)

    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder")
    os.makedirs(png_dir, exist_ok=True)

    component_label = f"F^{{{mu}{nu}}}_{{\\mathrm{{{range_type}}}}}"

    for i_omega, subset in data.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_value = subset['omega'].iloc[0]

        re_values = subset[re_col].to_numpy()
        im_values = subset[im_col].to_numpy()
        abs_values = np.hypot(re_values, im_values)
        phase_values = np.arctan2(im_values, re_values)

        # Re/Im are signed (diverging cmap, auto-scaled per panel); magnitude is
        # non-negative (sequential cmap, floored at 0); phase wraps at +-pi (cyclic cmap).
        panels = [
            (re_values, f"$\\mathrm{{Re}}({component_label})$", 'RdBu_r', None, None),
            (im_values, f"$\\mathrm{{Im}}({component_label})$", 'RdBu_r', None, None),
            (abs_values, f"$|{component_label}|$", 'viridis', 0, None),
            (phase_values, f"$\\arg({component_label})$", 'twilight', -np.pi, np.pi),
        ]

        fig, axes = plt.subplots(2, 2, figsize=(12, 11))
        for ax, (values, title, cmap, vmin, vmax) in zip(axes.flat, panels):
            if grid_shape is not None:
                sc = ax.pcolormesh(x, y, values.reshape(grid_shape), cmap=cmap, vmin=vmin, vmax=vmax)
            else:
                sc = ax.scatter(x, y, c=values, cmap=cmap, s=15, edgecolors='none', vmin=vmin, vmax=vmax)
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(title)
            fig.colorbar(sc, ax=ax)

        fig.suptitle(f"Radiated field, $\\omega$ index {i_omega} ($\\omega/\\omega_1$={omega_value:.4g}), "
                     f"{detector_type} detector, {detector_geometry_label}",
                     fontsize=13, fontweight='bold')
        plt.tight_layout()

        output_name = f"radiation_field_{range_type}_F{mu}{nu}_omega{i_omega}.png"
        output_img = os.path.join(png_dir, output_name)
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: python3 {sys.argv[0]} <long|short> <mu> <nu> [path_to_radiation_field.dat]")
        sys.exit(1)

    range_type = sys.argv[1]
    if range_type not in ('long', 'short'):
        print("Error: first argument must be 'long' or 'short'")
        sys.exit(1)

    try:
        mu = int(sys.argv[2])
        nu = int(sys.argv[3])
    except ValueError:
        print("Error: mu and nu must be integers")
        sys.exit(1)

    if not (0 <= mu <= 3 and 0 <= nu <= 3):
        print("Error: mu and nu must each be in [0, 3]")
        sys.exit(1)

    if len(sys.argv) >= 5:
        input_file = sys.argv[4]
    else:
        try:
            input_file = find_latest_output_file("radiation_field.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            sys.exit(1)
        print(f"No file given; using latest run's radiation field: {input_file}")

    plot_radiation_component(range_type, mu, nu, input_file)
