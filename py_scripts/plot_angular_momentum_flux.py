"""
Computes the spectral electromagnetic angular-momentum flux density along the detector's own
normal direction, per theory/angular_momentum_flux_density.md, entirely in Python from an
already-completed run's radiation_field.dat -- no new C++ output is needed, this reuses the
long-range/short-range Faraday tensor components the solver already exports.

Scope and restrictions (read before trusting a result):
- Only 'rectangular' and 'circular' detectors are supported. The flux formula assumes a flat
  transverse plane with a single, shared normal; a spherical detector's points don't generally lie
  in such a plane (unless collapsed to a single point, in which case use a 1x1 rectangular/circular
  detector instead -- see CLAUDE.md's dense_frequency_spectrum bullet, which explicitly allows any
  detector type collapsed to one point, not just spherical).
- The theory doc derives the flux through a plane whose normal is 'Oz'; here that's generalized to
  "the detector's own local normal" (detector_direction_theta/phi), not necessarily the laser's
  canonical Oz -- see get_canonical_to_local_rotation. The two coincide exactly when
  detector_direction_theta/phi = (0, *) (detector facing the laser) or (pi, *) (facing directly
  backward, the repo's default config) at the theta=0/pi identity/flip degeneracies
  Core::MathUtils::rotation_matrix_from_direction hardcodes; for any other angle the detector plane
  is tilted, and the flux is now that tilted plane's own out-of-plane (z_local) component, not the
  canonical Oz component.
- Requires print_field_in_canonical_frame=true (the config default) OR laser_nx/ny/nz to be
  present so the lab-frame case can be un-rotated -- both are handled below.

Decomposition: the theory doc splits the total field into long-range ('l', the ~1/R radiation
field) and short-range ('s', the ~1/R^2 velocity field) contributions and works out the flux's
long-long ('ll'), long-short ('ls'), and short-long ('sl') cross-terms; the short-short ('ss') term
is the obvious omitted fourth term needed for ll+ls+sl+ss to equal the flux computed directly from
the physical total field (E_l+E_s, B_l+B_s) -- included here for that consistency, not because the
doc names it.
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file
from plot_radiation_field import (
    read_config_value,
    get_screen_coordinates,
    get_rectangular_cell_edges,
    get_circular_cell_edges,
    get_detector_geometry_label,
)

# Core::PhysUtils::AtomicUnits (phys_utils.hpp): c and epsilon_0 in the solver's own atomic units.
C_LIGHT = 137.036
EPSILON_0 = 1.0 / (4.0 * np.pi)


def convert_unit_to_number(unit, config_path):
    """Python mirror of Core::IoUtils::convert_unit_to_number (io_utils.hpp) -- only the branches
    detector-geometry/direction keys actually use ('lambda', 'pi', 'a.u.', bare)."""
    unit = unit.lower()
    if unit == 'pi':
        return np.pi
    if unit == 'lambda':
        omega = float(read_config_value('laser_frequency', config_path)[0])
        return 2.0 * np.pi / omega * C_LIGHT
    if unit == 'mc':
        return C_LIGHT
    if unit == 'cycles_adim':
        return 2.0 * np.pi
    if unit == 'omega_laser':
        return float(read_config_value('laser_frequency', config_path)[0])
    if unit in ('a.u.', ''):
        return 1.0
    print(f"Warning: Unknown unit '{unit}'. Assuming 1.0.")
    return 1.0


def rotation_matrix_from_direction(nx, ny, nz):
    """
    Python mirror of Core::MathUtils::rotation_matrix_from_direction (math_utils.hpp): the 3x3
    rotation mapping local +z onto the given direction, including its two hardcoded near-degenerate
    cases (exactly +z / exactly -z), needed here bit-for-bit since the angular-momentum flux is not
    invariant under an arbitrary choice within those degeneracies (unlike a plain axis alignment).
    """
    norm = np.sqrt(nx * nx + ny * ny + nz * nz)
    if norm < 1e-12:
        return np.eye(3)
    nx, ny, nz = nx / norm, ny / norm, nz / norm
    if nz > 0.999999:
        return np.eye(3)
    if nz < -0.999999:
        return np.diag([1.0, -1.0, -1.0])
    k = 1.0 / (1.0 + nz)
    return np.array([
        [1.0 - nx * nx * k, -nx * ny * k, nx],
        [-nx * ny * k, 1.0 - ny * ny * k, ny],
        [-nx, -ny, nz],
    ])


def get_detector_local_rotation(config_path):
    """
    Returns the 3x3 rotation R such that v_canonical = R @ v_local -- a Python mirror of
    Core::Detector::Detector_2D's own local_rotation, built from detector_direction_theta/phi.
    """
    theta_str, theta_unit = read_config_value('detector_direction_theta', config_path)
    phi_str, phi_unit = read_config_value('detector_direction_phi', config_path)
    theta = float(theta_str) * convert_unit_to_number(theta_unit, config_path)
    phi = float(phi_str) * convert_unit_to_number(phi_unit, config_path)
    dir_x, dir_y, dir_z = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)
    return rotation_matrix_from_direction(dir_x, dir_y, dir_z)


def get_laser_lab_rotation(config_path):
    """
    Returns the 3x3 rotation R such that v_lab = R @ v_canonical -- a Python mirror of
    Core::Laser::LaserField's own rotation_matrix, built from laser_nx/ny/nz.
    """
    nx = float(read_config_value('laser_nx', config_path)[0])
    ny = float(read_config_value('laser_ny', config_path)[0])
    nz = float(read_config_value('laser_nz', config_path)[0])
    return rotation_matrix_from_direction(nx, ny, nz)


def get_field_to_canonical_rotation(config_path):
    """
    Returns the 3x3 rotation R such that v_canonical = R @ v_field, where v_field is however
    radiation_field.dat's Faraday tensor is actually stored: the canonical frame itself (identity)
    if this run's config has print_field_in_canonical_frame=true (the default), or the lab frame
    otherwise, in which case the laser's own rotation (get_laser_lab_rotation) is undone.
    """
    canonical_frame = read_config_value('print_field_in_canonical_frame', config_path)[0]
    if canonical_frame.lower() == 'true':
        return np.eye(3)
    return get_laser_lab_rotation(config_path).T


def get_canonical_to_local_rotation(config_path):
    """
    Returns the 3x3 rotation R such that v_local = R @ v_field (see get_field_to_canonical_rotation
    for what frame v_field is actually in), undoing Core::Detector::Detector_2D's own
    local_rotation on top of that -- needed because this script's screen coordinates (x_local,
    y_local from get_screen_coordinates_local_au) are the detector's own local, pre-rotation
    coordinates, while radiation_field.dat's Faraday tensor is not.
    """
    return get_detector_local_rotation(config_path).T @ get_field_to_canonical_rotation(config_path)


def get_screen_coordinates_local_au(detector_type, config_path):
    """
    Returns (x_local, y_local): the detector's own local-frame transverse coordinates of each
    screen point, in atomic length units and in the same row-major i_screen order as
    radiation_field.dat -- unlike plot_radiation_field.get_screen_coordinates (display-only, left
    in raw config-file units for readable axis ticks), these are unit-converted, since they get
    combined dimensionally with Faraday tensor field values (already in atomic units) below.
    """
    if detector_type == 'rectangular':
        Nx = int(read_config_value('rectangular_detector_Nx', config_path)[0])
        Ny = int(read_config_value('rectangular_detector_Ny', config_path)[0])
        x_min_s, x_unit = read_config_value('rectangular_detector_x_min', config_path)
        x_max_s, _ = read_config_value('rectangular_detector_x_max', config_path)
        y_min_s, y_unit = read_config_value('rectangular_detector_y_min', config_path)
        y_max_s, _ = read_config_value('rectangular_detector_y_max', config_path)
        x_min = float(x_min_s) * convert_unit_to_number(x_unit, config_path)
        x_max = float(x_max_s) * convert_unit_to_number(x_unit, config_path)
        y_min = float(y_min_s) * convert_unit_to_number(y_unit, config_path)
        y_max = float(y_max_s) * convert_unit_to_number(y_unit, config_path)

        x_vals = np.linspace(x_min, x_max, Nx)
        y_vals = np.linspace(y_min, y_max, Ny)
        x_grid, y_grid = np.meshgrid(x_vals, y_vals, indexing='ij')
        return x_grid.ravel(), y_grid.ravel()

    elif detector_type == 'circular':
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        R_min_s, r_unit = read_config_value('circular_detector_R_min', config_path)
        R_max_s, _ = read_config_value('circular_detector_R_max', config_path)
        R_min = float(R_min_s) * convert_unit_to_number(r_unit, config_path)
        R_max = float(R_max_s) * convert_unit_to_number(r_unit, config_path)

        r_vals = np.sqrt(np.linspace(R_min**2, R_max**2, N_R))
        phi_vals = np.linspace(0.0, 2.0 * np.pi, N_phi)
        r_grid, phi_grid = np.meshgrid(r_vals, phi_vals, indexing='ij')
        return (r_grid * np.cos(phi_grid)).ravel(), (r_grid * np.sin(phi_grid)).ravel()

    else:
        raise ValueError(
            f"detector_type '{detector_type}' is not a flat screen with a single well-defined "
            "local (x, y) plane -- the angular-momentum flux density formula (theory/"
            "angular_momentum_flux_density.md) assumes flux through one flat transverse plane, "
            "which only 'rectangular' and 'circular' detectors are. For a single-point spectrum "
            "(the dense_frequency_spectrum workflow), use a 1x1 rectangular or circular detector "
            "instead of a spherical one."
        )


def extract_rotated_faraday_fields(data, R):
    """
    Returns a dict of six complex numpy arrays per range (Ex/Ey/Ez/Bx/By/Bz, suffixed '_l' for
    long-range and '_s' for short-range), read out of radiation_field.dat's F^{mu nu} columns and
    rotated by R into whatever target frame the caller needs (this script's callers pass
    get_canonical_to_local_rotation's result to land in the detector's own local frame;
    plot_spherical_field_components.py instead passes get_field_to_canonical_rotation's result to
    land in the canonical frame).

    Column mapping mirrors Core::Laser::LaserField::get_faraday_tensor's F^{mu nu} <-> E/B sign
    convention (F^{i0}=E_i, F^{jk}=-eps_jkl B_l): Ex=F10, Ey=F20, Ez=F30, Bx=F32, By=F13, Bz=F21.
    """
    fields = {}
    for prefix, suffix in (('LR', 'l'), ('SR', 's')):
        def col(mu, nu):
            return (data[f'{prefix}_F{mu}{nu}_re'] + 1j * data[f'{prefix}_F{mu}{nu}_im']).to_numpy()

        E_stored = np.stack([col(1, 0), col(2, 0), col(3, 0)])
        B_stored = np.stack([col(3, 2), col(1, 3), col(2, 1)])
        Ex, Ey, Ez = R @ E_stored
        Bx, By, Bz = R @ B_stored
        fields[f'Ex_{suffix}'], fields[f'Ey_{suffix}'], fields[f'Ez_{suffix}'] = Ex, Ey, Ez
        fields[f'Bx_{suffix}'], fields[f'By_{suffix}'], fields[f'Bz_{suffix}'] = Bx, By, Bz
    return fields


def angular_momentum_flux_density(x, y, Ex, Ey, Ez, Bx, By, Bz):
    """
    d(F_Jz)/d(omega) = (epsilon_0/pi) * Re[x*(Ez*Bx* - Ex*Bz*) - y*(Ey*Bz* - Ez*By*)]
    (theory/angular_momentum_flux_density.md, section 2), vectorized over any matching-shape
    complex field-component arrays and real local-frame position arrays.
    """
    return (EPSILON_0 / np.pi) * np.real(
        x * (Ez * np.conj(Bx) - Ex * np.conj(Bz)) - y * (Ey * np.conj(Bz) - Ez * np.conj(By))
    )


def compute_angular_momentum_flux(radiation_filepath):
    """
    Returns (result, detector_type, config_path). `result` is a DataFrame with one row per
    radiation_field.dat row (i_omega, omega, i_screen), plus the local-frame (x_local, y_local)
    position and the flux density's long-long/long-short/short-long/short-short decomposition and
    total (their sum, equal to what plugging in the full field E_l+E_s, B_l+B_s would give).
    """
    config_path = os.path.join(os.path.dirname(radiation_filepath), 'config.cfg')
    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"'{config_path}' not found -- this run predates per-run config snapshots "
            "(Core::IoUtils::copy_config_to_run_directory); re-run the solver to regenerate it.")

    detector_type = read_config_value('detector_type', config_path)[0]
    x_local, y_local = get_screen_coordinates_local_au(detector_type, config_path)

    data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    x = x_local[data['i_screen'].to_numpy()]
    y = y_local[data['i_screen'].to_numpy()]

    R_to_local = get_canonical_to_local_rotation(config_path)
    f = extract_rotated_faraday_fields(data, R_to_local)

    flux_ll = angular_momentum_flux_density(x, y, f['Ex_l'], f['Ey_l'], f['Ez_l'], f['Bx_l'], f['By_l'], f['Bz_l'])
    flux_ls = angular_momentum_flux_density(x, y, f['Ex_l'], f['Ey_l'], f['Ez_l'], f['Bx_s'], f['By_s'], f['Bz_s'])
    flux_sl = angular_momentum_flux_density(x, y, f['Ex_s'], f['Ey_s'], f['Ez_s'], f['Bx_l'], f['By_l'], f['Bz_l'])
    flux_ss = angular_momentum_flux_density(x, y, f['Ex_s'], f['Ey_s'], f['Ez_s'], f['Bx_s'], f['By_s'], f['Bz_s'])

    result = data[['i_omega', 'omega', 'i_screen']].copy()
    result['x_local'] = x
    result['y_local'] = y
    result['flux_ll'] = flux_ll
    result['flux_ls'] = flux_ls
    result['flux_sl'] = flux_sl
    result['flux_ss'] = flux_ss
    result['flux_total'] = flux_ll + flux_ls + flux_sl + flux_ss
    return result, detector_type, config_path


def plot_angular_momentum_flux(radiation_filepath):
    """
    Plots all four of the theory doc's flux-density terms -- total, long-range only (ll), and the
    two long/short cross-terms (ls, sl) -- computed by compute_angular_momentum_flux, in one of two
    shapes depending on the detector's point count, mirroring the two shapes
    plot_radiation_field.py/plot_point_spectrum.py already use for the Faraday tensor itself:
    - a single screen point (a 1x1 rectangular/circular detector, the dense_frequency_spectrum
      workflow): one line plot vs. omega, all four terms overlaid.
    - multiple screen points: one 2x2-panel heatmap PNG per frequency, rendered on the detector's
      own native grid like plot_radiation_field.py's pcolormesh panels.
    """
    result, detector_type, config_path = compute_angular_momentum_flux(radiation_filepath)

    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder")
    os.makedirs(png_dir, exist_ok=True)

    # The four terms named in the theory doc: total (the full field), long-range only (ll), and
    # the two long/short cross-terms (ls, sl) -- the doc's own decomposition, short of the
    # implicit (ss) term folded into 'total' but not surfaced as its own panel/line.
    terms = (
        ('flux_total', 'Total (long+short)'),
        ('flux_ll', 'Long-range only'),
        ('flux_ls', 'Long E x Short B'),
        ('flux_sl', 'Short E x Long B'),
    )

    if result['i_screen'].nunique() == 1:
        subset = result.sort_values('omega')
        fig, ax = plt.subplots(figsize=(8, 5))
        for col, label in terms:
            ax.plot(subset['omega'], subset[col], marker='.', label=label)
        ax.set_xlabel("$\\omega / \\omega_1$ (units of the fundamental)")
        ax.set_ylabel("$d\\mathcal{F}_{J_z}/d\\omega$ (a.u.)")
        ax.grid(True)
        ax.legend()
        fig.suptitle("Spectral angular-momentum flux density", fontsize=13, fontweight='bold')
        plt.tight_layout()

        output_img = os.path.join(png_dir, "angular_momentum_flux_spectrum.png")
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)
        return

    if detector_type == 'rectangular':
        Nx = int(read_config_value('rectangular_detector_Nx', config_path)[0])
        Ny = int(read_config_value('rectangular_detector_Ny', config_path)[0])
        grid_shape = (Nx, Ny)
        x_edges, y_edges = get_rectangular_cell_edges(config_path)
    else:
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        grid_shape = (N_R, N_phi)
        x_edges, y_edges = get_circular_cell_edges(config_path)

    _, _, x_label, y_label = get_screen_coordinates(radiation_filepath, config_path)
    detector_geometry_label = get_detector_geometry_label(detector_type, config_path)

    for i_omega, subset in result.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_value = subset['omega'].iloc[0]

        fig, axes = plt.subplots(2, 2, figsize=(13, 11))
        for ax, (col, title) in zip(axes.flat, terms):
            values = subset[col].to_numpy().reshape(grid_shape)
            vmax = np.max(np.abs(values)) if np.any(values) else 1.0
            sc = ax.pcolormesh(x_edges, y_edges, values, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(title)
            fig.colorbar(sc, ax=ax)

        fig.suptitle(f"$d\\mathcal{{F}}_{{J_z}}/d\\omega$, $\\omega$ index {i_omega} "
                     f"($\\omega/\\omega_1$={omega_value:.4g}), {detector_type} detector, "
                     f"{detector_geometry_label}", fontsize=12, fontweight='bold')
        plt.tight_layout()

        output_img = os.path.join(png_dir, f"angular_momentum_flux_omega{i_omega}.png")
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) >= 2:
        input_file = sys.argv[1]
    else:
        try:
            input_file = find_latest_output_file("radiation_field.dat")
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            sys.exit(1)
        print(f"No file given; using latest run's radiation field: {input_file}")

    try:
        plot_angular_momentum_flux(input_file)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)
