"""
Computes and plots the z-component of the spin angular momentum spectral density,
$d{\\cal S}_z/d\\omega = (4/\\omega)\\,\\mathrm{Im}[\\tilde E_x^*\\tilde E_y]$, from a run's
coherently-summed Faraday tensor (radiation_field.dat) or the incident beam's own analytic field
(incident_field.dat, via --incident) -- see
theory/numerical_calculation_of_angular_momentum.md for the derivation (its "Primary definition of
OAM and SAM (Oz components) densities" section; only the S_z piece is implemented here, not L_z).

**Always uses the total (physical) Faraday tensor, LR+SR+BR summed** -- unlike plot_field.py's own
'total'/'long'/'short'/'boundary' choice, there is no long-range/short-range/boundary decomposition
of an angular-momentum observable: that split only has meaning for the derivation's intermediate
F^{mu nu} terms (see theory/FT_Faraday_tensor-direct_and_simplified_forms.md), and only their sum is
the actual radiated field E appearing in the S_z formula above. BR is treated as zero if the file
predates the boundary term; for incident_field.dat, SR/BR are already identically zero so the sum
reduces to LR alone.

Ex/Ey are recovered directly from the stored F^{10}/F^{20} columns (Ex=c*F10, Ey=c*F20 -- the same
F^{mu nu}<->E/B convention documented in faraday_frame_utils.py's extract_rotated_faraday_fields),
scaled by C_LIGHT.

**Frequency convention**: radiation_field.dat/incident_field.dat only ever export
omega/omega_1 (dimensionless, already normalized by the fundamental -- see
Radiation::plot_radiation_field's own 'omega' column comment); no Python script in this repo
recovers the raw atomic-unit omega from a .dat file (it would require re-deriving
fundamental_frequency, which is only written out, human-readable, in run_log.txt). So the
$1/\\omega$ prefactor above uses that same dimensionless ratio, making this
$d{\\cal S}_z/d(\\omega/\\omega_1)$ rather than $d{\\cal S}_z/d\\omega$ in raw atomic units --
consistent with every other frequency-axis quantity already plotted by this project (e.g.
plot_point_spectrum.py), and differing from the true $d{\\cal S}_z/d\\omega$ only by the constant
factor $\\omega_1$ (fundamental_frequency), which does not affect the pattern's shape at any given
$\\omega$ index or its relative comparison across harmonics.

Reuses plot_field.py's screen-geometry/cell-edge reconstruction directly (same CLI shape:
[path_to_run_folder] [--incident]), and the same is_incident split into
png_folder/radiation/{emitted,incident}/ subfolders, so a real run's spin-density plots and the
incident beam's never collide and sit alongside the Faraday-tensor component plots.
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.run_output_utils import find_latest_output_file
from plot_field import (
    read_config_value,
    get_screen_coordinates,
    get_rectangular_cell_edges,
    get_spherical_plot_grid,
    get_circular_cell_edges,
    get_detector_geometry_label,
)
from faraday_frame_utils import C_LIGHT
from utils.w0_axes_utils import get_laser_lg_w0_in_axes_units, add_w0_secondary_axes

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))


def _complex_column(data, prefix, mu, nu):
    return (data[f'{prefix}_F{mu}{nu}_re'] + 1j * data[f'{prefix}_F{mu}{nu}_im']).to_numpy()


def _get_Ex_Ey(data):
    """
    Returns (Ex, Ey), each a complex numpy array (one entry per row of `data`), from the total
    (LR+SR+BR) Faraday tensor -- see this module's own doc comment for why only the total, never a
    single range, is used here. Ex=c*F10, Ey=c*F20.
    """
    Ex = _complex_column(data, 'LR', 1, 0) + _complex_column(data, 'SR', 1, 0)
    Ey = _complex_column(data, 'LR', 2, 0) + _complex_column(data, 'SR', 2, 0)
    if {'BR_F10_re', 'BR_F10_im'}.issubset(data.columns):
        Ex = Ex + _complex_column(data, 'BR', 1, 0)
        Ey = Ey + _complex_column(data, 'BR', 2, 0)
    return C_LIGHT * Ex, C_LIGHT * Ey


def plot_spin_angular_momentum(radiation_filepath):
    """
    Plots $d{\\cal S}_z/d(\\omega/\\omega_1)$ (single-panel heatmap, diverging colormap centered on
    zero since the quantity is real-valued -- unlike plot_field.py's Re/Im/Abs/Phase 2x2 grid, there
    is nothing else to show) over the detector screen, one figure per configured frequency, computed
    from the total Faraday tensor only (see this module's own doc comment).

    `radiation_filepath` may point at either a run's radiation_field.dat or its incident_field.dat
    sibling, detected by filename -- see this module's own doc comment for how each is handled.
    """
    is_incident = os.path.basename(radiation_filepath) == "incident_field.dat"
    title_suffix = " (incident beam)" if is_incident else ""

    try:
        data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    except Exception as e:
        print(f"Error reading file '{radiation_filepath}': {e}")
        sys.exit(1)

    required = {f'{p}_F{m}{n}_re' for p in ('LR', 'SR') for m, n in ((1, 0), (2, 0))} | \
        {f'{p}_F{m}{n}_im' for p in ('LR', 'SR') for m, n in ((1, 0), (2, 0))}
    if not required.issubset(data.columns):
        print(f"Error: File must contain columns {sorted(required)}")
        sys.exit(1)

    # Read back the config Core::IoUtils::copy_config_to_run_directory saved alongside this
    # specific run's .dat files, rather than config/coherent_thomson.cfg -- see plot_field.py's
    # plot_radiation_component for the full rationale.
    config_path = os.path.join(os.path.dirname(radiation_filepath), 'config.cfg')
    if not os.path.exists(config_path):
        print(f"Error: '{config_path}' not found -- this run predates per-run config snapshots "
              "(Core::IoUtils::copy_config_to_run_directory); re-run the solver to regenerate it.")
        sys.exit(1)

    detector_type, _ = read_config_value('detector_type', config_path)
    _, _, x_label, y_label = get_screen_coordinates(radiation_filepath, config_path)
    detector_geometry_label = get_detector_geometry_label(detector_type, config_path)

    # Same cell-edge/grid-shape reconstruction as plot_field.py's plot_radiation_component -- see
    # that function's own doc comment for why pcolormesh needs explicit cell corners for all three
    # detector types.
    aspect = 'equal'
    axes_unit = None
    if detector_type == 'rectangular':
        Nx = int(read_config_value('rectangular_detector_Nx', config_path)[0])
        Ny = int(read_config_value('rectangular_detector_Ny', config_path)[0])
        grid_shape = (Nx, Ny)
        x, y = get_rectangular_cell_edges(config_path)
        axes_unit = read_config_value('rectangular_detector_x_min', config_path)[1]
    elif detector_type == 'spherical':
        N_theta = int(read_config_value('spherical_detector_N_theta', config_path)[0])
        N_phi = int(read_config_value('spherical_detector_N_phi', config_path)[0])
        grid_shape = (N_theta, N_phi)
        x, y, x_label, y_label, aspect = get_spherical_plot_grid(config_path)
        if aspect == 'equal':
            axes_unit = read_config_value('spherical_detector_radius', config_path)[1]
    else:
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        grid_shape = (N_R, N_phi)
        x, y = get_circular_cell_edges(config_path)
        axes_unit = read_config_value('circular_detector_R_min', config_path)[1]

    w0 = get_laser_lg_w0_in_axes_units(radiation_filepath, axes_unit) if axes_unit else None

    # Own 'emitted'/'incident' subfolders of png_folder/radiation/, exactly mirroring
    # plot_field.py's plot_radiation_component split -- so both sets of plots (Faraday tensor
    # components and spin density) sit side by side under the same two subfolders, distinguished
    # only by filename.
    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder", MODULE_NAME,
                           "incident" if is_incident else "emitted")
    os.makedirs(png_dir, exist_ok=True)

    for i_omega, subset in data.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_ratio = subset['omega'].iloc[0]

        Ex, Ey = _get_Ex_Ey(subset)
        if omega_ratio == 0.0:
            print(f"Warning: omega/omega_1 == 0 at i_omega={i_omega}; skipping (division by zero).")
            continue
        dSz = (4.0 / omega_ratio) * np.imag(np.conj(Ex) * Ey)

        abs_max = np.abs(dSz).max() if dSz.size else None

        fig, ax = plt.subplots(figsize=(6.5, 5.5), layout='constrained')
        sc = ax.pcolormesh(x, y, dSz.reshape(grid_shape), cmap='RdBu_r', vmin=-abs_max, vmax=abs_max)
        ax.set_aspect(aspect, adjustable='box')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title("$d{\\cal S}_z/d(\\omega/\\omega_1)$")
        ax.grid(True, alpha=0.25)
        fig.colorbar(sc, ax=ax, shrink=0.85)
        add_w0_secondary_axes(ax, w0)

        fig.suptitle(f"Spin angular momentum density{title_suffix}, $\\omega$ index {i_omega} "
                     f"($\\omega/\\omega_1$={omega_ratio:.4g}), {detector_type} detector, "
                     f"{detector_geometry_label}", fontsize=11)

        output_name = f"spin_angular_momentum_z_omega{i_omega}.png"
        output_img = os.path.join(png_dir, output_name)
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)


if __name__ == "__main__":
    args = sys.argv[1:]
    use_incident = "--incident" in args
    if use_incident:
        args.remove("--incident")

    field_filename = "incident_field.dat" if use_incident else "radiation_field.dat"
    if len(args) >= 1:
        input_file = os.path.join(args[0], field_filename)
        if not os.path.exists(input_file):
            print(f"Error: '{input_file}' not found")
            sys.exit(1)
    else:
        try:
            input_file = find_latest_output_file(field_filename)
        except (ValueError, FileNotFoundError) as e:
            print(f"Usage error: {e}")
            sys.exit(1)
        print(f"No folder given; using latest run: {input_file}")

    plot_spin_angular_momentum(input_file)
