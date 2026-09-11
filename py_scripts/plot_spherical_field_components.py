"""
Computes the canonical-frame (laser along Oz) spherical components of the coherently summed
Faraday tensor -- E_r, E_theta, E_phi, B_r, B_theta, B_phi -- for a spherical-detector run's
radiation_field.dat. Python-only post-processing, no new C++ output needed: reuses the rotation
machinery plot_angular_momentum_flux.py already built (detector/laser direction reconstruction,
Faraday-tensor <-> Cartesian E/B column mapping).

A spherical detector is the natural fit for this (unlike the flat rectangular/circular detectors
plot_angular_momentum_flux.py is restricted to): every screen point already has its own natural
observation direction (theta, phi), so there's no shared-plane assumption to violate. Each point's
own local (theta_local, phi_local) -- Core::Detector::SphericalDetector's own cone-point
construction, generally relative to the detector's own axis, not necessarily canonical Oz -- is
first rotated by the detector's local_rotation (detector_direction_theta/phi) into a canonical
observation direction (theta_c, phi_c); the field (already rotated into the canonical frame, or
rotated there from the lab frame if print_field_in_canonical_frame=false) is then projected onto
the standard orthonormal (r_hat, theta_hat, phi_hat) basis at that direction.
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from run_output_utils import find_latest_output_file
from plot_radiation_field import read_config_value, get_spherical_plot_grid, get_detector_geometry_label
from plot_angular_momentum_flux import (
    convert_unit_to_number,
    get_detector_local_rotation,
    get_field_to_canonical_rotation,
    extract_rotated_faraday_fields,
)
from w0_axes_utils import get_laser_lg_w0_in_axes_units, add_w0_secondary_axes


def get_spherical_local_angles(config_path):
    """
    Returns (theta_local, phi_local): each screen point's own local (pre-rotation) spherical
    angles, in the same row-major i_screen order as radiation_field.dat -- mirrors
    Core::Detector::SphericalDetector's own 'cos_theta_cone_min + i*d_cos_theta_cone' /
    'phi_cone_min + j*d_phi_cone' construction.
    """
    N_theta = int(read_config_value('spherical_detector_N_theta', config_path)[0])
    N_phi = int(read_config_value('spherical_detector_N_phi', config_path)[0])
    theta_min_s, theta_unit = read_config_value('spherical_detector_theta_min', config_path)
    theta_max_s, _ = read_config_value('spherical_detector_theta_max', config_path)
    phi_min_s, phi_unit = read_config_value('spherical_detector_phi_min', config_path)
    phi_max_s, _ = read_config_value('spherical_detector_phi_max', config_path)
    theta_min = float(theta_min_s) * convert_unit_to_number(theta_unit, config_path)
    theta_max = float(theta_max_s) * convert_unit_to_number(theta_unit, config_path)
    phi_min = float(phi_min_s) * convert_unit_to_number(phi_unit, config_path)
    phi_max = float(phi_max_s) * convert_unit_to_number(phi_unit, config_path)

    cos_theta_min, cos_theta_max = np.cos(theta_min), np.cos(theta_max)
    d_cos_theta = (cos_theta_max - cos_theta_min) / (N_theta - 1) if N_theta > 1 else 0.0
    d_phi = (phi_max - phi_min) / (N_phi - 1) if N_phi > 1 else 0.0

    cos_theta_vals = cos_theta_min + np.arange(N_theta) * d_cos_theta
    theta_vals = np.arccos(np.clip(cos_theta_vals, -1.0, 1.0))
    phi_vals = phi_min + np.arange(N_phi) * d_phi

    theta_grid, phi_grid = np.meshgrid(theta_vals, phi_vals, indexing='ij')
    return theta_grid.ravel(), phi_grid.ravel()


def get_canonical_observation_angles(config_path):
    """
    Returns (theta_c, phi_c): each screen point's observation direction expressed in the canonical
    frame (laser along Oz), by rotating its own local direction (get_spherical_local_angles) with
    the detector's own local_rotation (get_detector_local_rotation) -- equal to the local angles
    only when detector_direction_theta/phi places the detector's axis exactly along canonical
    +-Oz (the two cases Core::MathUtils::rotation_matrix_from_direction special-cases).
    """
    theta_local, phi_local = get_spherical_local_angles(config_path)
    n_local = np.stack([
        np.sin(theta_local) * np.cos(phi_local),
        np.sin(theta_local) * np.sin(phi_local),
        np.cos(theta_local),
    ])
    n_canonical = get_detector_local_rotation(config_path) @ n_local
    theta_c = np.arccos(np.clip(n_canonical[2], -1.0, 1.0))
    phi_c = np.arctan2(n_canonical[1], n_canonical[0])
    return theta_c, phi_c


def spherical_unit_vectors(theta, phi):
    """
    Returns (r_hat, theta_hat, phi_hat), each shape (3, N): the standard right-handed orthonormal
    spherical basis at each (theta, phi) (r_hat matches
    Core::MathUtils::create_unit_light_like_vector(theta, phi)'s own direction convention).
    """
    r_hat = np.stack([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    theta_hat = np.stack([np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)])
    phi_hat = np.stack([-np.sin(phi), np.cos(phi), np.zeros_like(phi)])
    return r_hat, theta_hat, phi_hat


def compute_spherical_field_components(radiation_filepath):
    """
    Returns (result, config_path). `result` is a DataFrame with one row per radiation_field.dat row
    (i_omega, omega, i_screen), plus six complex columns per range -- 'LR_Er'/'LR_Etheta'/
    'LR_Ephi'/'LR_Br'/'LR_Btheta'/'LR_Bphi' and the 'SR_'/'BR_' equivalents -- the canonical-frame
    spherical components of the long-range/short-range/boundary Faraday tensor. Unlike
    plot_angular_momentum_flux.py's compute_angular_momentum_flux, this function makes no
    long+short 'total' assumption that boundary would invalidate -- it just projects each range's
    tensor independently, so adding 'BR' here is a plain three-way extension of the existing
    two-way loop below.
    """
    config_path = os.path.join(os.path.dirname(radiation_filepath), 'config.cfg')
    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"'{config_path}' not found -- this run predates per-run config snapshots "
            "(Core::IoUtils::copy_config_to_run_directory); re-run the solver to regenerate it.")

    detector_type = read_config_value('detector_type', config_path)[0]
    if detector_type != 'spherical':
        raise ValueError(
            f"detector_type '{detector_type}' is not 'spherical' -- spherical field components "
            "are only meaningful for a spherical detector, where every screen point already has "
            "its own natural observation direction to project the field onto."
        )

    theta_c, phi_c = get_canonical_observation_angles(config_path)
    r_hat, theta_hat, phi_hat = spherical_unit_vectors(theta_c, phi_c)

    data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    i_screen = data['i_screen'].to_numpy()
    r_hat, theta_hat, phi_hat = r_hat[:, i_screen], theta_hat[:, i_screen], phi_hat[:, i_screen]

    R_to_canonical = get_field_to_canonical_rotation(config_path)
    f = extract_rotated_faraday_fields(data, R_to_canonical)

    result = data[['i_omega', 'omega', 'i_screen']].copy()
    for prefix, suffix in (('LR', 'l'), ('SR', 's'), ('BR', 'b')):
        E = np.stack([f[f'Ex_{suffix}'], f[f'Ey_{suffix}'], f[f'Ez_{suffix}']])
        B = np.stack([f[f'Bx_{suffix}'], f[f'By_{suffix}'], f[f'Bz_{suffix}']])
        for name, vec in (('E', E), ('B', B)):
            result[f'{prefix}_{name}r'] = np.sum(vec * r_hat, axis=0)
            result[f'{prefix}_{name}theta'] = np.sum(vec * theta_hat, axis=0)
            result[f'{prefix}_{name}phi'] = np.sum(vec * phi_hat, axis=0)

    return result, config_path


def plot_spherical_field_component(range_type, field, component, radiation_filepath):
    """
    Plots the real part, imaginary part, magnitude, and phase (2x2 grid) of the requested
    canonical-frame spherical field component (field in {'E', 'B'}, component in
    {'r', 'theta', 'phi'}, long-range, short-range, or boundary) -- the spherical-component
    counterpart of plot_radiation_field.py's plot_radiation_component / plot_point_spectrum.py's
    plot_point_spectrum (which instead read Cartesian F^{mu nu} components directly).

    'boundary' is identically zero when the run used radiation_formula="direct" -- see
    theory/FT_Faraday_tensor-direct_and_simplified_forms.md's "Form 2's boundary term F_b" section.

    Renders as a heatmap per frequency (on the spherical detector's own stereographic-projection
    grid, like plot_radiation_field.py) if there's more than one screen point, or a line plot vs.
    omega (like plot_point_spectrum.py) for a single-point (dense_frequency_spectrum-style) run.
    """
    prefix = {'long': 'LR', 'short': 'SR', 'boundary': 'BR'}[range_type]
    col = f'{prefix}_{field}{component}'

    result, config_path = compute_spherical_field_components(radiation_filepath)

    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder")
    os.makedirs(png_dir, exist_ok=True)

    component_label = f"{field}_{{{component}}}^{{\\mathrm{{{range_type}}}}}"

    if result['i_screen'].nunique() == 1:
        subset = result.sort_values('omega')
        omega = subset['omega'].to_numpy()
        re_values = np.real(subset[col].to_numpy())
        im_values = np.imag(subset[col].to_numpy())
        abs_values = np.hypot(re_values, im_values)
        phase_values = np.arctan2(im_values, re_values)

        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        ax_re, ax_im, ax_abs, ax_phase = axes.flat
        ax_re.plot(omega, re_values, marker='.')
        ax_im.plot(omega, im_values, marker='.')
        ax_abs.plot(omega, abs_values, marker='.')
        ax_phase.plot(omega, phase_values, marker='.')

        ax_re.set_title(f"$\\mathrm{{Re}}({component_label})$")
        ax_im.set_title(f"$\\mathrm{{Im}}({component_label})$")
        ax_abs.set_title(f"$|{component_label}|$")
        ax_phase.set_title(f"$\\arg({component_label})$")
        ax_phase.set_ylim(-np.pi, np.pi)
        for ax in axes.flat:
            ax.set_xlabel("$\\omega / \\omega_1$ (units of the fundamental)")
            ax.grid(True, which='major')

        fig.suptitle("Spherical field component vs. $\\omega/\\omega_1$ (1 screen point)",
                     fontsize=13, fontweight='bold')
        plt.tight_layout()

        output_img = os.path.join(png_dir, f"spherical_point_spectrum_{range_type}_{field}{component}.png")
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)
        return

    N_theta = int(read_config_value('spherical_detector_N_theta', config_path)[0])
    N_phi = int(read_config_value('spherical_detector_N_phi', config_path)[0])
    grid_shape = (N_theta, N_phi)
    x_edges, y_edges, x_label, y_label, aspect = get_spherical_plot_grid(config_path)
    detector_geometry_label = get_detector_geometry_label('spherical', config_path)
    # Only the stereographic-projection grid (aspect == 'equal') is in length units; the
    # (theta, phi) angular fallback has no length scale for a w0-multiple axis to supplement.
    w0 = None
    if aspect == 'equal':
        axes_unit = read_config_value('spherical_detector_radius', config_path)[1]
        w0 = get_laser_lg_w0_in_axes_units(radiation_filepath, axes_unit)

    for i_omega, subset in result.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_value = subset['omega'].iloc[0]

        re_values = np.real(subset[col].to_numpy())
        im_values = np.imag(subset[col].to_numpy())
        abs_values = np.hypot(re_values, im_values)
        phase_values = np.arctan2(im_values, re_values)

        # Re/Im share the modulus panel's own [0, abs_max] scale, symmetrized to [-abs_max, abs_max],
        # rather than each auto-scaling to its own range -- see plot_radiation_field.py's
        # plot_radiation_component for the same convention and its rationale.
        abs_max = abs_values.max() if abs_values.size else None
        panels = [
            (re_values, f"$\\mathrm{{Re}}({component_label})$", 'RdBu_r', -abs_max, abs_max),
            (im_values, f"$\\mathrm{{Im}}({component_label})$", 'RdBu_r', -abs_max, abs_max),
            (abs_values, f"$|{component_label}|$", 'viridis', 0, abs_max),
            (phase_values, f"$\\arg({component_label})$", 'twilight', -np.pi, np.pi),
        ]

        fig, axes = plt.subplots(2, 2, figsize=(10, 8.5), layout='constrained')
        cbars = []
        for ax, (values, title, cmap, vmin, vmax) in zip(axes.flat, panels):
            sc = ax.pcolormesh(x_edges, y_edges, values.reshape(grid_shape), cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_aspect(aspect, adjustable='box')
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title(title)
            ax.grid(True, alpha=0.25)
            # shrink=0.85 (matching the Superradiant_Thomson reference plots' style) keeps the
            # colorbar from spanning the full axis height, which otherwise crowds the secondary
            # right y/w0 axis added below. No cbar.set_label: the panel title above already names
            # the quantity, and a second, redundant vertical label ate into the width available to
            # the right column's colorbar, misaligning it against the left column's.
            cbar = fig.colorbar(sc, ax=ax, shrink=0.85)
            if cmap == 'twilight':
                cbar.set_ticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
                cbar.ax.set_yticklabels(['$-\\pi$', '$-\\pi/2$', '$0$', '$\\pi/2$', '$\\pi$'])
            add_w0_secondary_axes(ax, w0)
            cbars.append(cbar)

        fig.suptitle(f"${component_label}$, $\\omega$ index {i_omega} ($\\omega/\\omega_1$={omega_value:.4g}), "
                     f"spherical detector, {detector_geometry_label}", fontsize=12)

        # constrained layout doesn't reliably column-align per-axes colorbars -- see
        # plot_radiation_field.py's plot_radiation_component for the full explanation (a "1e-N"
        # scientific offset label above the Re/Im colorbars vs. the twilight phase colorbar's plain
        # '$\pi$' tick label throws off the solver's per-column spacing). Force one layout pass,
        # freeze it, then snap each column's two colorbars to a shared x0.
        fig.canvas.draw()
        fig.set_layout_engine(None)
        for cbar in cbars:
            # fig.colorbar(ax=...) attaches an automatic locator to the colorbar axes that
            # recomputes its position relative to the parent ax on every draw -- including the one
            # savefig triggers -- which would silently undo the set_position() calls below unless
            # cleared first.
            cbar.ax.set_axes_locator(None)
        for top_cbar, bottom_cbar in zip(cbars[:2], cbars[2:]):
            top_pos = top_cbar.ax.get_position()
            bottom_pos = bottom_cbar.ax.get_position()
            bottom_cbar.ax.set_position([top_pos.x0, bottom_pos.y0, bottom_pos.width, bottom_pos.height])

        output_img = os.path.join(png_dir, f"spherical_field_{range_type}_{field}{component}_omega{i_omega}.png")
        # bbox_inches='tight' is safe here -- see plot_radiation_field.py's plot_radiation_component
        # for why: the layout is already frozen and every axes position fixed by hand, so 'tight'
        # only crops the outer margin and can no longer re-trigger a layout pass that un-aligns the
        # colorbars. It's needed again to keep the phase colorbar's own tick labels (pushed out to
        # the wider x0 set by the "1e-N" offset above Re/Im) from hanging past the figure's edge.
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: python3 {sys.argv[0]} <long|short|boundary> <E|B> <r|theta|phi> [path_to_radiation_field.dat]")
        sys.exit(1)

    range_type = sys.argv[1]
    if range_type not in ('long', 'short', 'boundary'):
        print("Error: first argument must be 'long', 'short', or 'boundary'")
        sys.exit(1)

    field = sys.argv[2]
    if field not in ('E', 'B'):
        print("Error: second argument must be 'E' or 'B'")
        sys.exit(1)

    component = sys.argv[3]
    if component not in ('r', 'theta', 'phi'):
        print("Error: third argument must be 'r', 'theta', or 'phi'")
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

    try:
        plot_spherical_field_component(range_type, field, component, input_file)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)
