"""
Computes and plots the screen-level angular-momentum volume density along Oz, split into a direct
(Poynting-moment) representation and a spin/orbital decomposition, per the reorganized
theory/angular-momentum-density-and-flux.md -- a from-scratch re-derivation, NOT a port of
plot_angular_momentum_density_screen.py's theory/angular_momentum_density_screen_from_Faraday.md.

This is a deliberately separate, from-scratch script (not a replacement) while
theory/angular-momentum-density-and-flux.md is still being reorganized -- plot_angular_momentum_flux.py,
plot_angular_momentum_flux_screen.py, and plot_angular_momentum_density_screen.py are all kept
unchanged and continue to implement their own theory docs exactly as before, so results here can be
cross-checked against them rather than silently superseding them. Once the new doc settles, this
script (or its successor) is expected to replace the old density script -- see CLAUDE.md.

Reuses plot_angular_momentum_flux.py's shared machinery (rotation, local screen coordinates, E/B
extraction, constants) and plot_angular_momentum_flux_screen.py's own (raw-omega reconstruction,
area-quadrature weights) rather than duplicating any of it -- the geometry/rotation machinery is
identical across every one of these scripts; only the final formula combining E/B into an
angular-momentum quantity differs. The azimuthal derivative, however, is NOT reused from
plot_angular_momentum_density_screen.py's FFT-based _azimuthal_derivative_rectangular: this doc's own
Section 4 specifies plain central finite differences on the Cartesian grid (not a spectral method),
so this script implements that literally instead (_azimuthal_derivative_rectangular below) --
_azimuthal_derivative_circular IS reused unchanged, since a periodic centered difference in phi at
fixed r is exactly what the doc's own polar-grid recipe asks for too.

Component conventions (theory doc Section 2.3 -- identical to plot_angular_momentum_density_screen.py's
own, and to extract_rotated_faraday_fields's E/B mapping): E_i = c*F^{i0}, B_1=-F^{23}, B_2=F^{13},
B_3=-F^{12}. extract_rotated_faraday_fields's returned Ex/Ey/Ez/Bx/By/Bz are already the genuine,
correctly-scaled physical field components matching this convention directly -- no extra factor of c
needs to be (re)introduced when reading E/B out of that function; c only appears where the doc's own
formulas have an explicit c or c^2 factor.

Four quantities are computed at every screen point and every frequency (theory doc Sections 3-4),
using the doc's own Parseval one-sided-spectrum prefactor 4*pi (Section 2.2):
- j3_direct  = 4*pi*epsilon_0 * Re[x*(E3 B1* - E1 B3*) - y*(E2 B3* - E3 B2*)]
- j3_spin    = s_z = (pi*epsilon_0/omega) * Im[(E1* E2 - E2* E1) + c^2*(B1* B2 - B2* B1)]
- j3_orbital = (pi*epsilon_0/omega) * Im[E1* dphi(E1) + E2* dphi(E2) + E3* dphi(E3)
                                          + c^2*(B1* dphi(B1) + B2* dphi(B2) + B3* dphi(B3))]
- j3_sum     = j3_spin + j3_orbital

dphi = x*d/dy - y*d/dx (Section 4).

**The exact pointwise identity** (theory doc Section 4) is
    j3_direct = j3_spin + j3_orbital + dD_z/dz + div_perp(M_perp),
where dD_z/dz = (1/2)*d/dz(x*s_x + y*s_y) is a genuine z-derivative -- needing the field at
neighboring z-planes, which no single-z export (radiation_field.dat/incident_field.dat included) can
provide -- and div_perp(M_perp)'s own SCREEN INTEGRAL is a boundary contour term the theory doc's
Section 5 shows vanishes exactly for a large-enough screen (fields ~0 at the edge); per the user's own
convention this project always sizes screens that way, so it is dropped from the integrated identity
below rather than computed pointwise. That leaves the practically checkable, purely INTEGRATED form
(Section 5, item 2):
    J3_direct(omega) ~= J3_spin(omega) + J3_orbital(omega) + dJ3(omega),
    dJ3(omega) = d/dz[ Integral (1/2)*(x*s_x + y*s_y) dx dy ]_{z=0}
               ~= [S(z=+delta_z) - S(z=-delta_z)] / (2*delta_z),
implemented by verify_direct_spin_orbital_identity below, using
Core::Radiation::export_incident_field_fourier's two extra z-shifted exports
(incident_field_zminus.dat/incident_field_zplus.dat, config key "incident_field_delta_z") -- the ONLY
field this project can evaluate off a single fixed z-plane without a full re-run, so this identity
check only ever runs in --incident mode.

Scope and caveats (same as every sibling script here):
- Only 'rectangular'/'circular' detectors (a flat screen with one well-defined local (x, y) plane).
- Same un-verified FT-normalization caveat as every other angular-momentum script: the theory doc's
  boxed formulas are applied literally, without independently re-checking that radiation_field.dat's
  actual export normalization matches the doc's assumed Fourier convention.
- Accepts --incident (Core::Radiation::export_incident_field_fourier's known-analytic incident-beam
  reference field) exactly like its sibling scripts, with the same _UNCALIBRATED_MAGNITUDE_CAVEAT.
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
    get_rectangular_cell_edges,
    get_circular_cell_edges,
    get_detector_geometry_label,
)
from plot_angular_momentum_flux import (
    C_LIGHT,
    EPSILON_0,
    convert_unit_to_number,
    get_canonical_to_local_rotation,
    get_screen_coordinates_local_au,
    extract_rotated_faraday_fields,
)
from plot_angular_momentum_flux_screen import (
    get_fundamental_frequency_au,
    get_screen_area_weights_au,
)
from plot_angular_momentum_density_screen import _azimuthal_derivative_circular
from utils.w0_axes_utils import get_laser_lg_w0_in_axes_units, add_w0_secondary_axes

MODULE_NAME = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

QUANTITIES = (
    ('j3_spin', 'spin', r'$dj_3^{\mathrm{spin}}/d\omega$ (a.u.)'),
    ('j3_orbital', 'orbital', r'$dj_3^{\mathrm{orbital}}/d\omega$ (a.u.)'),
    ('j3_sum', 'sum', r'$dj_3^{\mathrm{sum}}/d\omega$ (a.u.)'),
    ('j3_direct', 'direct', r'$dj_3^{\mathrm{direct}}/d\omega$ (a.u.)'),
)


def _azimuthal_derivative_rectangular(values, x1_grid, x2_grid):
    """
    dphi_F = x*dF/dy - y*dF/dx (theory doc Section 4), via plain central finite differences
    (np.gradient) in x and y -- the doc's own specified numerical recipe for a Cartesian grid
    (explicit centered-difference stencils with Delta_x/Delta_y denominators), deliberately NOT the
    FFT-based spectral derivative plot_angular_momentum_density_screen.py uses for the analogous
    quantity under the older theory doc (see this module's own docstring for why). `values`/
    `x1_grid`/`x2_grid` all shape (Nx, Ny), meshgrid('ij') convention. A size-1 axis has no defined
    derivative (a single-row/column detector has no spatial variation to differentiate along that
    axis) -- np.gradient requires at least 2 points along an axis, so that case is treated as zero.
    """
    dF_dx = (np.gradient(values, x1_grid[:, 0], axis=0) if values.shape[0] > 1
             else np.zeros_like(values))
    dF_dy = (np.gradient(values, x2_grid[0, :], axis=1) if values.shape[1] > 1
             else np.zeros_like(values))
    return x1_grid * dF_dy - x2_grid * dF_dx


def _j3_direct(x, y, E1, E2, E3, B1, B2, B3):
    """
    dj3_direct/domega = 4*pi*epsilon_0 * Re[x*(E3 B1* - E1 B3*) - y*(E2 B3* - E3 B2*)]
    (theory doc Section 3.1), vectorized over any matching-shape complex field-component arrays and
    real local-frame position arrays. Note the 4*pi*epsilon_0 prefactor (Section 2.2's Parseval
    one-sided factor) -- NOT plot_angular_momentum_flux.py's own angular_momentum_flux_density, which
    shares the identical bracket but carries an epsilon_0/pi prefactor instead, from a different
    (older, since-revised) theory-doc derivation; the two are proportional but not interchangeable.
    """
    return 4.0 * np.pi * EPSILON_0 * np.real(
        x * (E3 * np.conj(B1) - E1 * np.conj(B3)) - y * (E2 * np.conj(B3) - E3 * np.conj(B2))
    )


def _spin_vector_components(E1, E2, E3, B1, B2, B3, omega_raw):
    """
    Returns (sx, sy, sz): the full physical spin-density vector s (theory doc Section 3.2's
    s = (4*pi*epsilon_0/(4*omega)) * Im(E* x E + c^2 B* x B), i.e. pi*epsilon_0/omega, componentwise
    via cyclic (1,2,3) -> (2,3,1) -> (3,1,2) index permutation). sz is j3_spin (Section 3.2); sx/sy
    (Section 4.1) are only needed for the z-derivative divergence-correction term
    verify_direct_spin_orbital_identity computes -- not otherwise part of the four plotted quantities.
    """
    prefactor = np.pi * EPSILON_0 / omega_raw
    sx = prefactor * np.imag((np.conj(E2) * E3 - np.conj(E3) * E2)
                             + C_LIGHT**2 * (np.conj(B2) * B3 - np.conj(B3) * B2))
    sy = prefactor * np.imag((np.conj(E3) * E1 - np.conj(E1) * E3)
                             + C_LIGHT**2 * (np.conj(B3) * B1 - np.conj(B1) * B3))
    sz = prefactor * np.imag((np.conj(E1) * E2 - np.conj(E2) * E1)
                             + C_LIGHT**2 * (np.conj(B1) * B2 - np.conj(B2) * B1))
    return sx, sy, sz


def compute_angular_momentum_density_and_flux(radiation_filepath):
    """
    Returns (result, detector_type, config_path). `result` has one row per radiation_field.dat row
    (i_omega, omega, i_screen), the local-frame (x_local, y_local) position, and the four
    angular-momentum volume densities (j3_spin, j3_orbital, j3_sum, j3_direct) from
    theory/angular-momentum-density-and-flux.md, evaluated on the physical total field (LR+SR+BR --
    BR is identically zero for radiation_formula="direct", so this is correct regardless of which
    formula produced the file).
    """
    config_path = os.path.join(os.path.dirname(radiation_filepath), 'config.cfg')
    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"'{config_path}' not found -- this run predates per-run config snapshots "
            "(Core::IoUtils::copy_config_to_run_directory); re-run the solver to regenerate it.")

    run_dir = os.path.dirname(radiation_filepath)
    omega_1 = get_fundamental_frequency_au(run_dir)

    detector_type = read_config_value('detector_type', config_path)[0]
    if detector_type not in ('rectangular', 'circular'):
        raise ValueError(
            f"detector_type '{detector_type}' is not a flat screen with a single well-defined local "
            "(x1, x2) plane -- theory/angular-momentum-density-and-flux.md assumes one flat "
            "transverse screen, which only 'rectangular' and 'circular' detectors are.")
    x_local, y_local = get_screen_coordinates_local_au(detector_type, config_path)

    if detector_type == 'rectangular':
        grid_shape = (int(read_config_value('rectangular_detector_Nx', config_path)[0]),
                      int(read_config_value('rectangular_detector_Ny', config_path)[0]))
    else:
        grid_shape = (int(read_config_value('circular_detector_N_R', config_path)[0]),
                      int(read_config_value('circular_detector_N_phi', config_path)[0]))
    x1_grid = x_local.reshape(grid_shape)
    x2_grid = y_local.reshape(grid_shape)

    data = pd.read_csv(radiation_filepath, sep=' ', comment='#')
    R_to_local = get_canonical_to_local_rotation(config_path)

    n_skipped = 0
    rows = []
    for i_omega, subset in data.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_ratio = subset['omega'].iloc[0]
        omega_raw = omega_ratio * omega_1
        i_screen = subset['i_screen'].to_numpy()
        x = x_local[i_screen]
        y = y_local[i_screen]

        f = extract_rotated_faraday_fields(subset, R_to_local)
        E1 = f['Ex_l'] + f['Ex_s'] + f['Ex_b']
        E2 = f['Ey_l'] + f['Ey_s'] + f['Ey_b']
        E3 = f['Ez_l'] + f['Ez_s'] + f['Ez_b']
        B1 = f['Bx_l'] + f['Bx_s'] + f['Bx_b']
        B2 = f['By_l'] + f['By_s'] + f['By_b']
        B3 = f['Bz_l'] + f['Bz_s'] + f['Bz_b']

        j3_direct = _j3_direct(x, y, E1, E2, E3, B1, B2, B3)

        if omega_raw > 0:
            _, _, j3_spin = _spin_vector_components(E1, E2, E3, B1, B2, B3, omega_raw)

            derivative = _azimuthal_derivative_rectangular if detector_type == 'rectangular' \
                else _azimuthal_derivative_circular
            dphi = {}
            for name, values in (('E1', E1), ('E2', E2), ('E3', E3), ('B1', B1), ('B2', B2), ('B3', B3)):
                grid = values.reshape(grid_shape)
                dphi_grid = (derivative(grid, x1_grid, x2_grid) if detector_type == 'rectangular'
                            else derivative(grid))
                dphi[name] = dphi_grid.ravel()

            orbital_complex = (np.conj(E1) * dphi['E1'] + np.conj(E2) * dphi['E2'] + np.conj(E3) * dphi['E3']
                               + C_LIGHT**2 * (np.conj(B1) * dphi['B1'] + np.conj(B2) * dphi['B2']
                                              + np.conj(B3) * dphi['B3']))
            j3_orbital = (np.pi * EPSILON_0 / omega_raw) * np.imag(orbital_complex)
        else:
            n_skipped += 1
            j3_spin = np.full_like(j3_direct, np.nan)
            j3_orbital = np.full_like(j3_direct, np.nan)

        j3_sum = j3_spin + j3_orbital

        rows.append(pd.DataFrame({
            'i_omega': i_omega, 'omega': omega_ratio, 'i_screen': i_screen,
            'x_local': x, 'y_local': y,
            'j3_spin': j3_spin, 'j3_orbital': j3_orbital,
            'j3_sum': j3_sum, 'j3_direct': j3_direct,
        }))

    if n_skipped:
        print(f"Warning: {n_skipped} frequency/frequencies have omega<=0; spin/orbital densities left "
              "as NaN there (undefined at omega=0).")

    result = pd.concat(rows, ignore_index=True)
    return result, detector_type, config_path


def integrate_angular_momentum_density_and_flux(result, detector_type, config_path):
    """
    Screen integration (theory doc Section 5, step 4): J3^q(omega) = sum_p Delta A_p * j3^q_p, for
    q in {spin, orbital, sum, direct} -- reuses plot_angular_momentum_flux_screen.py's
    get_screen_area_weights_au. Returns a DataFrame with one row per unique i_omega.
    """
    area_weights = get_screen_area_weights_au(detector_type, config_path)

    rows = []
    for i_omega, subset in result.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        if len(subset) != len(area_weights):
            raise ValueError(
                f"i_omega={i_omega}: {len(subset)} screen points in radiation_field.dat but "
                f"{len(area_weights)} expected from the detector's own configured grid -- geometry "
                "mismatch (e.g. the config was edited after this run).")
        rows.append({
            'i_omega': i_omega,
            'omega': subset['omega'].iloc[0],
            'j3_spin': np.sum(area_weights * subset['j3_spin'].to_numpy()),
            'j3_orbital': np.sum(area_weights * subset['j3_orbital'].to_numpy()),
            'j3_sum': np.sum(area_weights * subset['j3_sum'].to_numpy()),
            'j3_direct': np.sum(area_weights * subset['j3_direct'].to_numpy()),
        })
    return pd.DataFrame(rows).sort_values('i_omega').reset_index(drop=True)


def _transverse_spin_moment_integral(field_filepath, config_path, detector_type, omega_1):
    """
    Reads one z-shifted incident-field file (incident_field_zminus.dat or incident_field_zplus.dat)
    and returns a DataFrame with one row per i_omega: i_omega, omega, and S = sum_p Delta_A_p *
    (1/2)*(x_p*sx_p + y_p*sy_p) -- the screen-integrated transverse spin moment
    theory/angular-momentum-density-and-flux.md's Section 4.1/5 divergence-correction term needs a
    central-difference d/dz of. Only sx/sy (not sz/j3_spin, not j3_orbital/j3_direct) are needed here.
    """
    x_local, y_local = get_screen_coordinates_local_au(detector_type, config_path)
    data = pd.read_csv(field_filepath, sep=' ', comment='#')
    R_to_local = get_canonical_to_local_rotation(config_path)
    area_weights = get_screen_area_weights_au(detector_type, config_path)

    rows = []
    for i_omega, subset in data.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_ratio = subset['omega'].iloc[0]
        omega_raw = omega_ratio * omega_1
        i_screen = subset['i_screen'].to_numpy()
        x = x_local[i_screen]
        y = y_local[i_screen]

        f = extract_rotated_faraday_fields(subset, R_to_local)
        E1 = f['Ex_l'] + f['Ex_s'] + f['Ex_b']
        E2 = f['Ey_l'] + f['Ey_s'] + f['Ey_b']
        E3 = f['Ez_l'] + f['Ez_s'] + f['Ez_b']
        B1 = f['Bx_l'] + f['Bx_s'] + f['Bx_b']
        B2 = f['By_l'] + f['By_s'] + f['By_b']
        B3 = f['Bz_l'] + f['Bz_s'] + f['Bz_b']

        if omega_raw > 0:
            sx, sy, _ = _spin_vector_components(E1, E2, E3, B1, B2, B3, omega_raw)
            S = np.sum(area_weights * 0.5 * (x * sx + y * sy))
        else:
            S = np.nan

        rows.append({'i_omega': i_omega, 'omega': omega_ratio, 'S': S})
    return pd.DataFrame(rows).sort_values('i_omega').reset_index(drop=True)


def verify_direct_spin_orbital_identity(run_dir, integrated):
    """
    Verifies the theory doc's integrated exact identity (Section 5, item 2, boundary contour term
    dropped per the module docstring's own reasoning):
        J3_direct(omega) ~= J3_spin(omega) + J3_orbital(omega) + dJ3(omega),
        dJ3(omega) = [S(z=+delta_z) - S(z=-delta_z)] / (2*delta_z).

    Only possible for the incident beam -- requires incident_field_zminus.dat/incident_field_zplus.dat
    (Core::Radiation::export_incident_field_fourier's two z-shifted exports) next to incident_field.dat;
    returns None (with a printed explanation, not an error) if either is missing, e.g. a run predating
    this feature, or `incident_field_delta_z` absent from that run's own config.cfg snapshot.

    `integrated` is compute_angular_momentum_density_and_flux/integrate_angular_momentum_density_and_flux's
    own result for incident_field.dat itself (the z=0 plane) -- reused rather than recomputed.
    """
    zminus_path = os.path.join(run_dir, "incident_field_zminus.dat")
    zplus_path = os.path.join(run_dir, "incident_field_zplus.dat")
    if not (os.path.exists(zminus_path) and os.path.exists(zplus_path)):
        print("Skipping the direct = spin+orbital+dJ3/dz identity check: incident_field_zminus.dat/"
              "incident_field_zplus.dat not found next to incident_field.dat (this run predates the "
              "z-offset incident-field export -- re-run the solver to regenerate them).")
        return None

    config_path = os.path.join(run_dir, 'config.cfg')
    detector_type = read_config_value('detector_type', config_path)[0]
    omega_1 = get_fundamental_frequency_au(run_dir)
    delta_z_val, delta_z_unit = read_config_value('incident_field_delta_z', config_path)
    delta_z_au = float(delta_z_val) * convert_unit_to_number(delta_z_unit, config_path)

    S_minus = _transverse_spin_moment_integral(zminus_path, config_path, detector_type, omega_1)
    S_plus = _transverse_spin_moment_integral(zplus_path, config_path, detector_type, omega_1)

    merged = S_minus.merge(S_plus, on=['i_omega', 'omega'], suffixes=('_minus', '_plus'))
    merged['dJ3'] = (merged['S_plus'] - merged['S_minus']) / (2.0 * delta_z_au)

    result = integrated.merge(merged[['i_omega', 'dJ3']], on='i_omega')
    result['spin_plus_orbital_plus_dJ3'] = result['j3_spin'] + result['j3_orbital'] + result['dJ3']
    result['residual'] = result['j3_direct'] - result['spin_plus_orbital_plus_dJ3']
    with np.errstate(divide='ignore', invalid='ignore'):
        result['relative_residual'] = result['residual'] / result['j3_direct']
    return result


_RUN_LOG_SECTION_HEADER = "Angular-momentum density-and-flux integration (plot_angular_momentum_density_and_flux.py)"

# Same caveat as every other angular-momentum script's own _UNCALIBRATED_MAGNITUDE_CAVEAT -- see
# plot_angular_momentum_density_screen.py's copy for the full reasoning.
_UNCALIBRATED_MAGNITUDE_CAVEAT = (
    "NOTE: incident_field.dat's absolute field magnitude has no Fourier-transform normalization "
    "applied (see this script's _UNCALIBRATED_MAGNITUDE_CAVEAT) -- only RATIOS between the columns "
    "below (e.g. j3_orbital/j3_spin) are physically meaningful; the absolute values are on an "
    "arbitrary scale and are NOT comparable to a real scattered-radiation run's own j3_* values."
)


def append_integration_to_run_log(run_dir, detector_type, integrated, label_suffix=""):
    """
    Appends (or, on a re-run, replaces) a j3 summary table in the run's own run_log.txt -- same
    pattern as plot_angular_momentum_density_screen.py's own append_integration_to_run_log, with its
    own distinct section header so re-running either script never clobbers the other's section. Also
    prints (and logs) the sum-vs-direct integrated comparison the theory doc invites (Section 5, step
    4: "Verify convergence: J3_direct(omega) ~= J3_sum(omega)").
    """
    header = _RUN_LOG_SECTION_HEADER + label_suffix
    log_path = os.path.join(run_dir, "run_log.txt")
    if not os.path.exists(log_path):
        raise FileNotFoundError(
            f"'{log_path}' not found -- this run predates run_log.txt (Logging::write_run_log); "
            "re-run the solver to regenerate it.")

    with open(log_path) as f:
        content = f.read()

    marker = f"\n{header}\n"
    marker_index = content.find(marker)
    if marker_index != -1:
        content = content[:marker_index]

    lines = [content.rstrip("\n"), "", header, "-" * len(header)]
    if label_suffix:
        lines += [f"  {_UNCALIBRATED_MAGNITUDE_CAVEAT}", ""]
    lines += [
             "  Rough sum_p Delta_A_p * j3_p screen integration (theory/angular-momentum-density-and-flux.md",
             f"  Section 5), {detector_type} detector, trapezoidal quadrature in "
             f"{'(x, y)' if detector_type == 'rectangular' else '(r^2, phi)'}. j3_sum and j3_direct are"
             "  NOT expected to agree pointwise, but per the theory doc their screen integrals should"
             "  converge to the same value -- see the 'sum/direct' column below.",
             f"  {'i_omega':>8} {'omega/omega_1':>14} {'J3_spin [a.u.]':>16} {'J3_orbital [a.u.]':>18} "
             f"{'J3_sum [a.u.]':>16} {'J3_direct [a.u.]':>18} {'sum/direct':>12}"]
    for _, row in integrated.iterrows():
        ratio = row['j3_sum'] / row['j3_direct'] if row['j3_direct'] != 0 else float('nan')
        lines.append(f"  {int(row['i_omega']):>8d} {row['omega']:>14.6f} {row['j3_spin']:>16.6e} "
                     f"{row['j3_orbital']:>18.6e} {row['j3_sum']:>16.6e} "
                     f"{row['j3_direct']:>18.6e} {ratio:>12.4f}")

    with open(log_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Appended angular-momentum density-and-flux integration to {log_path}")


_IDENTITY_RUN_LOG_SECTION_HEADER = ("Angular-momentum exact-identity check (plot_angular_momentum_"
                                    "density_and_flux.py, incident beam only)")


def append_identity_check_to_run_log(run_dir, identity_result):
    """
    Appends (or, on a re-run, replaces) the direct = spin+orbital+dJ3 identity table
    (verify_direct_spin_orbital_identity's result) in the run's own run_log.txt, under its own
    section header distinct from append_integration_to_run_log's -- only ever called in --incident
    mode (see plot_angular_momentum_density_and_flux), so this section only ever appears once per run,
    never duplicated/overwritten by a non-incident invocation.
    """
    header = _IDENTITY_RUN_LOG_SECTION_HEADER
    log_path = os.path.join(run_dir, "run_log.txt")
    if not os.path.exists(log_path):
        raise FileNotFoundError(
            f"'{log_path}' not found -- this run predates run_log.txt (Logging::write_run_log); "
            "re-run the solver to regenerate it.")

    with open(log_path) as f:
        content = f.read()

    marker = f"\n{header}\n"
    marker_index = content.find(marker)
    if marker_index != -1:
        content = content[:marker_index]

    lines = [content.rstrip("\n"), "", header, "-" * len(header),
            f"  {_UNCALIBRATED_MAGNITUDE_CAVEAT}", "",
             "  theory/angular-momentum-density-and-flux.md Section 4/5's exact pointwise identity,",
             "  screen-integrated with the boundary contour term dropped (assumed zero -- see this",
             "  script's module docstring):",
             "    J3_direct ~= J3_spin + J3_orbital + dJ3,  dJ3 = [S(z+dz)-S(z-dz)]/(2*dz)",
             f"  {'i_omega':>8} {'omega/omega_1':>14} {'J3_direct':>14} {'spin+orb+dJ3':>14} "
             f"{'residual':>14} {'rel. resid.':>12}"]
    for _, row in identity_result.iterrows():
        lines.append(f"  {int(row['i_omega']):>8d} {row['omega']:>14.6f} {row['j3_direct']:>14.6e} "
                     f"{row['spin_plus_orbital_plus_dJ3']:>14.6e} {row['residual']:>14.6e} "
                     f"{row['relative_residual']:>12.4f}")

    with open(log_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Appended angular-momentum exact-identity check to {log_path}")


def _plot_one_quantity(result, detector_type, config_path, radiation_filepath, png_dir, column, name, axis_label,
                       label_suffix="", title_suffix=""):
    """
    Renders `column`: a line plot vs. omega/omega_1 for a 1x1 detector (dense_frequency_spectrum
    workflow), or one heatmap PNG per frequency on the detector's own native grid otherwise -- same
    rendering convention as every sibling angular-momentum script's own _plot_one_quantity.
    """
    if result['i_screen'].nunique() == 1:
        subset = result.sort_values('omega')
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(subset['omega'], subset[column], marker='.')
        ax.set_xlabel("$\\omega / \\omega_1$ (units of the fundamental)")
        ax.set_ylabel(axis_label)
        ax.grid(True)
        fig.suptitle(f"Spectral angular-momentum density ({name}){title_suffix}", fontsize=13, fontweight='bold')
        plt.tight_layout()

        output_img = os.path.join(png_dir, f"angular_momentum_density_and_flux_{name}{label_suffix}_spectrum.png")
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)
        return

    if detector_type == 'rectangular':
        Nx = int(read_config_value('rectangular_detector_Nx', config_path)[0])
        Ny = int(read_config_value('rectangular_detector_Ny', config_path)[0])
        grid_shape = (Nx, Ny)
        x_edges, y_edges = get_rectangular_cell_edges(config_path)
        axes_unit = read_config_value('rectangular_detector_x_min', config_path)[1]
    else:
        N_R = int(read_config_value('circular_detector_N_R', config_path)[0])
        N_phi = int(read_config_value('circular_detector_N_phi', config_path)[0])
        grid_shape = (N_R, N_phi)
        x_edges, y_edges = get_circular_cell_edges(config_path)
        axes_unit = read_config_value('circular_detector_R_min', config_path)[1]

    x_label = "$x$" + (f" [{axes_unit}]" if axes_unit else "")
    y_label = "$y$" + (f" [{axes_unit}]" if axes_unit else "")
    detector_geometry_label = get_detector_geometry_label(detector_type, config_path)
    w0 = get_laser_lg_w0_in_axes_units(radiation_filepath, axes_unit)

    for i_omega, subset in result.groupby('i_omega'):
        subset = subset.sort_values('i_screen')
        omega_value = subset['omega'].iloc[0]

        values = subset[column].to_numpy().reshape(grid_shape)
        finite = values[np.isfinite(values)]
        vmax = np.max(np.abs(finite)) if finite.size and np.any(finite) else 1.0

        fig, ax = plt.subplots(figsize=(7, 6))
        sc = ax.pcolormesh(x_edges, y_edges, values, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        fig.colorbar(sc, ax=ax, label=axis_label, pad=0.15)
        add_w0_secondary_axes(ax, w0)

        fig.suptitle(f"Angular-momentum density ({name}){title_suffix}, $\\omega$ index {i_omega} "
                     f"($\\omega/\\omega_1$={omega_value:.4g}), {detector_type} detector, "
                     f"{detector_geometry_label}", fontsize=11, fontweight='bold')
        plt.tight_layout()

        output_img = os.path.join(png_dir, f"angular_momentum_density_and_flux_{name}{label_suffix}_omega{i_omega}.png")
        plt.savefig(output_img, dpi=200, bbox_inches='tight')
        print(f"Successfully saved plot to {output_img}")
        plt.close(fig)


def plot_angular_momentum_density_and_flux(radiation_filepath):
    """
    Computes and renders the four screen-level angular-momentum volume densities (spin, orbital, sum,
    direct) from theory/angular-momentum-density-and-flux.md. Also integrates each over the screen
    and appends the resulting J3(omega) table to the run's own run_log.txt.

    `radiation_filepath` may point at either a run's radiation_field.dat or its incident_field.dat
    sibling (Core::Radiation::export_incident_field_fourier), exactly like every sibling
    angular-momentum script -- detected by filename so incident-beam runs get their own PNG names and
    run_log section instead of overwriting the real run's.
    """
    is_incident = os.path.basename(radiation_filepath) == "incident_field.dat"
    label_suffix = "_incident" if is_incident else ""
    title_suffix = " (incident beam)" if is_incident else ""
    run_log_header_suffix = " -- incident beam" if is_incident else ""

    result, detector_type, config_path = compute_angular_momentum_density_and_flux(radiation_filepath)

    png_dir = os.path.join(os.path.dirname(radiation_filepath), "png_folder", MODULE_NAME)
    os.makedirs(png_dir, exist_ok=True)

    for column, name, axis_label in QUANTITIES:
        _plot_one_quantity(result, detector_type, config_path, radiation_filepath, png_dir, column, name, axis_label,
                           label_suffix=label_suffix, title_suffix=title_suffix)

    integrated = integrate_angular_momentum_density_and_flux(result, detector_type, config_path)
    if is_incident:
        print(_UNCALIBRATED_MAGNITUDE_CAVEAT)
    print(integrated.to_string(index=False))
    run_dir = os.path.dirname(radiation_filepath)
    append_integration_to_run_log(run_dir, detector_type, integrated, label_suffix=run_log_header_suffix)

    # The exact direct=spin+orbital+dJ3 identity needs the incident beam's two z-shifted siblings
    # (incident_field_zminus.dat/incident_field_zplus.dat) -- only ever available in --incident mode,
    # since the real scattered field has no equivalent off-plane evaluation. See
    # verify_direct_spin_orbital_identity's own doc comment.
    if is_incident:
        identity_result = verify_direct_spin_orbital_identity(run_dir, integrated)
        if identity_result is not None:
            print(identity_result.to_string(index=False))
            append_identity_check_to_run_log(run_dir, identity_result)


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

    try:
        plot_angular_momentum_density_and_flux(input_file)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)
