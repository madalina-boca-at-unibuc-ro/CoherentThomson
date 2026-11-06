#!/usr/bin/env python3
"""b with the exact q against the field: at H = Lambda_1 both are ~7e-6, and b's phase is the conjugate of the field's.

Run from the repository root: .venv/bin/python theory/lattice_stack_tests/form_factor_exact_q.py
"""
import numpy as np

from lattice_stack_common import fit_scale, form_factor_b, laser_kL_units, read_field, run_variant

POINT = {"circular_detector_N_R": 1, "circular_detector_N_phi": 1, "N_harmonics": 1, "beam_lattice_spacing": "10 lambda"}

flat = run_variant("form_factor_exact_q", "flat", POINT)
k1, q1 = laser_kL_units(flat)
d = 1 / q1 / 14
run = run_variant("form_factor_exact_q", "M14", dict(POINT, beam_stack_sheets=14, beam_stack_spacing=f"{d:.12e} lambda"))
c, _ = fit_scale(read_field(flat)[1][0], read_field(run)[1][0])
print(f"head-on, H = 2 pi / |k_s - k_L|: b = {form_factor_b(run):.3e}, field c = {c:.3e}")

SMALL = dict(POINT, laser_type="plane_wave", beam_cylinder_radius="0.1 lambda", beam_lattice_spacing="0.04 lambda",
             laser_direction_theta="0.9 pi", laser_a0=0.2)
tilted_flat = run_variant("form_factor_exact_q", "tilted_flat", SMALL)
for key in ("beam_center_x", "beam_center_z"):
    run = run_variant("form_factor_exact_q", f"tilted_{key}", dict(SMALL, **{key: "1e-3 lambda"}))
    c, _ = fit_scale(read_field(tilted_flat)[1][0], read_field(run)[1][0])
    ratio = form_factor_b(run) / form_factor_b(tilted_flat)
    print(f"tilted, {key}: arg(b ratio) = {np.angle(ratio):+.6f}, -arg(c) = {-np.angle(c):+.6f}, "
          f"difference {np.angle(ratio) + np.angle(c):+.1e} rad")
