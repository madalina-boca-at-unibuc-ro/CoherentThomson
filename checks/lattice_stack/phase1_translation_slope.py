#!/usr/bin/env python3
"""Sheet-translation test: translate the flat lattice along z and fit the complex field ratio to exp(-i s z).

Needs no stack code, only beam_center_z. Compares the fitted slope s with the analytic q_N = N k_1 + k_L.
Run: .venv/bin/python theory/lattice_stack_tests/phase1_translation_slope.py
"""
import numpy as np

from lattice_stack_common import fit_scale, laser_kL_units, read_field, run_variant

OFFSETS = [0.0, 1e-5, -1e-5, 1e-4, -1e-4, 3e-4, -3e-4, 7e-4, -7e-4, 1.4e-3, -1.4e-3, 1.4e-2, -1.4e-2]  # in lambda_L
fields = {}
for dz in OFFSETS:
    run = run_variant("phase1_slope", f"dz_{dz:+.1e}", {"beam_center_z": f"{dz:.10e} lambda"})
    if dz == 0.0:
        k1, q1 = laser_kL_units(run)
    fields[dz] = read_field(run)[1]
print(f"analytic: k_1 = {k1:.6f} k_L, q_1 = {q1:.6f} k_L\n")

for N, i_omega in ((1, 0), (2, 1)):
    q_analytic = N * k1 + 1.0
    print(f"N = {N}: analytic q_N = {q_analytic:.6f} k_L")
    print("   dz [lambda]    phase of c [rad]  predicted [rad]   diff [rad]     |c|-1       eta")
    xs, ys = [], []
    for dz in sorted(OFFSETS):
        c, eta = fit_scale(fields[0.0][i_omega], fields[dz][i_omega])
        predicted = -q_analytic * 2 * np.pi * dz
        # Unwrap against the prediction: the fit only needs the difference from it.
        diff = np.angle(c * np.exp(-1j * predicted))
        print(f"  {dz:+.3e}    {predicted + diff:+14.6f}  {predicted:+14.6f}  {diff:+.3e}  {abs(c) - 1:+.2e}  {eta:.2e}")
        xs.append(2 * np.pi * dz)
        ys.append(diff)
    slope_diff = np.polyfit(xs, ys, 1)[0]
    print(f"  fitted slope s = {q_analytic - slope_diff:.6f} k_L, s - q_N = {-slope_diff:+.3e} k_L, "
          f"relative {-slope_diff / q_analytic:+.2e}\n")
