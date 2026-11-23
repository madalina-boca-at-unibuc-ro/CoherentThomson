#!/usr/bin/env python3
"""Phase 1 checkpoint: uniform stacks at fixed total charge against the discrete prediction B_M(q d).

Laser, transverse lattice and screen are those of run A (a0 = 0.1) on a 16 x 8 screen, harmonics N = 1, 2. c is the
best complex factor in F_M = c F_1 over the screen and all components; `resid` = |F_M - B F_1| / |F_1|, absolute so
that it stays meaningful at the minima of B. Takes about 30 min on 24 threads; finished runs are reused.

Run from the repository root: .venv/bin/python theory/lattice_stack_tests/phase1_stack_scan.py
"""
import numpy as np

from lattice_stack_common import fit_scale, form_factor_b, laser_kL_units, read_field, run_variant, stack_factor

TAU = 2 * np.pi


def stack(group, name, sheets, spacing, extra=None):
    # Since the 2026-11-09 rework a stack is given by its width H = M d (beam_cylinder_height).
    overrides = {"beam_stack_sheets": sheets, "beam_cylinder_height": f"{sheets * spacing:.12e} lambda"}
    overrides.update(extra or {})
    return run_variant(group, name, overrides)


def norm(a):
    return np.sqrt(np.vdot(a, a).real)


reference = run_variant("phase1_scan", "flat", {})
k1, q1 = laser_kL_units(reference)
F1 = read_field(reference)[1]
# Slopes: q_N = N q_1 with q_1 analytic, and the translation prediction 2 N k_1 / (1 + beta0) that the slope test found.
pz = 13.483
beta0 = pz / np.hypot(1.0, pz)
s1 = 2 * k1 / (1 + beta0)
print(f"k_1 = {k1:.6f} k_L, analytic q_1 = {q1:.6f} k_L, translation slope s_1 = {s1:.6f} k_L")
print(f"Lambda_1 = 2 pi / q_1 = {1 / q1:.6e} lambda_L\n")


def report(label, run, sheets, spacing):
    F = read_field(run)[1]
    b = form_factor_b(run)
    b_pred = stack_factor(TAU * q1 * spacing, sheets)
    parts = [f"{label:<22} H/Lam={sheets * spacing * q1:6.3f} | b={b.real:+.6f}{b.imag:+.0e}j err={abs(b - b_pred):.1e} |"]
    for N, i_omega in ((1, 0), (2, 1)):
        pred = stack_factor(TAU * N * s1 * spacing, sheets)
        c, _ = fit_scale(F1[i_omega], F[i_omega])
        resid = norm(F[i_omega] - pred * F1[i_omega]) / norm(F1[i_omega])
        parts.append(f"N={N}: B={pred.real:+.6f} c={c.real:+.6f}{c.imag:+.1e}j |c-B|={abs(c - pred):.1e} resid={resid:.1e} |")
    print(" ".join(parts), flush=True)
    return F


print("=== Scan through the minimum at H = Lambda_1, d = 1e-4 lambda_L, fixed total charge")
for sheets in (2, 4, 7, 10, 12, 14, 16, 18, 21, 24, 28):
    report(f"M={sheets}", stack("phase1_scan", f"M{sheets}_d1", sheets, 1e-4), sheets, 1e-4)

print("\n=== Exact minimum, H = Lambda_1 (d = Lambda_1 / M), and its refinement")
exact = {}
for sheets in (14, 28, 56):
    d = 1 / q1 / sheets
    exact[sheets] = report(f"M={sheets} d=Lam/{sheets}", stack("phase1_scan", f"exactmin_M{sheets}", sheets, d), sheets, d)

print("\n=== Spacing refinement at fixed H and charge (d = 1e-4 against 0.5e-4 lambda_L)")
for coarse in (7, 14, 21):
    F_c = read_field(stack("phase1_scan", f"M{coarse}_d1", coarse, 1e-4))[1]
    fine = 2 * coarse
    run = stack("phase1_scan", f"M{fine}_d05", fine, 0.5e-4)
    F_f = report(f"M={fine} d=0.5e-4", run, fine, 0.5e-4)
    for N, i_omega in ((1, 0), (2, 1)):
        predicted = abs(stack_factor(TAU * N * s1 * 0.5e-4, fine) - stack_factor(TAU * N * s1 * 1e-4, coarse))
        measured = norm(F_f[i_omega] - F_c[i_omega]) / norm(F1[i_omega])
        print(f"    H = {coarse}e-4 lambda_L, N={N}: |F_fine - F_coarse| / |F_1| = {measured:.3e}, discrete prediction {predicted:.3e}")

print("\n=== Off-centre stack: beam_center_z shifts the phase by exp(-i s z)")
dz = 1e-3
run = stack("phase1_scan", "M7_shifted", 7, 1e-4, {"beam_center_z": f"{dz:.10e} lambda"})
F = read_field(run)[1]
for N, i_omega in ((1, 0), (2, 1)):
    pred = stack_factor(TAU * N * s1 * 1e-4, 7) * np.exp(-1j * TAU * N * s1 * dz)
    c, _ = fit_scale(F1[i_omega], F[i_omega])
    print(f"N={N}: predicted c = {pred:.6f}, measured c = {c:.6f}, |c-B| = {abs(c - pred):.1e}")

print("\n=== Control, diagnostic only: q d = 2 pi rephases the sheets (d = Lambda_1), one screen point")
point = {"circular_detector_N_R": 1, "circular_detector_N_phi": 1, "N_harmonics": 1}
for sheets in (2, 3, 5):
    b = form_factor_b(stack("phase1_scan", f"rephase_M{sheets}", sheets, 1 / q1, point))
    print(f"M={sheets}: |b| = {abs(b):.9f}, H = {sheets} Lambda_1")

print("\n=== Radiation formulas at the exact minimum and at H = Lambda_1 / 2, one screen point (totals)")
for formula in ("direct", "simplified", "long_distance_direct_approx", "long_distance_simplified_approx"):
    settings = dict(point, radiation_formula=formula)
    flat = read_field(run_variant("phase1_scan", f"point_flat_{formula}", settings))[1][0]
    for sheets, label in ((14, "H=Lam"), (7, "H=Lam/2")):
        d = (1 / q1 / sheets) * (1 if sheets == 14 else 0.5)
        F = read_field(stack("phase1_scan", f"point_{label.replace('/', '')}_{formula}", sheets, d, settings))[1][0]
        pred = stack_factor(TAU * s1 * d, sheets)
        c, _ = fit_scale(flat, F)
        print(f"{formula:<33} {label:<8} B={pred.real:+.2e} c={c.real:+.6e}{c.imag:+.1e}j "
              f"|c-B|={abs(c - pred):.1e} resid={norm(F - pred * flat) / norm(flat):.1e}", flush=True)
