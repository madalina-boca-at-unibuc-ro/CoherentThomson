#!/usr/bin/env python3
"""Phase 3a checkpoint: moved-sheet stacks through the solver, against a flat disk and the cold-beam prediction.

c is the best complex factor in F = c F_flat. The prediction is (1/M) sum_j exp(-i q_N z_j) = (-1)^n J_n(n zeta) exp(+i n phase)
for whole periods at integer n = N ratio: the field's translation phase, the conjugate of b's.

Run from the repository root: .venv/bin/python theory/lattice_stack_tests/phase3a_checkpoint.py
It writes phase3a_results.txt. Finished runs are reused, so delete ~/output_tests/lattice_stack/phase3a after changing the solver.
"""
import numpy as np
from scipy.special import jv

from lattice_stack_common import HERE, fit_scale, form_factor_b, read_field, read_log_value, run_variant
from phase2_density_model import b_from_run

GROUP = "phase3a"
SCATTER = {"plot_beam_scatter": "true", "N_harmonics": 2}
POINT = dict(SCATTER, circular_detector_N_R=1, circular_detector_N_phi=1, beam_lattice_spacing="10 lambda")


def modulated(zeta, sheets, periods, ratio=1, phase=None):
    keys = {"beam_stack_modulation": "true", "beam_modulation_zeta": zeta, "beam_modulation_period_ratio": ratio,
            "beam_modulation_sheets_per_period": sheets, "beam_modulation_periods": periods}
    if phase is not None:
        keys["beam_modulation_phase"] = f"{phase / np.pi:.15f} pi"
    return keys


def predicted(N, zeta, ratio=1, phase=0.0):
    n = N * ratio
    return (-1.0) ** n * jv(n, n * zeta) * np.exp(1j * n * phase) if zeta > 0 else 0.0


def compare(out, label, flat, run, zeta, ratio=1, phase=0.0):
    """One line per harmonic: fitted c against the prediction, and the residual over the screen."""
    field_flat = read_field(flat)[1]
    field = read_field(run)[1]
    for i, N in enumerate((1, 2)):
        c, _ = fit_scale(field_flat[i], field[i])
        pred = predicted(N, zeta, ratio, phase)
        norm = np.linalg.norm(field_flat[i])
        scale_resid = np.linalg.norm(field[i] - c * field_flat[i]) / norm
        pred_resid = np.linalg.norm(field[i] - pred * field_flat[i]) / norm
        out(f"{label:30s} N={N} pred {pred:+.6f}  c {c:+.6f}  |c-pred| {abs(c - pred):.1e}  |c|^2 {abs(c) ** 2:.5f}  "
            f"resid(fit) {scale_resid:.1e}  resid(pred) {pred_resid:.1e}")


def log_check(out, label, run):
    mine, log, count = b_from_run(run)
    out(f"{label:30s} {count} electrons: b from beam file {mine:+.7f}, log {log:+.7f}, difference {abs(mine - log):.1e}")


def main():
    path = HERE / "phase3a_results.txt"
    with open(path, "w") as f:
        def out(line):
            print(line, flush=True)
            f.write(line + "\n")

        out("=== Run A, a0 = 0.1, 1,076 sites per sheet, 16 x 8 screen, long_distance_direct_approx. Total charge fixed.")
        flat = run_variant(GROUP, "flat", SCATTER)
        uniform = run_variant(GROUP, "uniform_S14_P3", dict(SCATTER, **modulated(0, 14, 3)))
        out("\n--- Uniform control and b_1 stacks, 3 periods; Lambda_1 from the log: "
            + read_log_value(uniform, r"Lambda_1 = 2 pi / \|q_axial\|\s+\S+ a\.u\.\s+(\S+) lambda") + " lambda")
        compare(out, "uniform zeta=0 S=14 P=3", flat, uniform, 0)
        m14 = run_variant(GROUP, "mod_zeta1_S14_P3", dict(SCATTER, **modulated(1, 14, 3)))
        m28 = run_variant(GROUP, "mod_zeta1_S28_P3", dict(SCATTER, **modulated(1, 28, 3)))
        compare(out, "zeta=1 S=14 P=3", flat, m14, 1)
        compare(out, "zeta=1 S=28 P=3", flat, m28, 1)
        f14, f28 = read_field(m14)[1], read_field(m28)[1]
        out("Grid change, |F(S=28) - F(S=14)| / |F_flat|: " + ", ".join(
            f"N={i + 1}: {np.linalg.norm(f28[i] - f14[i]) / np.linalg.norm(read_field(flat)[1][i]):.1e}" for i in (0, 1)))

        out("\n--- One period: other zeta, and the phase origin")
        cases = {"zeta=1 S=14 P=1": (1, 14, 1, None), "zeta=1.5 S=14 P=1": (1.5, 14, 1, None), "zeta=1.5 S=28 P=1": (1.5, 28, 1, None),
                 "zeta=1 S=14 P=1 phase=0.7": (1, 14, 1, 0.7)}
        runs = {}
        for label, (zeta, sheets, periods, phase) in cases.items():
            name = label.replace(" ", "_").replace("=", "").replace(".", "p")
            runs[label] = run_variant(GROUP, name, dict(SCATTER, **modulated(zeta, sheets, periods, phase=phase)))
            compare(out, label, flat, runs[label], zeta, phase=phase or 0.0)

        out("\n--- Beam file and q in the log against the log's b")
        log_check(out, "uniform S=14 P=3", uniform)
        log_check(out, "zeta=1 S=14 P=3", m14)
        log_check(out, "zeta=1 S=28 P=3", m28)
        for label, run in runs.items():
            log_check(out, label, run)

        out("\n--- One screen point, 64 sites per sheet")
        flat_point = run_variant(GROUP, "point_flat", POINT)
        for formula in ("direct", "simplified", "long_distance_direct_approx", "long_distance_simplified_approx"):
            f0 = run_variant(GROUP, f"point_flat_{formula}", dict(POINT, radiation_formula=formula))
            run = run_variant(GROUP, f"point_{formula}", dict(POINT, radiation_formula=formula, **modulated(1, 14, 1)))
            compare(out, formula, f0, run, 1)
        run20 = run_variant(GROUP, "point_ratio20_S120", dict(POINT, **modulated(1, 120, 1, ratio=20)))
        compare(out, "b_20: ratio=20 zeta=1 S=120", flat_point, run20, 1, ratio=20)
        run20b = run_variant(GROUP, "point_ratio20_zeta1p111_S120", dict(POINT, **modulated(1.111, 120, 1, ratio=20)))
        compare(out, "b_20: ratio=20 zeta=1.111 S=120", flat_point, run20b, 1.111, ratio=20)

        out("\n--- The period follows a0: a0 = 0.2, same keys (Lambda_1 is computed)")
        flat_a2 = run_variant(GROUP, "point_a02_flat", dict(POINT, laser_a0=0.2))
        run_a2 = run_variant(GROUP, "point_a02_mod", dict(POINT, laser_a0=0.2, **modulated(1, 14, 1)))
        out("Lambda_1 = " + read_log_value(run_a2, r"Lambda_1 = 2 pi / \|q_axial\|\s+\S+ a\.u\.\s+(\S+) lambda") + " lambda")
        compare(out, "a0=0.2 zeta=1 S=14 P=1", flat_a2, run_a2, 1)

        out("\n--- Tilted collision, small plane-wave disk, theta = 0.9 pi, a0 = 0.2: only q_axial can be matched")
        small = dict(POINT, laser_type="plane_wave", beam_cylinder_radius="0.1 lambda", beam_lattice_spacing="0.04 lambda",
                     laser_direction_theta="0.9 pi", laser_a0=0.2)
        flat_t = run_variant(GROUP, "tilted_flat", small)
        run_t = run_variant(GROUP, "tilted_mod", dict(small, **modulated(1, 14, 1)))
        out("q_transverse / |q_axial| = " + read_log_value(run_t, r"q_transverse / \|q_axial\|\s+(\S+)"))
        compare(out, "tilted zeta=1 S=14 P=1", flat_t, run_t, 1)
        log_check(out, "tilted zeta=1 S=14 P=1", run_t)

        out("\n--- Warnings in the log")
        warn = run_variant(GROUP, "point_warn_P1p5_S10", dict(POINT, **modulated(1, 5, 1.5)))
        text = (warn / "run_log.txt").read_text()
        for line in text.splitlines():
            if "WARNING" in line or "periods H / Lambda_mod" in line or line.startswith("  stack sheets"):
                out(line.strip())
        out(f"\nWritten to {path}")


if __name__ == "__main__":
    main()
