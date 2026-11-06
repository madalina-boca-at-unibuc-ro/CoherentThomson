"""Shared helpers for the lattice-stack tests: run the solver on config variants and read the field back.

Every run is a variant of base_run_A.cfg (run A of the progress report on a small screen), writing to
~/output_tests/lattice_stack/<group>/<name>/. Fields are compared as totals (LR + SR + BR), the only physical sum.
"""
import re
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOLVER = ROOT / "bin" / "coherent_thomson_solver"
BASE = HERE / "base_run_A.cfg"
SCRATCH = Path.home() / "output_tests" / "lattice_stack" / "_configs"
# The independent components of the antisymmetric F^{mu nu}.
COMPONENTS = ["F01", "F02", "F03", "F12", "F13", "F23"]


def run_variant(group, name, overrides):
    """Runs the solver with `overrides` ({key: value}, added if absent) and returns the run folder."""
    overrides = dict(overrides, output_folder=f"output_tests/lattice_stack/{group}/{name}")
    lines = BASE.read_text().splitlines()
    remaining = dict(overrides)
    for i, line in enumerate(lines):
        key = line.split()[0] if line.split() else ""
        if key in remaining:
            lines[i] = f"{key} {remaining.pop(key)}"
    lines += [f"{k} {v}" for k, v in remaining.items()]
    SCRATCH.mkdir(parents=True, exist_ok=True)
    cfg = SCRATCH / f"{group}__{name}.cfg"
    cfg.write_text("\n".join(lines) + "\n")
    parent = Path.home() / overrides["output_folder"]
    # Reuse a finished run with exactly this config, so the analysis can be repeated without recomputing. Delete the
    # folder after changing the solver.
    if parent.exists():
        for old in sorted(p for p in parent.iterdir() if p.is_dir()):
            if (old / "radiation_field.dat").exists() and (old / "config.cfg").read_text() == cfg.read_text():
                return old
    out = subprocess.run([str(SOLVER), str(cfg)], capture_output=True, text=True, cwd=ROOT)
    if out.returncode != 0:
        raise RuntimeError(f"{name} failed:\n{out.stdout[-2000:]}\n{out.stderr[-2000:]}")
    return sorted(p for p in parent.iterdir() if p.is_dir())[-1]


def read_field(run):
    """Returns (omega_list, {omega_index: array [N_screen, 6]}) of the total field in radiation_field.dat."""
    path = Path(run) / "radiation_field.dat"
    with open(path) as f:
        header = next(line for line in f if line.startswith("i_omega")).split()
    data = np.loadtxt(path, comments="#", dtype=str)
    data = data[data[:, 0] != "i_omega"].astype(float)
    col = {name: i for i, name in enumerate(header)}
    total = np.zeros((data.shape[0], 6), dtype=complex)
    for j, comp in enumerate(COMPONENTS):
        for term in ("LR", "SR", "BR"):
            total[:, j] += data[:, col[f"{term}_{comp}_re"]] + 1j * data[:, col[f"{term}_{comp}_im"]]
    i_omega = data[:, 0].astype(int)
    omegas = [data[i_omega == i, 1][0] for i in np.unique(i_omega)]
    return omegas, {int(i): total[i_omega == i] for i in np.unique(i_omega)}


def read_log_value(run, pattern):
    """Returns the first group of `pattern` in the run's run_log.txt."""
    match = re.search(pattern, (Path(run) / "run_log.txt").read_text())
    if not match:
        raise RuntimeError(f"{pattern} not found in {run}/run_log.txt")
    return match.group(1)


def fit_scale(reference, field):
    """Best complex c in field = c * reference, and eta = |field - c ref|^2 / |c ref|^2."""
    c = np.vdot(reference, field) / np.vdot(reference, reference)
    residual = field - c * reference
    return c, np.vdot(residual, residual).real / max(np.vdot(c * reference, c * reference).real, 1e-300)


def laser_kL_units(run):
    """Returns k_s/k_L and (k_s - k_L)/k_L at harmonic 1, head-on: the run log's k_s, and k_s + 1."""
    k_s = float(read_log_value(run, r"k_s\s+\(.*?\[\(\S+, \S+, ([-0-9.e+]+)\) k_L\]"))
    return k_s, k_s + 1.0


def form_factor_b(run):
    """Returns the run log's complex form factor b at harmonic 1."""
    re_part, im_part = read_log_value(run, r"b \(re, im\)\s+\((\S+, \S+)\)").split(", ")
    return complex(float(re_part), float(im_part))


def stack_factor(q_d, sheets):
    """B_M = (1/M) sum_j exp(i q z_j) of M sheets of spacing d, from the finite sum (no removable singularities)."""
    return np.mean(np.exp(1j * q_d * (np.arange(sheets) - (sheets - 1) / 2)))
