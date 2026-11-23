#!/usr/bin/env python3
"""Gaussian lattice stacks (beam_stack_sigma_z, theory/lattice-stack-implemented.md "Gaussian stack"): equally spaced
sheets over L = H + 2 cut sigma_z, charges from the smoothed density, then the modulation moving them.

Run from the repository root: python3 theory/lattice_stack_tests/gaussian_stack_checks.py <old_solver> <new_solver>
<old_solver> is a binary built before the change (commit c5c90ee), <new_solver> one built after it. Results are in
gaussian_stack_results.txt. Runs go to ~/output_tests/lattice_stack/gaussian_stack_checks/.

1. sigma_z = 0 is unchanged: every .dat byte-identical to the old solver (lattices, stacks, random beams).
2. Unmodulated Gaussian stacks: the sheet charges in electron_beam_scatter.dat against the density formula, and the
   log's sheet sum against sinc(q H/2) exp(-q^2 sigma_z^2/2).
3. The weights reach the field: on a small screen, the stack's field over the flat disk's, fitted, against the log's
   sheet sum (its conjugate: a translated electron's field gets exp(-i q.r)). A uniform stack is the control.
4. Modulated Gaussian stacks: the log's sheet sum against a fine quadrature of the moved continuous density.
5. Notes, warnings and refusals.
"""
import filecmp
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.special import erfc

from lattice_stack_common import BASE, ROOT, fit_scale, read_field

OUT = Path.home() / "output_tests" / "lattice_stack" / "gaussian_stack_checks"
OLD, NEW = sys.argv[1], sys.argv[2]
# A 2 x 2 screen with a coarse lattice keeps runs to seconds; 16 threads, as these may run on a shared node.
POINT = {"circular_detector_N_R": 2, "circular_detector_N_phi": 2, "N_harmonics": 2, "beam_lattice_spacing": "10 lambda",
         "plot_beam_scatter": "true", "num_threads": 16}
RANDOM = {"beam_distribution": "random", "beam_particle_count": 60, "beam_cylinder_height": "5e-3 lambda"}


def write_config(name, overrides, remove=()):
    lines, left = [], dict(overrides)
    for line in BASE.read_text().splitlines():
        key = line.split()[0] if line.split() else ""
        if key in remove:
            continue
        lines.append(f"{key} {left.pop(key)}" if key in left else line)
    lines += [f"{k} {v}" for k, v in left.items()]
    path = OUT / "configs" / f"{name}.cfg"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def run(solver, tag, name, overrides, remove=(), expect_failure=False):
    """Runs the solver; returns the run folder and stderr, or (None, the error text) for an expected failure."""
    folder = f"output_tests/lattice_stack/gaussian_stack_checks/{tag}/{name}"
    cfg = write_config(f"{tag}_{name}", dict(overrides, output_folder=folder), remove)
    done = subprocess.run([str(solver), str(cfg)], capture_output=True, text=True, cwd=ROOT)
    if expect_failure:
        if done.returncode == 0:
            raise RuntimeError(f"{name}: expected a refusal, but the run finished")
        return None, (done.stdout + done.stderr).strip().splitlines()[-1]
    if done.returncode != 0:
        raise RuntimeError(f"{name}: {done.stdout[-800:]}{done.stderr[-800:]}")
    return sorted((Path.home() / folder).iterdir())[-1], done.stderr


def log(run_dir):
    return (run_dir / "run_log.txt").read_text()


def log_length(run_dir, label):
    """A length from the log's "label  <a.u.> a.u.  <lambda> lambda" line, in lambda."""
    return float(re.search(r"^  " + re.escape(label) + r"\s+\S+ a\.u\.\s+(\S+) lambda", log(run_dir), re.M).group(1))


def density(z, H, sigma):
    """The smoothed density up to a constant, as Particle::sheet_weights evaluates it (independently written)."""
    a = np.abs(z)
    if H == 0.0:
        return np.exp(-a**2 / (2 * sigma**2))
    s = np.sqrt(2) * sigma
    return 0.5 * (erfc((a - H / 2) / s) - erfc((a + H / 2) / s))


def scatter(run_dir):
    """Axial positions (lambda) and charge weights from electron_beam_scatter.dat (beam along +z here)."""
    lines = (run_dir / "electron_beam_scatter.dat").read_text().splitlines()
    header = lines[1].split()
    data = np.array([[float(v) for v in line.split()] for line in lines[2:]])
    w = data[:, header.index("w")] if "w" in header else np.full(len(data), np.nan)
    return data[:, header.index("z")], w


def sheet_sums(run_dir):
    """The log's rows: {N: (sheet sum, smooth value or None)} for unmodulated stacks, or modulated sums."""
    rows = {}
    for N, re_, im_, smooth in re.findall(r"^  (\d)  \((\S+), (\S+)\)  (\S+)  \S+$", log(run_dir), re.M):
        rows[int(N)] = (complex(float(re_), float(im_)), float(smooth))
    for N, re_, im_ in re.findall(r"^  (\d)  \S+  \((\S+), (\S+)\)  none", log(run_dir), re.M):
        rows[int(N)] = (complex(float(re_), float(im_)), None)
    return rows


def stack(H, sigma, M, extra=None):
    """A Gaussian stack in units of lambda_1."""
    return dict(POINT, beam_stack_sheets=M, beam_cylinder_height=f"{H} lambda_1", beam_stack_sigma_z=f"{sigma} lambda_1",
                **(extra or {}))


print("=== 1. sigma_z = 0: old solver against new, every .dat")
SAME = {
    "flat_lattice": POINT,
    "uniform_stack_M14": dict(POINT, beam_stack_sheets=14, beam_cylinder_height="1.4e-3 lambda"),
    "bunched_M14": dict(POINT, beam_stack_sheets=14, beam_cylinder_height="1 lambda_1", beam_modulation="true",
                        beam_modulation_zeta=1, beam_modulation_period="1 lambda_1"),
    "random": dict(POINT, **RANDOM),
    "random_sigma_z": dict(POINT, **RANDOM, beam_sigma_z="2e-3 lambda"),
    "random_modulated": dict(POINT, **RANDOM, beam_sigma_z="2e-3 lambda", beam_modulation="true",
                             beam_modulation_zeta=1, beam_modulation_period="1 lambda_1"),
}
for name, overrides in SAME.items():
    remove = ("beam_lattice_spacing",) if name.startswith("random") else ()
    overrides = {k: v for k, v in overrides.items() if k not in remove}
    old, _ = run(OLD, "old", name, overrides, remove)
    new, _ = run(NEW, "new", name, overrides, remove)
    files = sorted(p.name for p in old.glob("*.dat"))
    same = sum(filecmp.cmp(old / f, new / f, shallow=False) for f in files)
    print(f"  {name:28s} {same} of {len(files)} .dat identical")
# A random beam ignores the stack keys: identical to the old solver's run without them.
old = sorted((OUT.parent / "gaussian_stack_checks" / "old" / "random").iterdir())[-1]
new, _ = run(NEW, "new", "random_with_stack_keys",
             {k: v for k, v in dict(POINT, **RANDOM, beam_stack_sheets=30, beam_stack_sigma_z="0.3 lambda_1",
                                    beam_stack_sigma_z_cut=4).items() if k != "beam_lattice_spacing"},
             ("beam_lattice_spacing",))
files = sorted(p.name for p in old.glob("*.dat"))
print(f"  {'random_with_stack_keys':28s} {sum(filecmp.cmp(old / f, new / f, shallow=False) for f in files)} of "
      f"{len(files)} .dat identical to the old solver's random run")

print("\n=== 2. Unmodulated Gaussian stacks (lengths in lambda_1)")
print("  H     sigma  M    sites  charges vs formula   N  |log sum - smooth|  |log sum - scatter sum|")
CASES = [(0.0, 0.2, 30), (0.0, 0.5, 60), (1.3, 0.1, 30), (1.0, 0.3, 60)]
gaussian_runs = {}
for H, sigma, M in CASES:
    name = f"H{H}_sigma{sigma}_M{M}"
    run_dir, _ = run(NEW, "gaussian", name, stack(H, sigma, M))
    gaussian_runs[name] = run_dir
    lam1 = log_length(run_dir, "Lambda_1 = 2 pi / |q_axial|")
    z, w = scatter(run_dir)
    sites = len(z) // M
    # Sheet m's electrons are consecutive: one charge per sheet, against the formula at the sheet's position.
    z_sheet, w_sheet = z[::sites] / lam1, w[::sites]
    L = H + 2 * 6.0 * sigma
    z_expected = (np.arange(M) - (M - 1) / 2) * L / M
    w_expected = density(z_expected, H, sigma)
    w_expected /= w_expected.sum()
    charge_err = np.max(np.abs(w_sheet - w_expected) / w_expected.max())
    position_err = np.max(np.abs(z_sheet - z_expected))
    for N, (s, smooth) in sorted(sheet_sums(run_dir).items()):
        from_scatter = np.sum(w * np.exp(1j * 2 * np.pi * N * z / lam1)) / w.sum()
        print(f"  {H:<5} {sigma:<6} {M:<4} {sites:<6} {charge_err:.1e} (z {position_err:.0e})   {N}  "
              f"{abs(s - smooth):.2e}            {abs(s - from_scatter):.2e}")

print("\n=== 3. Field: stack over flat disk (fitted on the 2 x 2 screen) against the log's sheet sum")
print("  case                    N  fitted c                    conj(log sum)               |c - conj|  eta")
flat, _ = run(NEW, "gaussian", "flat", POINT)
_, field_flat = read_field(flat)
control, _ = run(NEW, "gaussian", "uniform_M30_H1.3", dict(POINT, beam_stack_sheets=30,
                                                            beam_cylinder_height="1.3 lambda_1"))
lam1 = log_length(control, "Lambda_1 = 2 pi / |q_axial|")
z, _ = scatter(control)
uniform_sum = {N: np.mean(np.exp(1j * 2 * np.pi * N * z / lam1)) for N in (1, 2)}
for name, run_dir in [("uniform_M30_H1.3 (ctrl)", control)] + list(gaussian_runs.items()):
    _, field = read_field(run_dir)
    sums = uniform_sum if run_dir == control else {N: s for N, (s, _) in sheet_sums(run_dir).items()}
    for i, N in enumerate((1, 2)):
        c, eta = fit_scale(field_flat[i].ravel(), field[i].ravel())
        expected = np.conj(sums[N])
        print(f"  {name:23s} {N}  ({c.real:+.6e}, {c.imag:+.6e})  ({expected.real:+.6e}, {expected.imag:+.6e})  "
              f"{abs(c - expected):.1e}     {eta:.1e}")

print("\n=== 4. Modulated Gaussian stacks: log sum against a fine quadrature of the moved density (q_axial > 0 here)")
print("  H    sigma  zeta  phase  M     N  |log sum - quadrature|  |quadrature|")
for H, sigma, zeta, phase, M in [(1.0, 0.1, 1.0, 0.0, 60), (1.0, 0.1, 1.0, 0.5, 60), (0.0, 0.5, 1.0, 0.0, 120),
                                 (2.0, 0.25, 0.7, 0.0, 120)]:
    name = f"mod_H{H}_sigma{sigma}_zeta{zeta}_phase{phase}_M{M}"
    run_dir, _ = run(NEW, "gaussian", name, stack(H, sigma, M, dict(beam_modulation="true", beam_modulation_zeta=zeta,
                                                                    beam_modulation_period="1 lambda_1",
                                                                    beam_modulation_phase=f"{phase} pi")))
    span = H + 2 * 9 * sigma
    u = np.linspace(-span / 2, span / 2, 800001)
    rho = density(u, H, sigma)
    moved = u + zeta / (2 * np.pi) * np.sin(2 * np.pi * u + phase * np.pi)
    for N, (s, _) in sorted(sheet_sums(run_dir).items()):
        reference = np.sum(rho * np.exp(1j * 2 * np.pi * N * moved)) / rho.sum()
        print(f"  {H:<4} {sigma:<6} {zeta:<5} {phase:<6} {M:<5} {N}  {abs(s - reference):.2e}               "
              f"{abs(reference):.3f}")

print("\n=== 5. Notes, warnings and refusals (new solver)")
FLAG = {
    "H0_sigma_M14": stack(0.0, 0.2, 14),
    "one_sheet_with_sigma": stack(0.0, 0.2, 1),
    "beam_sigma_z_on_lattice": dict(POINT, beam_stack_sheets=14, beam_cylinder_height="1 lambda_1",
                                    beam_sigma_z="1e-4 lambda"),
    "coarse_spacing": stack(1.0, 0.1, 12),
    "low_cut": stack(0.0, 0.2, 30, dict(beam_stack_sigma_z_cut=3)),
    "uniform_H0_M14": dict(POINT, beam_stack_sheets=14, beam_cylinder_height="0 lambda"),
}
for name, overrides in FLAG.items():
    run_dir, stderr = run(NEW, "flagged", name, overrides)
    count = re.search(r"particle count\s+(\d+)", log(run_dir)).group(1)
    sheets = re.search(r"stack sheets M\s+(\d+)", log(run_dir)).group(1)
    warnings = [line.removeprefix("Warning: ") for line in stderr.splitlines() if line.startswith("Warning:")]
    print(f"  {name:24s} M = {sheets:>3s}, {count:>5s} electrons; " + (" | ".join(warnings) or "no warning"))
REFUSE = {
    "stack_sigma_negative": stack(0.0, -0.2, 14),
    "cut_zero": stack(0.0, 0.2, 14, dict(beam_stack_sigma_z_cut=0)),
    "stack_sigma_bad_unit_elsewhere": dict(stack(0.0, 0.2, 14), beam_center_z="1 lambda_1"),
}
for name, overrides in REFUSE.items():
    _, error = run(NEW, "refused", name, overrides, expect_failure=True)
    print(f"  {name:32s} {error}")

print("\n=== 6. plot_beam_density_z.py on Gaussian stacks")
for name in ("H0.0_sigma0.5_M60", "H1.3_sigma0.1_M30"):
    done = subprocess.run([sys.executable, str(ROOT / "py_scripts" / "particle" / "plot_beam_density_z.py"),
                           str(gaussian_runs[name])], capture_output=True, text=True, cwd=ROOT,
                          env=dict(os.environ, MPLBACKEND="Agg"))
    lines = [line for line in done.stdout.splitlines() if line.startswith(("N =", "mean", "var", "fraction"))]
    print(f"  {name}: exit {done.returncode}; " + " | ".join(lines) + (done.stderr[-300:] if done.returncode else ""))
