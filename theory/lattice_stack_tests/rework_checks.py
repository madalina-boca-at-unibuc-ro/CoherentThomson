#!/usr/bin/env python3
"""Beam distribution rework (theory/lattice_stack-work.md): the new layout against the old solver, the modulated random
beam against the cold-beam value, and the refused or corrected configs.

Run from the repository root: python3 theory/lattice_stack_tests/rework_checks.py <old_solver> <new_solver>
<old_solver> is a binary built before the rework (commit 35b3a1a), <new_solver> one built after it. Results are in
rework_results.txt. Runs go to ~/output_tests/lattice_stack/rework_checks/.
"""
import filecmp
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

from lattice_stack_common import BASE, ROOT, read_field

OUT = Path.home() / "output_tests" / "lattice_stack" / "rework_checks"
OLD, NEW = sys.argv[1], sys.argv[2]
# One screen point would hide a wrong transverse layout; a 2 x 2 screen with a coarse lattice keeps runs to seconds.
POINT = {"circular_detector_N_R": 2, "circular_detector_N_phi": 2, "N_harmonics": 2, "beam_lattice_spacing": "10 lambda",
         "plot_beam_scatter": "true"}
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
    folder = f"output_tests/lattice_stack/rework_checks/{tag}/{name}"
    cfg = write_config(f"{tag}_{name}", dict(overrides, output_folder=folder), remove)
    done = subprocess.run([str(solver), str(cfg)], capture_output=True, text=True, cwd=ROOT)
    if expect_failure:
        if done.returncode == 0:
            raise RuntimeError(f"{name}: expected a refusal, but the run finished")
        return None, (done.stdout + done.stderr).strip().splitlines()[-1]
    if done.returncode != 0:
        raise RuntimeError(f"{name}: {done.stdout[-800:]}{done.stderr[-800:]}")
    return sorted((Path.home() / folder).iterdir())[-1], done.stderr


def compare(old, new):
    """Identical .dat files, and the largest |new - old| / max |old| of the total field."""
    files = sorted(p.name for p in old.glob("*.dat"))
    same = sum(filecmp.cmp(old / f, new / f, shallow=False) for f in files)
    _, field_old = read_field(old)
    _, field_new = read_field(new)
    worst = max(np.abs(field_new[i] - field_old[i]).max() / np.abs(field_old[i]).max() for i in field_old)
    return f"{same} of {len(files)} .dat identical, field max rel. diff {worst:.1e}"


def log(run_dir):
    return (run_dir / "run_log.txt").read_text()


print("=== 1. Old solver against new, same beam")
pairs = {
    "flat_lattice": (POINT, POINT),
    "uniform_stack_M14": (dict(POINT, beam_stack_sheets=14, beam_stack_spacing="1e-4 lambda"),
                          dict(POINT, beam_stack_sheets=14, beam_cylinder_height="1.4e-3 lambda")),
    "bunched_M14_P1": (dict(POINT, beam_stack_modulation="true", beam_modulation_zeta=1, beam_modulation_period_ratio=1,
                            beam_modulation_sheets_per_period=14, beam_modulation_periods=1),
                       dict(POINT, beam_stack_sheets=14, beam_cylinder_height="1 lambda_1", beam_modulation="true",
                            beam_modulation_zeta=1, beam_modulation_period="1 lambda_1")),
    "bunched_M28_P2_phase": (dict(POINT, beam_stack_modulation="true", beam_modulation_zeta=1.111,
                                  beam_modulation_period_ratio=1, beam_modulation_sheets_per_period=14,
                                  beam_modulation_periods=2, beam_modulation_phase="1.0 pi"),
                             dict(POINT, beam_stack_sheets=28, beam_cylinder_height="2 lambda_1", beam_modulation="true",
                                  beam_modulation_zeta=1.111, beam_modulation_period="1 lambda_1",
                                  beam_modulation_phase="1.0 pi")),
}
for formula in ("simplified", "direct", "long_distance_direct_approx"):
    both = dict(POINT, **RANDOM, radiation_formula=formula)
    pairs[f"random_{formula}"] = (both, both)
both = dict(POINT, **RANDOM, beam_sigma_z="2e-3 lambda")
pairs["random_sigma_z"] = (both, both)
for name, (old_overrides, new_overrides) in pairs.items():
    remove = ("beam_lattice_spacing",) if name.startswith("random") else ()
    old, _ = run(OLD, "old", name, {k: v for k, v in old_overrides.items() if k not in remove}, remove)
    new, _ = run(NEW, "new", name, {k: v for k, v in new_overrides.items() if k not in remove}, remove)
    print(f"  {name:28s} {compare(old, new)}")

print("\n=== 2. Modulated random beam: electron sum against the cold-beam value (zeta = 1, H = Lambda_mod = Lambda_1)")
print("  N_e    seed   N  |sum - cold|   1/sqrt(N_e)")
for count in (2000, 8000):
    for seed in (532, 533):
        overrides = dict(POINT, beam_distribution="random", beam_particle_count=count, beam_cylinder_height="1 lambda_1",
                         beam_modulation="true", beam_modulation_zeta=1, beam_modulation_period="1 lambda_1",
                         random_seed=seed, circular_detector_N_R=1, circular_detector_N_phi=1)
        overrides.pop("beam_lattice_spacing")
        new, _ = run(NEW, "random_modulated", f"N{count}_seed{seed}", overrides, ("beam_lattice_spacing",))
        rows = re.findall(r"^  (\d)  \S+  \((\S+), (\S+)\)  \((\S+), (\S+)\)", log(new), re.M)
        for N, s_re, s_im, c_re, c_im in rows:
            diff = abs(complex(float(s_re), float(s_im)) - complex(float(c_re), float(c_im)))
            print(f"  {count:5d}  {seed}   {N}  {diff:.3e}      {count ** -0.5:.3e}")

print("\n=== 3. Refused configs (new solver)")
GOOD = dict(POINT, beam_stack_sheets=14, beam_cylinder_height="1 lambda_1", beam_modulation="true",
            beam_modulation_zeta=1, beam_modulation_period="1 lambda_1")
REFUSE = {
    "removed_stack_spacing": dict(POINT, beam_stack_sheets=14, beam_stack_spacing="1e-4 lambda"),
    "removed_stack_modulation": dict(GOOD, beam_stack_modulation="true"),
    "removed_sheets_per_period": dict(GOOD, beam_modulation_sheets_per_period=14),
    "removed_periods": dict(GOOD, beam_modulation_periods=1),
    "removed_period_ratio": dict(GOOD, beam_modulation_period_ratio=1),
    "sigma_z_on_lattice": dict(GOOD, beam_sigma_z="1e-4 lambda"),
    "height_negative": dict(GOOD, beam_cylinder_height="-1 lambda_1"),
    "lambda_1_elsewhere": dict(GOOD, beam_center_z="1 lambda_1"),
    "period_without_unit": dict(GOOD, beam_modulation_period="1"),
    "period_zero": dict(GOOD, beam_modulation_period="0 lambda_1"),
    "zeta_negative": dict(GOOD, beam_modulation_zeta=-1),
    "switch_not_bool": dict(GOOD, beam_modulation="yes"),
    "zeta_without_switch": dict(POINT, beam_modulation_zeta=1),
    "sheets_not_integer": dict(GOOD, beam_stack_sheets=14.5),
    # Beam along +x, laser and screen along z: q has no component along the beam axis.
    "lambda_1_q_perpendicular": dict(POINT, average_px="13.483 mc", average_pz="0.0 mc",
                                     beam_cylinder_height="1 lambda_1", beam_stack_sheets=2),
}
for name, overrides in REFUSE.items():
    _, error = run(NEW, "refused", name, overrides, expect_failure=True)
    print(f"  {name:28s} {error}")

print("\n=== 4. Corrected or flagged configs (new solver): electrons, and the warning")
FLAG = {
    "one_sheet_with_width": dict(POINT, beam_stack_sheets=1, beam_cylinder_height="1 lambda"),
    "sheets_without_width": dict(POINT, beam_stack_sheets=14, beam_cylinder_height="0 lambda"),
    "modulated_single_sheet": dict(POINT, beam_modulation="true", beam_modulation_zeta=1,
                                   beam_modulation_period="1 lambda_1"),
    "modulated_flat_random": dict(POINT, beam_distribution="random", beam_particle_count=60,
                                  beam_cylinder_height="0 lambda", beam_modulation="true", beam_modulation_zeta=1,
                                  beam_modulation_period="1 lambda_1"),
    "non_whole_periods": dict(GOOD, beam_cylinder_height="1.5 lambda_1"),
    "too_few_sheets": dict(GOOD, beam_stack_sheets=4),
    # A random beam ignores the lattice keys (user, 2026-11-09): 60 electrons, no error.
    "stack_keys_on_random": dict(POINT, **RANDOM, beam_stack_sheets=14),
}
for name, overrides in FLAG.items():
    remove = ("beam_lattice_spacing",) if "random" in name else ()
    new, stderr = run(NEW, "flagged", name, {k: v for k, v in overrides.items() if k not in remove}, remove)
    count = re.search(r"particle count\s+(\d+)", log(new)).group(1)
    warnings = [line.removeprefix("Warning: ") for line in stderr.splitlines() if line.startswith("Warning:")]
    print(f"  {name:28s} {count:>6s} electrons; " + (" | ".join(warnings) or "no warning"))
