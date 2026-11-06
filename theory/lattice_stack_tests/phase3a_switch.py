#!/usr/bin/env python3
"""beam_stack_modulation: each setting matches the run with only the keys it uses. Stack keys are compared across solvers.

Run from the repository root: .venv/bin/python theory/lattice_stack_tests/phase3a_switch.py <solver built before the switch>
"""
import filecmp
import subprocess
import sys
from pathlib import Path

from lattice_stack_common import BASE, ROOT, SOLVER

OUT = Path.home() / "output_tests" / "lattice_stack" / "phase3a_switch"
POINT = {"circular_detector_N_R": 2, "circular_detector_N_phi": 2, "N_harmonics": 2, "beam_lattice_spacing": "10 lambda",
         "plot_beam_scatter": "true"}
STACK = {"beam_stack_sheets": 14, "beam_stack_spacing": "1e-4 lambda"}
MOD = {"beam_modulation_zeta": 1, "beam_modulation_period_ratio": 1, "beam_modulation_sheets_per_period": 14,
       "beam_modulation_periods": 1, "beam_modulation_phase": "0.25 pi"}


def run(solver, name, keys):
    folder = f"output_tests/lattice_stack/phase3a_switch/{name}"
    left, lines = dict(keys, output_folder=folder), []
    for line in BASE.read_text().splitlines():
        key = line.split()[0] if line.split() else ""
        lines.append(f"{key} {left.pop(key)}" if key in left else line)
    lines += [f"{k} {v}" for k, v in left.items()]
    cfg = OUT / "configs" / f"{name}.cfg"
    cfg.parent.mkdir(parents=True, exist_ok=True)
    cfg.write_text("\n".join(lines) + "\n")
    done = subprocess.run([str(solver), str(cfg)], capture_output=True, text=True, cwd=ROOT)
    if done.returncode != 0:
        raise RuntimeError(f"{name}: {done.stdout[-800:]}{done.stderr[-800:]}")
    return sorted((Path.home() / folder).iterdir())[-1]


def same(a, b):
    files = sorted(p.name for p in a.glob("*.dat"))
    return f"{sum(filecmp.cmp(a / f, b / f, shallow=False) for f in files)} of {len(files)} identical"


old = sys.argv[1]
print("true, both sets, against the old solver with only beam_modulation_*:",
      same(run(old, "old_mod", dict(POINT, **MOD)), run(SOLVER, "true_both", dict(POINT, **STACK, **MOD, beam_stack_modulation="true"))))
print("false, both sets, against the old solver with only beam_stack_*:   ",
      same(run(old, "old_stack", dict(POINT, **STACK)), run(SOLVER, "false_both", dict(POINT, **STACK, **MOD, beam_stack_modulation="false"))))
print("no switch and no beam_modulation_*, against the old solver:        ",
      same(run(old, "old_stack2", dict(POINT, **STACK)), run(SOLVER, "none", dict(POINT, **STACK))))
log = next(p for p in (OUT / "true_both").iterdir()) / "run_log.txt"
print("log note:", [line.strip() for line in log.read_text().splitlines() if "are ignored" in line])
