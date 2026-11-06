#!/usr/bin/env python3
"""Phase 3a: beams without modulation stay byte-identical, and bad modulation configs are refused.

Run from the repository root: .venv/bin/python theory/lattice_stack_tests/phase3a_noop_and_rejects.py <baseline_solver>
<baseline_solver> is a binary built before the change.
"""
import filecmp
import subprocess
import sys
from pathlib import Path

from lattice_stack_common import BASE, ROOT, SOLVER

OUT = Path.home() / "output_tests" / "lattice_stack" / "phase3a_noop"
POINT = {"circular_detector_N_R": 2, "circular_detector_N_phi": 2, "N_harmonics": 2, "beam_lattice_spacing": "10 lambda",
         "plot_beam_scatter": "true"}


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


def run(solver, name, tag, overrides, remove=()):
    folder = f"output_tests/lattice_stack/phase3a_noop/{tag}/{name}"
    cfg = write_config(f"{tag}_{name}", dict(overrides, output_folder=folder), remove)
    done = subprocess.run([str(solver), str(cfg)], capture_output=True, text=True, cwd=ROOT)
    if done.returncode != 0:
        raise RuntimeError(f"{name}: {done.stdout[-800:]}{done.stderr[-800:]}")
    return sorted((Path.home() / folder).iterdir())[-1]


cases = {"flat_lattice": POINT, "stack_M14": dict(POINT, beam_stack_sheets=14, beam_stack_spacing="1e-4 lambda")}
for formula in ("simplified", "direct", "long_distance_direct_approx"):
    cases[f"random_{formula}"] = dict(POINT, beam_distribution="random", beam_particle_count=60,
                                      beam_cylinder_height="5e-3 lambda", radiation_formula=formula)
print("Byte-identity against the pre-change binary (.dat files):")
for name, overrides in cases.items():
    remove = ("beam_lattice_spacing",) if "random" in name else ()
    overrides = {k: v for k, v in overrides.items() if k not in remove}
    old, new = run(sys.argv[1], name, "old", overrides, remove), run(SOLVER, name, "new", overrides, remove)
    files = sorted(p.name for p in old.glob("*.dat"))
    same = sum(filecmp.cmp(old / f, new / f, shallow=False) for f in files)
    print(f"  {name:36s} {same} of {len(files)} identical")

GOOD = dict(POINT, beam_stack_modulation="true", beam_modulation_zeta=1, beam_modulation_period_ratio=1, beam_modulation_sheets_per_period=14,
            beam_modulation_periods=1)
REJECT = {
    "zeta_negative": dict(GOOD, beam_modulation_zeta=-1),
    "ratio_zero": dict(GOOD, beam_modulation_period_ratio=0),
    "sheets_zero": dict(GOOD, beam_modulation_sheets_per_period=0),
    "sheets_not_integer": dict(GOOD, beam_modulation_sheets_per_period=14.5),
    "periods_zero": dict(GOOD, beam_modulation_periods=0),
    "switch_not_boolean": dict(GOOD, beam_stack_modulation="yes"),
    "zeta_without_switch": {k: v for k, v in GOOD.items() if k != "beam_stack_modulation"},
    "random_beam": dict(GOOD, beam_distribution="random", beam_particle_count=60),
    "ratio_without_zeta": dict(POINT, beam_modulation_period_ratio=1),
    "phase_without_zeta": dict(POINT, beam_modulation_phase="0 pi"),
}
print("Configs that must stop with an error:")
for name, overrides in REJECT.items():
    cfg = write_config(f"reject_{name}", dict(overrides, output_folder="output_tests/lattice_stack/phase3a_noop/reject"))
    done = subprocess.run([str(SOLVER), str(cfg)], capture_output=True, text=True, cwd=ROOT)
    text = (done.stdout + done.stderr).splitlines()
    message = next((line for line in reversed(text) if line.startswith("Error")), text[-1] if text else "")
    print(f"  {name:22s} exit {done.returncode}: {message[:110]}")
missing = dict(GOOD)
for key in ("beam_modulation_period_ratio", "beam_modulation_sheets_per_period", "beam_modulation_periods"):
    cfg = write_config(f"reject_missing_{key}", {k: v for k, v in dict(GOOD, output_folder="output_tests/lattice_stack/phase3a_noop/reject").items() if k != key})
    done = subprocess.run([str(SOLVER), str(cfg)], capture_output=True, text=True, cwd=ROOT)
    print(f"  missing {key:34s} exit {done.returncode}")
