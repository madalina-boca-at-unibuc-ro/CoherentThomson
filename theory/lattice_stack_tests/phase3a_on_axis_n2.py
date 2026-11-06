#!/usr/bin/env python3
"""Why one on-axis point misleads at N = 2: the flat field there is close to zero. Compares it with the screen, for a b_20 stack.

Run from the repository root: .venv/bin/python theory/lattice_stack_tests/phase3a_on_axis_n2.py
"""
import numpy as np

from lattice_stack_common import fit_scale, read_field, run_variant
from phase3a_checkpoint import GROUP, POINT, SCATTER, modulated, predicted

SPARSE = dict(SCATTER, beam_lattice_spacing="10 lambda")
for name, keys in (("point", POINT), ("screen", SPARSE)):
    flat = read_field(run_variant(GROUP, f"{name}_flat_check", keys))[1]
    print(f"--- {name}: |F_flat|/sqrt(values), N=1 and N=2: {np.sqrt(np.mean(np.abs(flat[0]) ** 2)):.3e}, {np.sqrt(np.mean(np.abs(flat[1]) ** 2)):.3e}")
    for zeta, sheets in ((1, 120), (1.111, 120), (1, 240), (1.111, 240)):
        label = f"{name}_r20_z{str(zeta).replace('.', 'p')}_S{sheets}"
        field = read_field(run_variant(GROUP, label, dict(keys, **modulated(zeta, sheets, 1, ratio=20))))[1]
        for i, N in enumerate((1, 2)):
            pred = predicted(N, zeta, 20)
            c, _ = fit_scale(flat[i], field[i])
            error = np.linalg.norm(field[i] - pred * flat[i])
            print(f"zeta={zeta} S={sheets} N={N}: pred {pred:+.6f} c {c:+.6f} |c-pred| {abs(c - pred):.1e}  "
                  f"|F - pred F_flat|: {error:.2e} = {error / np.linalg.norm(flat[i]):.1e} of |F_flat|")
