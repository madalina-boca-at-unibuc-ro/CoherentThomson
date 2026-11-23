# checks/

Reusable verification scripts and their configs. Run from the repository root. Outputs go to `~/output_tests/...`,
so the plotting scripts' "latest run" in `~/output` is not affected. The checks that need no script (simplified vs
direct, rotation covariance, byte-identity, Python reference) are in CLAUDE.md, Verification.

| Check | Run | Results |
|---|---|---|
| `thomson_dipole/`: the only absolute check (pattern and amplitude), under a second | `./bin/coherent_thomson_solver checks/thomson_dipole/thomson_dipole.cfg`, then `python3 checks/thomson_dipole/thomson_dipole_check.py <run>` | Printed; agrees to ~1e-4 |
| `one_electron/`: one electron on axis (`one_electron_on_axis.cfg`) and its dense on-axis spectrum (`_spectrum.cfg`); seconds | `./bin/coherent_thomson_solver checks/one_electron/<cfg>` | Run folder |
| `nt_scan/`: `trajectory_NT` convergence vs γ, both formulas; ~3 min on 16 threads | `python3 checks/nt_scan/nt_scan.py run` (or without `run` to re-analyze) | Printed; index in `~/output_tests/nt_scan/` |
| `lattice_stack/`: lattice-stack checkpoints (spec 09). They **reuse finished runs in `~/output_tests/lattice_stack/`: delete those after changing the solver** | `python3 checks/lattice_stack/<script>.py`; `rework_checks.py` and `gaussian_stack_checks.py` take `<old_solver> <new_solver>` | `checks/lattice_stack/*_results.txt` |
