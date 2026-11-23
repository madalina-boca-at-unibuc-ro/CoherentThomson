# 06: Form-factor diagnostic (implemented)

Stub: this spec (`theory/specs/06_form_factor_diagnostic.md`) was implemented without its own write-up.

- **Status:** done.
- **Commit:** `bed34c5` (2026-10-23).
- **Code:** `Simulation::compute_form_factor` (`src/core/simulation/form_factor.cpp`), written to `run_log.txt`.
- **Documented in:** `theory/README.md` (Open issues: `b(q)`). The code uses `q = k_s − (K_s·P / K_L·P) k_L` (measured), not the spec's `k_s − k_L`: the change and its measurement are in `theory/work/dev_notes_2026-07-28_to_2026-11-23.md` (2026-11-06) and `theory/implemented/09_lattice_stack.md`.
