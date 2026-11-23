# 02: Gaussian-smoothed electron cylinder (implemented)

Stub: this spec (`theory/specs/02_gaussian_smoothed_cylinder.md`) was implemented without its own write-up.

- **Status:** done.
- **Commit:** `ba49aba` (2026-10-09).
- **Code:** `Particle::generate_cylinder_beam` (`src/core/particle/electron_factory.cpp`), key `beam_sigma_z`.
- **Documented in:** CLAUDE.md (Beam along its axis); `py_scripts/particle/plot_beam_density_z.py` plots the analytic density.
- **Checked:** `beam_sigma_z 0` gives byte-identical output to the plain cylinder.
