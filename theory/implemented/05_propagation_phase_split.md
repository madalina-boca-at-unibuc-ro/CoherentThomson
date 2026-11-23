# 05: Propagation phase split (implemented)

Stub: this spec (`theory/specs/05_propagation_phase_split.md`) was implemented without its own write-up.

- **Status:** done.
- **Commit:** `bed34c5` (2026-10-23).
- **Code:** `radiation.cpp` (phases carry only `e^{ik(R - |X|)}`) and `Simulation::run_simulation` (applies `e^{ik|X|}`).
- **Documented in:** `theory/reference/implementation_details.md` §4.5 (current description, kept in sync with the code); the spec itself holds the measurements.
- **Checked:** All modes converge to 6–7 digits at far screens (spec, measurements).
