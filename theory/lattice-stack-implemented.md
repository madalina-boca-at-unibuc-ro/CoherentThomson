# Lattice stack: what was implemented

Record of the finished work. The proposal is
[CoherentThomson-lattice-stack-implementation.md](CoherentThomson-lattice-stack-implementation.md); discussion and open
questions are in [lattice_stack-work.md](lattice_stack-work.md).

| Part | Status |
|---|---|
| Preparation: trajectory memory | Done, commit `3052872` |
| Phase 1: uniform stack | Not started |
| Phase 2: density model | Not started |
| Phase 3a: `b₁` | Not started |
| Phase 3b: `b₂₀` | Not started |

## Preparation: trajectory memory (2026-11-06)

**Problem.** Every electron kept its full trajectory until the end of the run, about 2 MB per electron at 20,000 points.
A 280-sheet stack (301,280 electrons) would have needed about 630 GB.

**Change** (commit `3052872`, branch `lattice-stack`):

- `Electron::compute_trajectory` now reserves the storage, instead of the constructor.
- New `Electron::release_trajectory()` frees it. `run_simulation` calls it after each electron's `compute_radiation`.
- Memory now scales with threads, not electrons. Nothing may read a trajectory after `run_simulation`. `electron.dat` is
  written before it, so it is unaffected.

**Check.** The old and new binaries ran the same configs on 4 threads with GCC 15.2. The configs are based on
`config/config_lattice.cfg`, with a 6 × 6 circular screen and 2 harmonics.

| Case | Electrons | `.dat` files identical |
|---|---|---|
| `simplified`, random beam | 60 | 5 of 5 |
| `direct`, random beam | 60 | 5 of 5 |
| `long_distance_simplified`, random beam | 60 | 5 of 5 |
| `long_distance_direct`, random beam | 60 | 5 of 5 |
| `long_distance_simplified_approx`, random beam | 60 | 5 of 5 |
| `long_distance_direct_approx`, random beam | 60 | 5 of 5 |
| `long_distance_simplified_approx`, polar lattice | 3,319 | 5 of 5 |

Peak memory of the lattice case, from `/usr/bin/time`:

| | Peak memory | Wall time |
|---|---|---|
| Before | 3.4 GB | 7.5 s |
| After | 22 MB | 7.0 s |

The outputs are in `~/output_tests/lattice_stack/memory_fix/{baseline,fixed}/`. The fix has not yet been run on the
cluster or at the size of a real stack.
