# Repository reorganization plan

> **Status: APPROVED 2026-11-23, carried out the same day** (commits 1–4 below). Kept in `theory/archive/` as the
> record of the moves.

Written 2026-11-13; revised 2026-11-23 against the repository at `4741d58` (three draft specs and the 2π fix note
added since, `AGENTS.md`, `reference_config/` and the rejected phase proposal already gone). No solver behaviour
changes: only files move, and paths in comments, docs and scripts follow.

## Decisions (user)

| Question | Decision |
|---|---|
| Rename `theory/`? | **No**, keep `theory/` (2026-11-13) |
| `dev_notes.md` | **Freeze** it as an archive; start a new log (2026-11-13) |
| Specs | **Number them**, `specs/NN_*.md` + `implemented/NN_*.md`, as in the microbunching project (2026-11-13) |
| Scripts that need solvers that no longer exist | **Delete** (`phase3a_switch.py`, `phase3a_noop_and_rejects.py`; git keeps them) (2026-11-13) |
| `reference_config/` | **Remove** (2026-11-13). Done: it no longer exists. The known-to-work archive comes back in scaled units once `physics/scaling_laws.md` is settled |
| Specs written after the plan | **Numbered 11–13** in dependency order (2026-11-23) |
| `spectral_density_2pi_fix.md` | **`work/`** (2026-11-23) |
| Paper run campaign (spec 12) | **Under `reports/`**, not a top-level `paper/` (rule 4) (2026-11-23) |

## Rules the layout follows

1. **One kind of document per folder,** with the status in the index (`theory/README.md`), not in folder names, so a
   status change never breaks a link.
2. **A spec and its write-up share a number.** Numbers follow the date the spec was first committed (same-day specs:
   dependency order). A spec that was implemented without a write-up gets a short `implemented/NN` stub: status,
   commit, where the code is, where it is documented. Rejected specs are not numbered; they go to `archive/` (or are
   deleted, as the rejected phase proposal was on 2026-11-23).
3. **Living references stay unnumbered** (`implementation_details.md`, `lg_decomposition_reference.md`): they describe
   the code as it is now and are kept in sync with it.
4. **An experiment that fed a report lives with that report** (configs, batch scripts, comparison scripts). `checks/`
   holds only reusable checks.
5. **Frozen configs stay frozen:** report configs are moved, never edited, even where their comments name old paths.
   The report's README lists the path changes once.

## Target layout

```
theory/
  README.md                ← NEW: index of every document with its status
  physics/                 ← formulas and derivations
  reference/               ← living descriptions of the code
  specs/NN_*.md            ← task specs, done or draft
  implemented/NN_*.md      ← what was done for each spec, and how it was checked
  work/                    ← logs, open investigations, fix write-ups
  archive/                 ← superseded, rejected, one-off
  cluster.md               ← srun.txt + slurm_head_node_report.md
checks/                    ← NEW, top level: reusable verification scripts and their configs
reports/2026-10-29_progress/   ← was theory/scientific_report/
reports/2026_thomson_paper/    ← NOT created here: spec 12's campaign folder, once that spec is implemented
config/                    ← unchanged
```

## Moves

All with `git mv`, in one commit with no content changes (§ Commits).

### `theory/physics/`

| From | To | Note |
|---|---|---|
| `FT_Faraday_tensor-direct_and_simplified_forms.md` | `physics/FT_faraday_tensor_v1.md` | **Not archived:** `radiation.cpp` and `simulation.cpp` cite it for the notation, the direct form, the boundary term and the Jacobian. Add a header line: the simplified short-range formula is superseded by v2 |
| `FT_Faraday_tensor-direct_and_simplified_forms-v2.md` | `physics/FT_faraday_tensor_v2.md` | |
| `numerical_calculation_of_angular_momentum.md` | `physics/angular_momentum_numerics.md` | |
| `scaling_laws.md` | `physics/scaling_laws.md` | Keeps its draft notice |

### `theory/reference/`

| From | To |
|---|---|
| `implementation_details.md` | `reference/implementation_details.md` |
| `lg_decomposition_reference.md` | `reference/lg_decomposition_reference.md` |

### `theory/specs/` and `theory/implemented/`

| NN | Spec (from → `specs/NN_*.md`) | Implemented (`implemented/NN_*.md`) | First committed |
|---|---|---|---|
| 01 | `long_distance_direct_simplified_coding_guide.md` → `01_long_distance_formulas.md` | `long_distance_implementation_review.md` → `01_long_distance_formulas.md` | 2026-09-28 |
| 02 | `gaussian_smoothed_electron_cylinder.md` → `02_gaussian_smoothed_cylinder.md` | stub | 2026-10-09 |
| 03 | `table of beam electron radius in units of w_0 for different p,l.md` → `03_lg_beam_radius_table.md` (the PNG `table-with-R_beam.png` → `specs/03_lg_beam_radius_table.png`) | stub (`IoUtils::lg_beam_radius_in_w0`) | 2026-10-23 |
| 04 | `post_processing.md` → `04_radial_phase_removal.md` | stub (`plot_processed.py radial_phase`) | 2026-10-23 |
| 05 | `stabilize_phase_split_proposal.md` → `05_propagation_phase_split.md` | stub (`implementation_details.md` §4.5) | 2026-10-23 |
| 06 | `thomson_form_factor_diagnostic.md` → `06_form_factor_diagnostic.md` | stub (`form_factor.cpp`; the `q` change in the frozen `work/dev_notes_2026-07-28_to_<D>.md`) | 2026-10-23 |
| 07 | `lg_screen_decomposition_simple.md` → `07_lg_screen_decomposition.md` | `lg_decomposition_implementation.md` → `07_lg_screen_decomposition.md` | 2026-10-27 |
| 08 | `gaussian_waist_postprocessing.md` → `08_gaussian_waist_estimate.md` | stub (`plot_observables.py`, `estimate_gaussian_waist`) | 2026-10-27 |
| 09 | `CoherentThomson-lattice-stack-implementation.md` → `09_lattice_stack.md` | `lattice-stack-implemented.md` → `09_lattice_stack.md` | 2026-11-06 |
| 10 | `effective_a0_spec.md` → `10_effective_a0.md` (draft) | — | 2026-11-10 |
| 11 | `scaled_parameters_spec.md` → `11_scaled_parameters.md` (draft) | — | 2026-11-21 |
| 12 | `paper_campaign_spec.md` → `12_paper_campaign.md` (draft; depends on 11) | — | 2026-11-21 |
| 13 | `helicity_plots_spec.md` → `13_helicity_plots.md` (draft) | — | 2026-11-21 |

Next free number: 14.

### `theory/work/`

`<D>` below is the date of the move commit.

| From | To | Note |
|---|---|---|
| `dev_notes.md` | `work/dev_notes_2026-07-28_to_<D>.md` | Frozen. One line added at the top: "Frozen `<D>`; new entries go to `dev_notes.md`". It was last edited 2026-11-21, after this plan was first written |
| — | `work/dev_notes.md` | NEW, empty except a pointer to the frozen file |
| `lattice_stack-work.md` | `work/lattice_stack-work.md` | Spec 09's work log |
| `LG_mode_decomposition.md` | `work/lg_mode_decomposition_plan.md` | The staged plan behind specs 04, 07, 08 |
| `spectral_density_2pi_fix.md` | `work/spectral_density_2pi_fix.md` | Fix write-up (2026-11-21), no spec |
| `future_dev/*.md` (6 files) | `work/future_dev/*.md` | Unchanged names |

### `theory/archive/`

| From | Why |
|---|---|
| `reorganization_plan.md` (this file) | After it is carried out |

### `theory/cluster.md`

`srun.txt` and `slurm_head_node_report.md` concatenated, each under its own heading, text unchanged.

### `checks/`

| From | To |
|---|---|
| `theory/lattice_stack_tests/*` (except the two deleted scripts; not `__pycache__/`) | `checks/lattice_stack/` (`ROOT = HERE.parents[1]` still resolves to the repository) |
| `theory/future_dev/tests/thomson_dipole.cfg`, `thomson_dipole_check.py` | `checks/thomson_dipole/` |
| `theory/future_dev/tests/one_electron_on_axis{,_spectrum}.cfg` | `checks/one_electron/` |
| `theory/future_dev/tests/nt_scan.py` | `checks/nt_scan/nt_scan.py` (its `REPO` goes from three `..` to two) |
| — | `checks/README.md` (NEW: one line per check, how to run it, where its results are) |

### `reports/2026-10-29_progress/`

| From | To |
|---|---|
| `theory/scientific_report/*` | `reports/2026-10-29_progress/*` (same depth for every `run_batch.sh` and `make_lens_figure.py`, so their relative paths still work) |
| `theory/future_dev/tests/bunch_length/{README.md, compare_bunch_length.py, run_batch.sh}` | `reports/2026-10-29_progress/configs/bunch_length/` (next to the frozen configs) |

`compare_lattice_random.py` builds its import path as `parents[1] / 'future_dev' / 'tests' / 'bunch_length'`; after
the move it is `parents[0] / 'configs' / 'bunch_length'`.

## Deletions

| What | Why |
|---|---|
| `theory/future_dev/tests/ponderomotive_lens/*.cfg` (8), then the empty folder | Byte-identical to `configs/lens_gaussian_pulse/` and `configs/lattice_check_flat_top/` in the report (re-checked 2026-11-23) |
| `theory/future_dev/tests/bunch_length/*.cfg` (4) | Byte-identical to the report's `configs/bunch_length/` (re-checked 2026-11-23) |
| `theory/lattice_stack_tests/phase3a_switch.py`, `phase3a_noop_and_rejects.py` | Need solvers that no longer exist (user decision) |

After the moves and deletions `theory/future_dev/` is empty and goes too.

## Edits

1. **Path references.** Update every reference to a moved file: 29 lines in 11 `src/` files (comments only), ~28 in
   `py_scripts/`, the moved scripts, about 200 in Markdown, and the non-frozen configs (`config/`, `checks/`).
   `bunch_length/run_batch.sh` uses `theory/future_dev/tests/bunch_length/*.cfg`, which becomes its own folder.
   The report's `runs_A-D*/run_batch.sh` name `theory/scientific_report/configs/...`, which becomes
   `reports/2026-10-29_progress/configs/...`. `work/future_dev/ponderomotive_lens.md` points at the report's copies
   of the deleted lens configs. Specs 11–13 cite each other by bare filename (`scaled_parameters_spec.md`, ...) and
   `theory/srun.txt`: these follow too.
2. **Spec 12 (draft) follows decision "under `reports/`":** its `paper/` folder becomes `reports/2026_thomson_paper/`
   (`make_configs.py`, `base/base.cfg`, `runs.csv`, `run_batch.sh`, `run_index.tsv`, `README.md`), including its §8
   documentation point. Its mention of `reference_config/` (§1) becomes past tense ("had", removed 2026-11). Nothing
   else in the spec changes; it stays a draft.
3. **`theory/README.md` (new):** one table per folder: file, one-line content, status (current / done / draft /
   superseded). It takes over CLAUDE.md's "Open issues" list.
4. **CLAUDE.md:**
   - replace "Open issues" with a pointer to `theory/README.md`;
   - add a "Spec-driven workflow" section in the microbunching wording: numbered specs, write-ups with the same number,
     drafts carry a notice and wait for approval;
   - update the paths of moved files; drop `ponderomotive_lens/` from the `theory/future_dev/tests/` list (its
     configs are the report's) and point the test-input bullets at `checks/`.
5. **`implemented/` stubs** for 02–06 and 08 (six files, a few lines each).
6. **`v1` header** (FT tensor) and the frozen-log line in `dev_notes`, as in the tables above.

## Commits

0. Done (`4741d58`, 2026-11-23): CLAUDE.md rewrite, Debug build folder, rejected phase proposal deleted, `AGENTS.md`
   removed (untracked; nothing to commit).
1. Pending work: `physics/scaling_laws.md` (still at `theory/scaling_laws.md` then), this plan.
2. `git mv` only (moves and renames), so `git log --follow` keeps every file's history.
3. Deletions.
4. Path references, new files (`theory/README.md`, `checks/README.md`, stubs, new `dev_notes.md`), spec 12's paths,
   CLAUDE.md.

## Verification before commit 4

- `grep` for every old path finds nothing, except inside the frozen report configs and the frozen `dev_notes`.
- Every relative Markdown link resolves (a short script; not committed).
- Clean build. `config/config.cfg` with a small electron count gives byte-identical `.dat` files before and after
  (only comments changed in `src/`).
- `checks/thomson_dipole` passes; `checks/lattice_stack/lattice_stack_common.py` and the report's
  `compare_lattice_random.py` import from their new places.

## Not in this plan (separate decisions)

- **Dropbox ignore for `.venv`, `build/`, `build-debug/`, `bin/`** (`attr -s com.dropbox.ignored -V 1 <dir>`).
  Ignoring a folder **removes it from Dropbox and from the other laptop**, which then has to rebuild its own `.venv`
  and binary. That is the intended state (`-march=native`), but it should be done deliberately, one laptop at a time.
- **The microbunching project's `AGENTS.md`,** which has drifted 55 lines from its CLAUDE.md. Same fix, other
  repository. This repository's `AGENTS.md` was deleted on 2026-11-23; if Dropbox restores it again, delete it again
  rather than syncing it with CLAUDE.md.
- **Known-to-work cases.** Once the scaling note is approved, a `config/reference/` of scaled cases, each with its
  recorded result, as in microbunching's `inputs/reference/`.
