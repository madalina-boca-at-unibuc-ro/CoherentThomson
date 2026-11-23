# Paper run campaign: generated, frozen configs: spec

> **Status: DRAFT, NOT APPROVED FOR IMPLEMENTATION.** Written 2026-11-21. Discuss first; implement only after the user
> approves, then stop at the checkpoint (§9). Depends on `theory/specs/11_scaled_parameters.md` (§7 below adds one
> function to it).

## 1. Purpose

The paper needs about 45 runs (`~/Dropbox/pub/ms/2026-Thomson/runs_sec4_plan.md`: groups M, S, L, T, C) that differ
in about 8 keys. Hand-edited configs had already drifted: `reference_config/` (removed since) had σ = 25 instead of
30, an anti-aligned helicity and a screen that was not at the paper's $Z/Z_*$.

The configs therefore become **generated and frozen**:
- one template holds the shared choices;
- one table holds what varies between runs;
- a script writes ordinary, complete config files;
- the files are committed and never edited by hand.

The solver and the config format do not change. The configs stay in physical units.

## 2. Layout

```text
reports/2026_thomson_paper/
├── README.md            # what each group is for, how to (re)generate and run, where outputs and the index live
├── base/base.cfg        # every key, physical units, the shared choices (§3)
├── runs.csv             # one row per run: what varies (§4)
├── make_configs.py      # base + row -> configs/<group>/<id>.cfg, then the checks (§5)
├── configs/{M,S,L,T,C}/ # generated, committed; each file says so in its header
├── run_batch.sh         # build once on the node, run a list of ids or groups, post-process, index (§6)
└── run_index.tsv        # committed copy of the output index after each batch (§6)
```

- **Placement:** `reports/2026_thomson_paper/`, with the paper it feeds (the 2026-11-23 reorganization: an experiment
  that feeds a report lives with that report), not in `theory/`.
- **`config/` stays as it is**, for generic and experimental configs.
- **`reference_config/`** (untracked) was to become `base/base.cfg` after the fixes in §3, but it was removed before
  this spec was implemented. **Open:** what `base/base.cfg` starts from instead (§3 takes every key from
  `config/config_lattice.cfg`, replaced on 2026-11-23 by `config/reference/lattice_2d.cfg`, a candidate).

## 3. `base/base.cfg`

- **Every key, as in `config/reference/lattice_2d.cfg` (was `config/config_lattice.cfg`), with its comments.** Keys that `runs.csv` sets get the placeholder
  value `@ROW` (§5), so a key that is forgotten in a row cannot silently take a stale value.
- **Shared choices**, from the run plan:
  - envelope: `laser_flat_duration 0`, cutoff 4;
  - beam: `beam_distribution polar_lattice`, `beam_cylinder_radius -1` (LG table);
  - formula: `radiation_formula long_distance_direct_approx`;
  - screen: `detector_type circular`, `N_R 128`, `N_phi 64`;
  - laser along −z, screen along +z;
  - `dense_frequency_spectrum false`, `plot_beam_scatter false`, `plot_field_heatmap false`;
  - `require_clean_build true`.
- **Sign convention: aligned spin and OAM** ($m$, $\sigma$ of the same sign in the lab), as in the report runs:
  `laser_zeta_1 (1, 0)`, `laser_zeta_2 (0, +1)`, with `laser_lg_l` ≥ 0. The anti-aligned rows flip `zeta_2` through
  their `sign` column.

## 4. `runs.csv`

One row per run. Comment lines start with `#`. The columns:

| Column | Meaning | Becomes |
|---|---|---|
| `id` | e.g. `M05` | file name and output folder |
| `group` | M, S, L, T, C | subfolder |
| `purpose` | free text, one line | header comment |
| `lambda_L_um` | 10 or 0.8 | `laser_frequency` (a.u., from λ) |
| `p_z_mc` | 13.483, 3.823, 0 | `average_pz` |
| `mode` | `LG00`, `LG01`, `LG02`, `LG22` | `laser_lg_p`, `laser_lg_l` |
| `sign` | `aligned` or `anti` | `laser_zeta_2_im` = +1 or −1 |
| `a0`, `N_c` | | `laser_a0`, `laser_wing_sigma … cycles_adim` |
| `screen` | `L`, `H`, or `L2` (C6: $Z/Z_*$ doubled) | `circular_detector_distance` (computed, §5), `R_max` |
| `d_lambda` | lattice spacing in λ | `beam_lattice_spacing` |
| `N_min`, `N_count` | harmonics | `N_harmonics_min`, `N_harmonics` |
| `NT` | | `trajectory_NT` |
| `point_theta` | S group only: screen point in $\theta/\theta_*$ (empty: full screen) | `N_R 1`, `N_phi 1`, `R_min = R_max` |
| `spectrum` | S group only: `lo:hi:n` in $\omega/\omega_{1,\rm ref}$ | `dense_frequency_spectrum true`, `omega_min/max`, `N_omega` |
| `extra` | `key=value;…` for anything else (e.g. C5: `beam_distribution=random;beam_particle_count=300000`) | applied last |

The **screen distance is never typed.** The generator computes it from the target $Z/Z_*$ and the row's physics, so
it cannot drift from the definitions.

The first version of `runs.csv` holds every run of the plan:
- M1–M20, plus two anti-aligned rows (LG01 and LG02 at $a_0=0.26$, $N_c=30$);
- S: the point spectra at 6–8 radii for S1–S4, as their first step;
- L1–L4, T1–T5, C1–C6.

## 5. `make_configs.py`

```
python3 reports/2026_thomson_paper/make_configs.py [--only ID,ID,...] [--check] [--list]
```

1. **Read** `base.cfg` and `runs.csv`. Validate every row: known columns, numeric fields, a known mode and screen,
   unique ids.
2. **Compute the derived keys:**
   - `laser_frequency` $=2\pi c/\lambda_L$ in a.u.;
   - the screen: $Z=(Z/Z_*)\,Z_*$ with $Z_*=\tfrac12k_{1,\rm ref}w_0^2$ and
     $k_{1,\rm ref}=k_LD^2/(1+a_0^2/2)$, written in `w0`; `R_max` from $\theta_{\max}/\theta_*$;
   - for point spectra, $R=Z\tan(\theta\,\theta_*)$.

   All of these come from `scaled_parameters.py` (§7), not from a second copy of the formulas.
3. **Substitute** every `@ROW` key and apply `extra`. Refuse if any `@ROW` is left, or if `extra` names an unknown
   key (one that `base.cfg` lacks).
4. **Write** `configs/<group>/<id>_<mode>_a<a0>_Nc<N_c>[_<lambda>][_anti].cfg`, e.g. `M05_LG01_a026_Nc30.cfg`, with:
   - `output_folder output/paper_2026/<group>/<id>`, so the timestamped run folder lands inside;
   - a header naming the id, the purpose, the source row and "generated by
     `reports/2026_thomson_paper/make_configs.py`, do not edit".

   The output is deterministic (no dates), so regenerating is a no-op unless the inputs changed.
5. **Check** each written file with `scaled_parameters.from_config` (§7) against the row's screen group
   (`--expect`), plus aliasing, envelope and sign. A failed check deletes the file and exits non-zero.
6. **Options:**
   - `--check` regenerates into a temporary folder and diffs against `configs/`. Any difference means a hand edit or
     an uncommitted change to the inputs.
   - `--list` prints the table of ids with their scaled parameters. It is the draft of the configuration table and
     of S4.

## 6. `run_batch.sh`

```
bash reports/2026_thomson_paper/run_batch.sh <ID|GROUP> [...]        # inside the compute-node job (theory/cluster.md)
```

1. **Build once,** clean, on the node. It refuses if the git state is dirty, as `require_clean_build` would anyway,
   but before any run starts.
2. **For each config:** run the solver; the run folder is the single new folder under
   `~/output/paper_2026/<group>/<id>/`.
3. **Post-process each run:**
   - `scaled_parameters.py --write-log --expect <screen>`;
   - `plot_observables.py`;
   - `lg_decomposition.py --scan-z0`;
   - the helicity plots (`theory/specs/13_helicity_plots.md`).

   A post-processing failure is logged but does not stop the batch.
4. **Append one line per run** to `~/output/paper_2026/index.tsv`: id, config file, config SHA-256, run folder, git
   commit, host, start, end, solver exit code, post-processing status. After each batch, copy it to
   `reports/2026_thomson_paper/run_index.tsv` and commit it. That makes the manifest from paper results to runs.

## 7. Addition to `theory/specs/11_scaled_parameters.md`

Add `from_config(config_path)`, which computes the same quantities from a config before any run. It uses:
- the requested lattice spacing, not the adjusted one;
- the LG-table radius from the same table the solver uses (`IoUtils::lg_beam_radius_in_w0`; Python copy kept in sync
  as the other C++/Python copies are).

`scaled_parameters(run_dir)` stays the authority after the run. Comparing the two flags any difference between what
was requested and what ran: for example, the lattice spacing adjustment, about 1%.

## 8. Documentation

- **`CLAUDE.md`:** the `reports/2026_thomson_paper/` folder and its rule (configs generated, never edited; change
  `runs.csv` or `base.cfg` and regenerate), next to the `config/` rules.
- **`reports/2026_thomson_paper/README.md`:** the workflow (edit the table, generate, commit, run the batch, commit the
  index) and the output layout.

## 9. Checkpoint

Implement §2–§7. Then stop and show:
- `runs.csv`;
- `make_configs.py --list`;
- the diff of one generated config against `base.cfg`;
- that `--check` is clean;
- one cheap end-to-end run (M10: LG00, $a_0=0.2$, $N_c=8$) through `run_batch.sh`, with its `index.tsv` line and
  post-processing.

## 10. Open points

- **The output root name:** `output/paper_2026/`.
- **Whether S rows** should later switch from points to the radial line (after the solver change in the run plan): a
  new column value, `point_theta` = `line`.
- **Whether `base.cfg`'s comments** are kept in full, or trimmed to keep the generated files readable.
