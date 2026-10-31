# Scientific report (progress report, 2026-10-29; Part II added 2026-10-31)

`progress_report.md` is the report for the group; `progress_report.pdf` is the version to send (A5 pages, readable on
a phone). Notes for rebuilding it:

## Runs and configs

Every run the report quotes is in `~/output/progress_report_20261029/` (moved there on 2026-10-31 from `~/output/`
and `~/output_tests/{ponderomotive_lens,bunch_length}/`; run folder names unchanged), in four groups. The configs
are copied into [configs/](configs/) with the same groups; each is byte-identical to the run's own `config.cfg` and
to the original in `config/expansion_a0_scan_server/` or `theory/future_dev/tests/{ponderomotive_lens,bunch_length}/`.
The output folder's `README.md` repeats this manifest with the batch logs, and its `SHA256SUMS` covers the outputs.

The Part I builds were dirty; each of those run folders' `source.patch` applies to the recorded commit. The Part II
builds are clean.

**`runs_A-D/`**: Sections 1–4.3, Figs. 1–7. Commit `5176781` (dirty only in `CLAUDE.md`), job 30438 on
`hpc-node02`, 65536 electrons, ~37 min each on 256 CPUs.

| Report label | Case | Run | Config (`configs/runs_A-D/`) |
|---|---|---|---|
| A | LG(0,0), a0 = 0.1, R_max 15 w0 | `20261029_132340` | `config_expansion_p0l0_a0_01_R15.cfg` |
| B | LG(0,0), a0 = 0.2, R_max 25 w0 | `20261029_143739` | `config_expansion_p0l0_a0_02_R25.cfg` |
| C | LG(2,2), a0 = 0.1, R_max 15 w0 | `20261029_140030` | `config_expansion_p2l2_a0_01_R15.cfg` |
| D | LG(2,2), a0 = 0.2, R_max 25 w0 | `20261029_151457` | `config_expansion_p2l2_a0_02_R25.cfg` |

**`lens_gaussian_pulse/`**: Sec. 1 (lattice of 1076 electrons, halving the spacing changes the field ~1%), the
Sec. 4.4 table and Fig. 8. Electrons at rest, polar lattice, the pulse of runs A–D. Commit `74539fd`, job 30464.

| Case | Run | Config (`configs/lens_gaussian_pulse/`) |
|---|---|---|
| a0 = 0.01 | `20261030_111837` | `gaussian_lattice_a0_001.cfg` |
| a0 = 0.1 (Fig. 8) | `20261030_111559` | `gaussian_lattice_a0_01.cfg` |
| a0 = 0.2 | `20261030_112109` | `gaussian_lattice_a0_02.cfg` |
| a0 = 0.2, spacing halved | `20261030_112351` | `gaussian_lattice_a0_02_fine.cfg` |

**`lattice_check_flat_top/`**: Sec. 1, "the lattice field agreed with a 65 536-electron random run" (`092541` vs
`093318`, a0 = 0.1). Electrons at rest, flat-top pulse. Commit `5176781` plus the then-uncommitted polar-lattice
code (later `465c0bd`), job 30462.

| Case | Run | Config (`configs/lattice_check_flat_top/`) |
|---|---|---|
| lattice, spacing 2.44 λ, a0 = 0.1 | `20261030_092541` | `lattice_a0_01.cfg` |
| lattice, spacing 1.22 λ, a0 = 0.1 | `20261030_092940` | `lattice_a0_01_fine.cfg` |
| lattice, spacing 2.44 λ, a0 = 0.01 | `20261030_092744` | `lattice_a0_001.cfg` |
| random, 65536 electrons, a0 = 0.1 | `20261030_093318` | `random_a0_01.cfg` |

Results and discussion of both lens groups: `theory/future_dev/ponderomotive_lens.md` (Results).

**`bunch_length/`**: Part II (Sec. 5, Fig. 9). Run B's config with a Gaussian bunch length `beam_sigma_z`, all four
on one transverse sample. Commit `088b7ad` (clean), `hpc-node02`, ~37 min each. Logs and the config → run index:
`logs/bunch_length/` in the output folder.

| Case | Run | Config (`configs/bunch_length/`) |
|---|---|---|
| reference, σ ≈ 0 | `20261030_143857` | `p0l0_a0_02_R25_sigma_ref.cfg` |
| q₁σ = 0.25 (σ = 0.56 nm) | `20261030_151536` | `p0l0_a0_02_R25_sigma_qs025.cfg` |
| q₁σ = 0.5 (σ = 1.1 nm) | `20261030_155222` | `p0l0_a0_02_R25_sigma_qs05.cfg` |
| q₁σ = 1 (σ = 2.2 nm) | `20261030_162911` | `p0l0_a0_02_R25_sigma_qs1.cfg` |

The Sec. 5.3 amplitude table is the output of `theory/future_dev/tests/bunch_length/compare_bunch_length.py` (no
arguments; it finds the runs here). The lens table comes from each run's `run_log.txt` (`lg_decomposition.py` with
a bounded `--z0` list). Full analysis: `theory/future_dev/tests/bunch_length/README.md`.

## Post-processing

Runs A–D: as in `config/expansion_a0_scan_server/README.md`, except that `lg_decomposition.py` was run with a wider
w0 list (the default range missed run D's N = 1 optimum):

```
W0=$(python3 -c "import numpy as np; print(','.join(f'{v:.4g}' for v in np.geomspace(1.5,75,101)))")
Z0=$(python3 -c "print(','.join(str(v) for v in range(-1500000, 1500001, 100000)))")
python3 py_scripts/radiation/lg_decomposition.py <run> --w0=$W0 --z0=$Z0 --z0-unit lambda
```

## Figures

`python3 make_figures.py` (in this folder) rebuilds `figures/fig1`–`fig7` from the PNGs of runs A–D and `fig9` from
those of the bunch-length runs (all under `~/output/progress_report_20261029/`).

`python3 make_lens_figure.py` builds Fig. 8 (Sec. 4.4) directly from `lens_gaussian_pulse/20261030_111559`.

## PDF

Needs pandoc and a LaTeX engine; `tectonic` (`conda create -p <env> -c conda-forge tectonic`) works with the
cluster's Anaconda pandoc 2.12:

```
pandoc progress_report.md -o progress_report.pdf --pdf-engine=tectonic -V papersize=a5 \
    -V geometry:margin=11mm -V fontsize=10pt -V colorlinks=true -H pdf_header.tex
```
