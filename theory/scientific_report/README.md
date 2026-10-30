# Scientific report (progress report, 2026-10-29)

`progress_report.md` is the report for the group; `progress_report.pdf` is the version to send (A5 pages, readable on
a phone). Notes for rebuilding it:

## Runs

| Report label | Case | Run folder (`~/output/`) |
|---|---|---|
| A | LG(0,0), a0 = 0.1, R_max 15 w0 | `20261029_132340` |
| B | LG(0,0), a0 = 0.2, R_max 25 w0 | `20261029_143739` |
| C | LG(2,2), a0 = 0.1, R_max 15 w0 | `20261029_140030` |
| D | LG(2,2), a0 = 0.2, R_max 25 w0 | `20261029_151457` |

Configs: `config/expansion_a0_scan_server/`. Post-processing as in that folder's `README.md`, except that
`lg_decomposition.py` was run with a wider w0 list (the default range missed run D's N = 1 optimum):

```
W0=$(python3 -c "import numpy as np; print(','.join(f'{v:.4g}' for v in np.geomspace(1.5,75,101)))")
Z0=$(python3 -c "print(','.join(str(v) for v in range(-1500000, 1500001, 100000)))")
python3 py_scripts/radiation/lg_decomposition.py <run> --w0=$W0 --z0=$Z0 --z0-unit lambda
```

## Figures

`python3 make_figures.py` (in this folder) rebuilds `figures/*.jpg` from the four run folders' PNGs. It needs the
run folders under `~/output/`.

`python3 make_lens_figure.py` builds Fig. 8 (Sec. 4.4) directly from the lattice run at `a0 = 0.1`,
`~/output_tests/ponderomotive_lens/20261030_111559` (pure Gaussian pulse). That run and the others quoted in
Sec. 1 and 4.4 are listed in `theory/future_dev/ponderomotive_lens.md` (Results); their configs are in
`theory/future_dev/tests/ponderomotive_lens/`.

## PDF

Needs pandoc and a LaTeX engine; `tectonic` (`conda create -p <env> -c conda-forge tectonic`) works with the
cluster's Anaconda pandoc 2.12:

```
pandoc progress_report.md -o progress_report.pdf --pdf-engine=tectonic -V papersize=a5 \
    -V geometry:margin=11mm -V fontsize=10pt -V colorlinks=true -H pdf_header.tex
```
