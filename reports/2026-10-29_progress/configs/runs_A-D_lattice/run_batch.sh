#!/bin/bash
# Progress report Sec. 5: runs A-D on the polar lattice (plus D with the spacing halved), then the same
# post-processing as runs A-D (theory/scientific_report/README.md, Post-processing). Run inside a compute-node job
# (theory/srun.txt), e.g. from the head node:
#   sbatch --partition=cpu-wide --time=02:00:00 --exclusive --exclude=hpc-ctrl01,hpc-db01 \
#          --output=$HOME/output/progress_report_20261029/logs/runs_A-D_lattice/slurm-%j.out \
#          theory/scientific_report/configs/runs_A-D_lattice/run_batch.sh
# Builds first on this node (-march=native). Runs go to ~/output/progress_report_20261029/runs_A-D_lattice/ (the
# configs' output_folder); progress and the config -> run folder mapping to logs/runs_A-D_lattice/index.txt there.
set -u
REPO=~/work/c++/CoherentThomson
CFG=theory/scientific_report/configs/runs_A-D_lattice
OUT=~/output/progress_report_20261029/runs_A-D_lattice
LOGS=~/output/progress_report_20261029/logs/runs_A-D_lattice
mkdir -p "$OUT" "$LOGS"
cd "$REPO" || exit 1

{ cmake -B build/ && cmake --build build/ --clean-first; } > "$LOGS/build.log" 2>&1 || { echo "build failed, see $LOGS/build.log"; exit 1; }
grep "Build type" "$LOGS/build.log"

W0=$(python3 -c "import numpy as np; print(','.join(f'{v:.4g}' for v in np.geomspace(1.5,75,101)))")
Z0=$(python3 -c "print(','.join(str(v) for v in range(-1500000, 1500001, 100000)))")

for name in p0l0_a0_01_R15 p2l2_a0_01_R15 p0l0_a0_02_R25 p2l2_a0_02_R25 p2l2_a0_02_R25_fine; do
  f=$CFG/config_expansion_${name%_fine}_lattice${name#${name%_fine}}.cfg
  before=$(ls -d "$OUT"/2*/ 2>/dev/null | sort | tail -1)
  echo "$(date '+%F %T') START $name on $(hostname)" | tee -a "$LOGS/index.txt"
  ./bin/coherent_thomson_solver "$f" > "$LOGS/$name.log" 2>&1
  rc=$?
  run=$(ls -d "$OUT"/2*/ 2>/dev/null | sort | tail -1)
  [ "$run" = "$before" ] && { echo "$(date '+%F %T') END   $name rc=$rc run=(no new run folder)" | tee -a "$LOGS/index.txt"; continue; }
  echo "$(date '+%F %T') END   $name rc=$rc run=$run" | tee -a "$LOGS/index.txt"
  {
    python3 py_scripts/radiation/plot_all_components.py "$run"
    python3 py_scripts/radiation/plot_observables.py "$run"
    python3 py_scripts/radiation/lg_decomposition.py "$run" --w0=$W0 --z0=$Z0 --z0-unit lambda
    python3 py_scripts/radiation/plot_processed.py radial_phase total 0 1 "$run"
  } > "$LOGS/post_$name.log" 2>&1
  echo "$(date '+%F %T') POST  $name rc=$?" | tee -a "$LOGS/index.txt"
  sleep 2   # run folders are timestamped to the second
done
echo "$(date '+%F %T') ALL DONE" | tee -a "$LOGS/index.txt"
