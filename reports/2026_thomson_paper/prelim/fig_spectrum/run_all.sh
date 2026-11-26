#!/bin/bash
# Build once on this (compute) node, run the 8 configs (cheap first), then post-process each run.
# Meant to be the whole Slurm job, so the node is released when it ends (theory/cluster.md):
#   cd ~/work/c++/CoherentThomson && git pull
#   tmux new -s spectra
#   srun --partition=cpu-wide --time=08:00:00 --cpus-per-task=96 --exclude=hpc-db01,hpc-ctrl01 \
#        bash reports/2026_thomson_paper/prelim/fig_spectrum/run_all.sh
#   (detach: Ctrl-b d; never type exit inside tmux)
# 96 CPUs, not --exclusive: num_threads 0 follows the allocation, and above ~130 threads the a0 = 1, N = 2 field
# buffers exceed the solver's 8 GB cap and the run is refused.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
cd "$REPO" || exit 1
echo "host $(hostname), $(nproc) CPUs, commit $(git rev-parse --short HEAD)"
cmake -B build/ && cmake --build build/ --clean-first -j "$(nproc)" || exit 1

IDS="10um_a0_0p1_N1 10um_a0_0p1_N2 800nm_a0_0p1_N1 800nm_a0_0p1_N2 10um_a0_1p0_N1 800nm_a0_1p0_N1 10um_a0_1p0_N2 800nm_a0_1p0_N2"
failed=""
for id in $IDS; do
  echo "=== $id  $(date)"
  ./bin/coherent_thomson_solver "$HERE/configs/$id.cfg" || failed="$failed $id"
done

PY="${PYTHON:-python3}"
for id in $IDS; do
  run="$(ls -d "$HOME/output/paper_2026/prelim/fig_spectrum/$id"/*/ 2>/dev/null | sort | tail -1)"
  [ -n "$run" ] && [ -f "$run/radiation_field.dat" ] || continue
  echo "=== plot $id: $run"
  "$PY" py_scripts/radiation/plot_spectrum_map.py "$run" || echo "plot failed: $id"
done
echo "done $(date); failed runs:${failed:- none}"
