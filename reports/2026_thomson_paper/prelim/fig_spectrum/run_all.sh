#!/bin/bash
# Build once on this (compute) node, run the given config ids (default: the a0 = 1 N_c = 8 and a0 = 0.5 N_c = 30
# batch, cheap first), then post-process each run.
# Meant to be the whole Slurm job, so the node is released when it ends (theory/cluster.md):
#   cd ~/work/c++/CoherentThomson && git pull
#   tmux new -s spectra
#   srun --partition=cpu-wide --time=08:00:00 --cpus-per-task=96 --exclude=hpc-db01,hpc-ctrl01 \
#        bash reports/2026_thomson_paper/prelim/fig_spectrum/run_all.sh [id ...]
#   (detach: Ctrl-b d; never type exit inside tmux)
# Not --exclusive: num_threads 0 follows the allocation, and with 256 threads the N = 2 field buffers (~0.04-0.06 GB
# per thread) exceed the solver's 8 GB cap and the run is refused. 96-128 CPUs.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
cd "$REPO" || exit 1
# On the cluster cmake and the plotting Python come from Anaconda, on PATH only in interactive shells
# (theory/cluster.md, "PATH in non-interactive shells").
if ! command -v cmake >/dev/null && [ -d "$HOME/anaconda3/bin" ]; then export PATH="$HOME/anaconda3/bin:$PATH"; fi
echo "host $(hostname), $(nproc) CPUs, commit $(git rev-parse --short HEAD)"
cmake -B build/ && cmake --build build/ --clean-first -j "$(nproc)" || exit 1

# First batch (2026-11-26, job 32515): 10um/800nm_a0_0p1_N{1,2} and _a0_1p0_N{1,2} (N_c = 30).
IDS="${*:-10um_a0_1p0_Nc8_N1 800nm_a0_1p0_Nc8_N1 10um_a0_0p5_Nc30_N1 800nm_a0_0p5_Nc30_N1 10um_a0_1p0_Nc8_N2 800nm_a0_1p0_Nc8_N2 10um_a0_0p5_Nc30_N2 800nm_a0_0p5_Nc30_N2}"
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
