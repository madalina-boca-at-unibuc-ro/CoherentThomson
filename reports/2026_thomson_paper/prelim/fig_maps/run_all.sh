#!/bin/bash
# Sec. IV.C field maps: build once on this (compute) node, run the given config ids (default: all four, N = 1 first),
# then post-process each run (scaled parameters into run_log.txt, observable maps, carrier-removed total field and
# helicity components, LG decomposition).
# Meant to be the whole Slurm job, so the node is released when it ends (theory/cluster.md):
#   cd ~/work/c++/CoherentThomson && git pull
#   tmux new -s maps
#   srun --partition=cpu-wide --time=04:00:00 --cpus-per-task=128 --exclude=hpc-db01,hpc-ctrl01 \
#        bash reports/2026_thomson_paper/prelim/fig_maps/run_all.sh [id ...]
#   (detach: Ctrl-b d; never type exit inside tmux)
# Field buffers ~19 MB per thread (66 560 points, one frequency): any CPU count is fine for memory.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
cd "$REPO" || exit 1
# On the cluster cmake and the plotting Python come from Anaconda, on PATH only in interactive shells
# (theory/cluster.md, "PATH in non-interactive shells").
if ! command -v cmake >/dev/null && [ -d "$HOME/anaconda3/bin" ]; then export PATH="$HOME/anaconda3/bin:$PATH"; fi
echo "host $(hostname), $(nproc) CPUs, commit $(git rev-parse --short HEAD)"
cmake -B build/ && cmake --build build/ --clean-first -j "$(nproc)" || exit 1

IDS="${*:-10um_a0_1p0_Nc8_N1 800nm_a0_1p0_Nc8_N1 10um_a0_1p0_Nc8_N2 800nm_a0_1p0_Nc8_N2}"
failed=""
for id in $IDS; do
  echo "=== $id  $(date)"
  ./bin/coherent_thomson_solver "$HERE/configs/$id.cfg" || failed="$failed $id"
done

PY="${PYTHON:-python3}"
for id in $IDS; do
  run="$(ls -d "$HOME/output/paper_2026/prelim/fig_maps/$id"/*/ 2>/dev/null | sort | tail -1)"
  [ -n "$run" ] && [ -f "$run/radiation_field.dat" ] || continue
  echo "=== post-process $id: $run  $(date)"
  "$PY" py_scripts/utils/scaled_parameters.py "$run" --write-log || echo "scaled_parameters failed: $id"
  "$PY" py_scripts/radiation/plot_observables.py "$run" || echo "plot_observables failed: $id"
  for what in total plus minus; do
    "$PY" py_scripts/radiation/plot_field_carrier_removed.py "$what" "$run" || echo "carrier_removed $what failed: $id"
  done
  "$PY" py_scripts/radiation/lg_decomposition.py "$run" || echo "lg_decomposition failed: $id"
done
echo "done $(date); failed runs:${failed:- none}"
