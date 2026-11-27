#!/bin/bash
# Sec. IV.D/E (LG content, ponderomotive lens), stage 1: build once on this (compute) node, run the given config ids
# (default: the a0 = 0.25 lines first, since they set the a0 = 0.25 map frequencies, then the a0 = 0.1 lines and the
# a0 = 0.1 maps), then post-process each run: line runs with plot_spectrum_map.py, map runs as fig_maps (scaled
# parameters into run_log.txt, observables, carrier-removed field and helicity components, LG decomposition).
# Meant to be the whole Slurm job, so the node is released when it ends (theory/cluster.md):
#   cd ~/work/c++/CoherentThomson && git pull
#   tmux new -s lg
#   srun --partition=cpu-wide --time=05:00:00 --cpus-per-task=128 --exclude=hpc-db01,hpc-ctrl01 \
#        bash reports/2026_thomson_paper/prelim/fig_lg/run_all.sh [id ...]
#   (detach: Ctrl-b d; never type exit inside tmux)
# Use 128 CPUs, not 256: the a0 = 0.25, N = 2 line buffers (~0.04 GB per thread) would pass the 8 GB cap.
# Estimate (make_configs.py): ~250 CPU h, ~2 h on 128 CPUs.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
cd "$REPO" || exit 1
if ! command -v cmake >/dev/null && [ -d "$HOME/anaconda3/bin" ]; then export PATH="$HOME/anaconda3/bin:$PATH"; fi
echo "host $(hostname), $(nproc) CPUs, commit $(git rev-parse --short HEAD)"
cmake -B build/ && cmake --build build/ --clean-first -j "$(nproc)" || exit 1

DEFAULT=""
for a in 0p25 0p1; do
  for m in LG02 LG22; do for N in 1 2; do DEFAULT="$DEFAULT 10um_${m}_a0_${a}_N${N}_line"; done; done
done
DEFAULT="$DEFAULT 10um_LG02_a0_0p1_map 10um_LG22_a0_0p1_map"
IDS="${*:-$DEFAULT}"
failed=""
for id in $IDS; do
  echo "=== $id  $(date)"
  ./bin/coherent_thomson_solver "$HERE/configs/$id.cfg" || failed="$failed $id"
done

PY="${PYTHON:-python3}"
for id in $IDS; do
  run="$(ls -d "$HOME/output/paper_2026/prelim/fig_lg/$id"/*/ 2>/dev/null | sort | tail -1)"
  [ -n "$run" ] && [ -f "$run/radiation_field.dat" ] || continue
  echo "=== post-process $id: $run  $(date)"
  "$PY" py_scripts/utils/scaled_parameters.py "$run" --write-log || echo "scaled_parameters failed: $id"
  case "$id" in
    *_line)
      "$PY" py_scripts/radiation/plot_spectrum_map.py "$run" || echo "plot_spectrum_map failed: $id" ;;
    *)
      "$PY" py_scripts/radiation/plot_observables.py "$run" || echo "plot_observables failed: $id"
      for what in total plus minus; do
        "$PY" py_scripts/radiation/plot_field_carrier_removed.py "$what" "$run" || echo "carrier_removed $what failed: $id"
      done
      "$PY" py_scripts/radiation/lg_decomposition.py "$run" || echo "lg_decomposition failed: $id" ;;
  esac
done
echo "done $(date); failed runs:${failed:- none}"
