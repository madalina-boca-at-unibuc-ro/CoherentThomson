#!/bin/bash
# Runs the four bunch-length configs (README.md here) one after another, each on all CPUs of the allocation
# (num_threads 0), ~37 min each on 256 CPUs. Start it inside the compute-node job (theory/srun.txt: srun inside tmux
# on the head node), from anywhere:   bash theory/future_dev/tests/bunch_length/run_batch.sh
# Builds first on this node (-march=native: the binary must be built where it runs). Progress and the
# config -> run folder mapping: ~/output_tests/bunch_length/logs/index.txt; solver output: one .log per config there.
set -u
REPO="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT=~/output_tests/bunch_length
LOGS=$OUT/logs
mkdir -p "$LOGS"
cd "$REPO" || exit 1

{ cmake -B build/ && cmake --build build/ --clean-first; } > "$LOGS/build.log" 2>&1 || { echo "build failed, see $LOGS/build.log"; exit 1; }
grep "Build type" "$LOGS/build.log"

for name in ref qs025 qs05 qs1; do
  f=theory/future_dev/tests/bunch_length/p0l0_a0_02_R25_sigma_$name.cfg
  before=$(ls -d "$OUT"/2*/ 2>/dev/null | sort | tail -1)
  echo "$(date '+%F %T') START $name on $(hostname)" | tee -a "$LOGS/index.txt"
  ./bin/coherent_thomson_solver "$f" > "$LOGS/$name.log" 2>&1
  rc=$?
  run=$(ls -d "$OUT"/2*/ 2>/dev/null | sort | tail -1)
  [ "$run" = "$before" ] && run="(no new run folder)"
  echo "$(date '+%F %T') END   $name rc=$rc run=$run" | tee -a "$LOGS/index.txt"
  sleep 2   # run folders are timestamped to the second
done
echo "$(date '+%F %T') ALL DONE" | tee -a "$LOGS/index.txt"
