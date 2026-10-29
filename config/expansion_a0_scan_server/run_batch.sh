#!/bin/bash
# Runs the four server configs of the a0 scan one after another, each on all CPUs of the allocation
# (num_threads 0). Start it inside the compute-node job (theory/srun.txt: srun inside tmux on the head node), from
# anywhere:   bash config/expansion_a0_scan_server/run_batch.sh
# Builds first on this node (-march=native: the binary must be built where it runs). Progress and the
# config -> run folder mapping: ~/output/a0_scan_server_logs/index.txt; solver output: one .log per config there.
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
LOGS=~/output/a0_scan_server_logs
mkdir -p "$LOGS"
cd "$REPO" || exit 1

{ cmake -B build/ && cmake --build build/ --clean-first; } > "$LOGS/build.log" 2>&1 || { echo "build failed, see $LOGS/build.log"; exit 1; }
grep "Build type" "$LOGS/build.log"

for f in config/expansion_a0_scan_server/config_expansion_p0l0_a0_01_R15.cfg \
         config/expansion_a0_scan_server/config_expansion_p2l2_a0_01_R15.cfg \
         config/expansion_a0_scan_server/config_expansion_p0l0_a0_02_R25.cfg \
         config/expansion_a0_scan_server/config_expansion_p2l2_a0_02_R25.cfg; do
  name=$(basename "$f" .cfg)
  before=$(ls -d ~/output/2*/ 2>/dev/null | sort | tail -1)
  echo "$(date '+%F %T') START $name on $(hostname)" | tee -a "$LOGS/index.txt"
  ./bin/coherent_thomson_solver "$f" > "$LOGS/$name.log" 2>&1
  rc=$?
  run=$(ls -d ~/output/2*/ 2>/dev/null | sort | tail -1)
  [ "$run" = "$before" ] && run="(no new run folder)"
  echo "$(date '+%F %T') END   $name rc=$rc run=$run" | tee -a "$LOGS/index.txt"
  sleep 2   # run folders are timestamped to the second
done
echo "$(date '+%F %T') ALL DONE" | tee -a "$LOGS/index.txt"
