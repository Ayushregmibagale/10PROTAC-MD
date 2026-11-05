#!/bin/bash
# run_queue_ernie.sh
# Waits for DDR2 (PID given as first argument) to finish, then runs
# RIPK2 → EPHB2 → MAPK14 sequentially on ernie.
#
# Usage:
#   bash ~/protac/scripts/rosetta_utils/run_queue_ernie.sh <DDR2_PID> 2>&1 | tee ~/queue.log &
#   echo "Queue PID: $!"

DDR2_PID=${1:?Usage: $0 <DDR2_PID>}
RUNS_DIR=~/protac/pRosettaC/runs/high_affinity_degraders
PYCACHE=~/protac/scripts/rosetta_utils/__pycache__
POIS=(RIPK2 EPHB2 MAPK14)

if [[ "$DDR2_PID" -eq 0 ]]; then
    echo "[queue] DDR2 already done — skipping wait."
else
    echo "[queue] Waiting for DDR2 (PID $DDR2_PID) to finish..."
    while kill -0 "$DDR2_PID" 2>/dev/null; do
        sleep 60
    done
    echo "[queue] DDR2 finished at $(date)"
fi
echo ""

for POI in "${POIS[@]}"; do
    echo "============================================"
    echo "[queue] Starting $POI at $(date)"
    echo "============================================"

    POI_DIR="$RUNS_DIR/$POI"
    PR="$POI_DIR/Patchdock_Results"

    if [[ ! -d "$PR" ]]; then
        echo "[queue] ERROR: $PR not found — skipping $POI"
        continue
    fi

    # Clean stale outputs from any previous (broken) run
    echo "[queue] Cleaning stale confs and combined PDBs for $POI..."
    find "$PR" -name "confs_*.sdf" -delete
    find "$PR" -name "combined_*_0001.pdb" -delete
    find "$PYCACHE" -name "pyrosetta_adapter*.pyc" -delete 2>/dev/null
    echo "[queue] Clean done."

    cd "$POI_DIR" || { echo "[queue] ERROR: cannot cd to $POI_DIR"; continue; }

    N_WORKERS_RESUME=16 python ~/protac/scripts/rosetta_utils/resume_pipeline.py \
        2>&1 | tee ~/poi_run_${POI}.log

    echo "[queue] $POI finished at $(date)"
    echo ""
done

echo "[queue] All queued POIs complete at $(date)."
echo "[queue] Logs: ~/poi_run_RIPK2.log, ~/poi_run_EPHB2.log, ~/poi_run_MAPK14.log"
