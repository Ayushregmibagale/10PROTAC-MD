#!/bin/bash
# Poll every 5 minutes; run analysis for any POI that has new Results/ but no cluster_features.csv.
# Launch this once in background alongside the main batch.
#
# Usage:
#   nohup bash watch_and_analyse.sh > /tmp/watch_analyse.log 2>&1 &

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ANALYSER="${PROJECT_ROOT}/scripts/analysis/run_poi_full_analysis.py"
RUNS="${PROJECT_ROOT}/pRosettaC/runs"

ALL_TARGETS=(
  "high_affinity_degraders/DDR2"
  "high_affinity_degraders/RIPK2"
  "high_affinity_degraders/EPHB2"
  "high_affinity_degraders/MAPK14"
  "high_affinity_no_degradation/AXL"
  "high_affinity_no_degradation/SLK"
  "high_affinity_no_degradation/ABL1"
  "high_affinity_no_degradation/EPHA2"
  "high_affinity_no_degradation/MAP4K5"
  "low_affinity_degraders/TNIK"
  "low_affinity_degraders/PIP4K2C"
  "low_affinity_degraders/CDK17"
)

echo "=== Watcher started $(date) ==="

while true; do
  for target_path in "${ALL_TARGETS[@]}"; do
    poi=$(basename "$target_path")
    run_dir="${RUNS}/${target_path}"
    results="${run_dir}/Results"
    features="${run_dir}/cluster_features.csv"

    # Results/ exists with clusters, but analysis hasn't been run yet
    if [[ -d "$results" ]] && ls "$results"/cluster* &>/dev/null && [[ ! -f "$features" ]]; then
      echo "$(date '+%H:%M:%S') Analysing ${poi}..."
      python "${ANALYSER}" "${poi}" 2>&1
      echo "$(date '+%H:%M:%S') Done ${poi}"
    fi
  done

  # Check if all done → run aggregation and exit
  n_done=0
  for target_path in "${ALL_TARGETS[@]}"; do
    poi=$(basename "$target_path")
    [[ -f "${RUNS}/${target_path}/cluster_features.csv" ]] && ((n_done++))
  done
  echo "$(date '+%H:%M:%S') ${n_done}/${#ALL_TARGETS[@]} POIs analysed"

  if [[ $n_done -eq ${#ALL_TARGETS[@]} ]]; then
    echo "All POIs done — running aggregation"
    python "${PROJECT_ROOT}/scripts/analysis/aggregate_features.py"
    echo "=== Watcher done $(date) ==="
    exit 0
  fi

  sleep 300  # check every 5 minutes
done
