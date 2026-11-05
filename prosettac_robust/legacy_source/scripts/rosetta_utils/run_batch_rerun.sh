#!/bin/bash
# Re-run Priority-1 (cleaned) and Priority-2 (new) POIs sequentially.
# Error in one POI is logged but doesn't stop the batch.
#
# Usage:
#   cd ~/PROTAC\ degradation\ Prediction
#   nohup bash scripts/rosetta_utils/run_batch_rerun.sh > pRosettaC/runs/batch_rerun.log 2>&1 &

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RUNNER="${PROJECT_ROOT}/scripts/rosetta_utils/run_local_prosettac.py"
ANALYSER="${PROJECT_ROOT}/scripts/analysis/run_poi_full_analysis.py"
RUNS="${PROJECT_ROOT}/pRosettaC/runs"

N_DOCK="${N_DOCK:-200}"
N_WORKERS="${N_WORKERS:-4}"

TARGETS=(
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

echo "=== pRosettaC batch re-run ==="
echo "N_DOCK=${N_DOCK}  N_WORKERS=${N_WORKERS}"
echo "Started: $(date)"
echo ""

for target_path in "${TARGETS[@]}"; do
  run_dir="${RUNS}/${target_path}"
  poi=$(basename "$target_path")
  echo ">>> $(date '+%Y-%m-%d %H:%M:%S')  ${target_path}"

  if [[ ! -f "${run_dir}/Protac_params.txt" ]]; then
    echo "  SKIP: Protac_params.txt not found at ${run_dir}"
    continue
  fi

  if [[ -d "${run_dir}/Results" && -n "$(ls ${run_dir}/Results/ 2>/dev/null)" ]]; then
    echo "  SKIP: Results/ already populated"
    continue
  fi

  # Skip if a run is already in progress (Patchdock_Results has files but no Results yet)
  if [[ -d "${run_dir}/Patchdock_Results" && -n "$(ls ${run_dir}/Patchdock_Results/ 2>/dev/null)" ]]; then
    echo "  SKIP: run already in progress (Patchdock_Results not empty)"
    continue
  fi

  (
    cd "${run_dir}" && \
    python "${RUNNER}" Protac_params.txt \
      --n_dock "${N_DOCK}" \
      --n_workers "${N_WORKERS}"
  )
  exit_code=$?
  if [[ $exit_code -ne 0 ]]; then
    echo "  WARNING: ${poi} exited with code ${exit_code} — continuing to next target"
  else
    echo "  Done: ${poi}"
    # Run per-POI analysis immediately after run completes
    if [[ -d "${run_dir}/Results" ]]; then
      echo "  Running analysis for ${poi}..."
      python "${ANALYSER}" "${poi}" 2>&1 | tail -10
    fi
  fi
  echo ""
done

echo "=== Batch complete: $(date) ==="
echo "Running final cross-POI aggregation..."
python "${PROJECT_ROOT}/scripts/analysis/aggregate_features.py"
