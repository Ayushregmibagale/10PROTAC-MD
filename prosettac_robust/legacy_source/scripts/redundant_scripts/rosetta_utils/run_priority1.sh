#!/bin/bash
# Run all Priority-1 (High Affinity + Strong Degrader) POIs sequentially.
# MET is already done; this script runs DDR2, RIPK2, EPHB2, MAPK14.
#
# Usage:
#   bash run_priority1.sh [--n_dock N]   (default N=50)
#
# Log files: <RUN_DIR>/pipeline_run.log per POI

set -euo pipefail

N_DOCK="${N_DOCK:-50}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUN_SCRIPT="${SCRIPT_DIR}/run_local_prosettac.py"

RUNS_BASE="/home/ayushregmibagale/PROTAC degradation Prediction/pRosettaC/runs/high_affinity_degraders"

TARGETS=(DDR2 RIPK2 EPHB2 MAPK14)

for TARGET in "${TARGETS[@]}"; do
    RUN_DIR="${RUNS_BASE}/${TARGET}"
    LOG="${RUN_DIR}/pipeline_run.log"

    if [[ -d "${RUN_DIR}/Results" ]]; then
        echo "=========================================="
        echo "  SKIP ${TARGET} — Results/ already exists"
        echo "=========================================="
        continue
    fi

    echo ""
    echo "=========================================="
    echo "  Starting ${TARGET}  (n_dock=${N_DOCK})"
    echo "  $(date)"
    echo "=========================================="

    cd "${RUN_DIR}"
    python "${RUN_SCRIPT}" Protac_params.txt --n_dock "${N_DOCK}" 2>&1 | tee "${LOG}"

    echo ""
    echo "  ${TARGET} done at $(date)"
done

echo ""
echo "All Priority-1 POI runs complete."
