#!/bin/bash
# Run pRosettaC for a single POI.
# Usage: ./run_prosettac.sh <TARGET> <CATEGORY>
#
# Example:
#   ./run_prosettac.sh MET high_affinity_degraders
#   ./run_prosettac.sh AXL high_affinity_no_degradation

set -euo pipefail

TARGET="${1:?Usage: $0 <TARGET> <CATEGORY>}"
CATEGORY="${2:?Usage: $0 <TARGET> <CATEGORY>}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

CONFIG="${PROJECT_ROOT}/pRosettaC/configs/${CATEGORY}/${TARGET}.flags"
LOG="${PROJECT_ROOT}/pRosettaC/results/logs/${TARGET}.log"
MODELS_DIR="${PROJECT_ROOT}/pRosettaC/results/models/${TARGET}"

mkdir -p "${MODELS_DIR}"

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: Config not found: ${CONFIG}"
    echo "Copy the template and fill in POI-specific fields:"
    echo "  cp ${PROJECT_ROOT}/pRosettaC/configs/template.flags ${CONFIG}"
    exit 1
fi

# ── Set this to your pRosettaC installation ──────────────────────────────────
PROSETTAC_SCRIPT="${PROSETTAC_HOME:-/path/to/pRosettaC}/pRosettaC.py"
ROSETTA_BIN="${ROSETTA_HOME:-/path/to/rosetta}/main/source/bin"
# ─────────────────────────────────────────────────────────────────────────────

echo "========================================"
echo " pRosettaC run: ${TARGET} (${CATEGORY})"
echo " Config : ${CONFIG}"
echo " Log    : ${LOG}"
echo " Models : ${MODELS_DIR}"
echo "========================================"

python "${PROSETTAC_SCRIPT}" \
    @"${CONFIG}" \
    -rosetta_bin "${ROSETTA_BIN}" \
    2>&1 | tee "${LOG}"

echo ""
echo "Run complete. Check scores at:"
echo "  ${PROJECT_ROOT}/pRosettaC/results/scores/${TARGET}_scores.sc"
