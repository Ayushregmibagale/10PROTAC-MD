#!/bin/bash
# Run pRosettaC locally for all 14 targets (sequential, one at a time).
# Adjust N_DOCK and N_WORKERS for your hardware.
#
# Usage:
#   cd ~/PROTAC\ degradation\ Prediction
#   bash scripts/rosetta_utils/batch_run_local.sh
#
# Each target creates its output under:
#   pRosettaC/runs/<category>/<target>/Results/
#
# Logs: pRosettaC/runs/<category>/<target>/log.txt (written by pRosettaC)

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RUNNER="${PROJECT_ROOT}/scripts/rosetta_utils/run_local_prosettac.py"
RUNS="${PROJECT_ROOT}/pRosettaC/runs"

N_DOCK="${N_DOCK:-200}"        # global docking poses per target (200 = practical default)
N_WORKERS="${N_WORKERS:-4}"    # parallel local-docking workers

# Priority order: degraders first, then non-degraders
TARGETS=(
  "high_affinity_degraders/MET"
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
  "low_affinity_degraders/MAPK14"
  "low_affinity_degraders/PIP4K2C"
  "low_affinity_degraders/CDK17"
)

echo "=== pRosettaC local batch run ==="
echo "N_DOCK=${N_DOCK}  N_WORKERS=${N_WORKERS}"
echo ""

for target_path in "${TARGETS[@]}"; do
  run_dir="${RUNS}/${target_path}"
  echo ">>> $(date '+%H:%M:%S')  ${target_path}"

  if [[ ! -f "${run_dir}/Protac_params.txt" ]]; then
    echo "  SKIP: Protac_params.txt not found"
    continue
  fi

  if [[ -d "${run_dir}/Results" ]]; then
    echo "  SKIP: Results/ already exists (delete to re-run)"
    continue
  fi

  (
    cd "${run_dir}"
    python "${RUNNER}" Protac_params.txt \
      --n_dock "${N_DOCK}" \
      --n_workers "${N_WORKERS}"
  )
  echo "  Done."
done

echo ""
echo "=== Batch complete ==="
