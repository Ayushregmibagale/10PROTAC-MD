#!/bin/bash
# Submit all pRosettaC runs in priority order.
# Priority 1 and 2 run first (most informative targets).
# Adjust PARALLEL_JOBS for your cluster/workstation.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUN_SCRIPT="${SCRIPT_DIR}/run_prosettac.sh"
PARALLEL_JOBS="${PARALLEL_JOBS:-4}"

run_group() {
    local label="$1"; shift
    echo ""
    echo "══════════════════════════════════════════"
    echo "  ${label}"
    echo "══════════════════════════════════════════"
    # Run jobs from stdin in parallel
    printf '%s\n' "$@" | xargs -P "${PARALLEL_JOBS}" -I{} bash -c "{}"
}

# Priority 1 — High affinity + Strong Degraders
run_group "PRIORITY 1: High Affinity + Strong Degraders" \
    "bash ${RUN_SCRIPT} MET     high_affinity_degraders" \
    "bash ${RUN_SCRIPT} DDR2    high_affinity_degraders" \
    "bash ${RUN_SCRIPT} RIPK2   high_affinity_degraders" \
    "bash ${RUN_SCRIPT} EPHB2   high_affinity_degraders" \
    "bash ${RUN_SCRIPT} MAPK14  high_affinity_degraders"

# Priority 2 — High affinity + No Degradation  (most informative for model)
run_group "PRIORITY 2: High Affinity + No Degradation" \
    "bash ${RUN_SCRIPT} AXL    high_affinity_no_degradation" \
    "bash ${RUN_SCRIPT} SLK    high_affinity_no_degradation" \
    "bash ${RUN_SCRIPT} ABL1   high_affinity_no_degradation" \
    "bash ${RUN_SCRIPT} EPHA2  high_affinity_no_degradation" \
    "bash ${RUN_SCRIPT} MAP4K5 high_affinity_no_degradation"

# Priority 3 — Low affinity + Degradation
run_group "PRIORITY 3: Low Affinity + Degradation" \
    "bash ${RUN_SCRIPT} TNIK    low_affinity_degraders" \
    "bash ${RUN_SCRIPT} MAPK14  low_affinity_degraders" \
    "bash ${RUN_SCRIPT} PIP4K2C low_affinity_degraders" \
    "bash ${RUN_SCRIPT} CDK17   low_affinity_degraders"

echo ""
echo "All pRosettaC runs submitted."
