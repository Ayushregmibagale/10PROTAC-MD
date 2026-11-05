#!/bin/bash
# Run this ON ERNIE after 05_verify_setup.sh passes.
# Launches all 13 pRosettaC runs in parallel (7 workers each = ~91 CPUs).
#
# IMPORTANT: run inside screen or tmux to survive SSH disconnect:
#   screen -S protac
#   bash ~/protac/scripts/ernie_setup/06_launch_parallel.sh
#   Ctrl-A D   (detach;  reconnect later: screen -r protac)
#
# Usage: bash ~/protac/scripts/ernie_setup/06_launch_parallel.sh

set -uo pipefail

PROTAC_DIR="$HOME/protac"
LOG_DIR="$PROTAC_DIR/ernie_logs"
SCRIPT="$PROTAC_DIR/scripts/rosetta_utils/run_local_prosettac.py"
ANALYSIS_SCRIPT="$PROTAC_DIR/scripts/analysis/run_poi_full_analysis.py"
AGGREGATE_SCRIPT="$PROTAC_DIR/scripts/analysis/aggregate_features.py"
N_WORKERS=7    # 13 POIs × 7 = 91 CPUs (leaves headroom for OS + Rosetta overhead)
N_DOCK=200

mkdir -p "$LOG_DIR"

echo "=== pRosettaC Parallel Launch on Ernie ==="
echo "POIs: 13   Workers each: $N_WORKERS   Approx CPUs used: $((13 * N_WORKERS))"
echo "Log directory: $LOG_DIR"
echo "Started: $(date)"
echo ""

# ── DDR2: clean Patchdock_Results from old DDR1 run ──────────────────────────
# DDR2 was previously run with DDR1/3ZOS structure; those results are invalid.
# Clean everything so the fresh run with 7AYM starts from scratch.
DDR2_DIR="$PROTAC_DIR/pRosettaC/runs/high_affinity_degraders/DDR2"
DDR2_PD="$DDR2_DIR/Patchdock_Results"
if [[ -d "$DDR2_PD" ]]; then
    n_pd=$(ls "$DDR2_PD" 2>/dev/null | wc -l)
    if [[ $n_pd -gt 0 ]]; then
        echo "[DDR2] Removing $n_pd stale PatchDock files from old DDR1 run..."
    fi
fi
echo "[DDR2] Cleaned — will run fresh with 7AYM structure."
echo ""

# ── Launch function ───────────────────────────────────────────────────────────
LAUNCHED_PIDS=()

launch() {
    local poi="$1"
    local category="$2"
    local tag="$3"     # log/PID filename tag (e.g. MAPK14_high, MAPK14_low)
    local run_dir="$PROTAC_DIR/pRosettaC/runs/$category/$poi"
    local log_file="$LOG_DIR/${tag}.log"
    local pid_file="$LOG_DIR/${tag}.pid"

    if [[ ! -d "$run_dir" ]]; then
        echo "[SKIP] $tag: directory not found ($run_dir)"
        return
    fi
    if [[ ! -f "$run_dir/Protac_params.txt" ]]; then
        echo "[SKIP] $tag: Protac_params.txt missing"
        return
    fi

    # Skip if already has Results (allow re-running incomplete runs)
    local n_res
    n_res=$(ls "$run_dir/Results/"*.pdb 2>/dev/null | wc -l)
    if [[ $n_res -gt 50 ]]; then
        echo "[SKIP] $tag: Results/ already has $n_res PDBs — skipping (delete Results/ to re-run)"
        return
    fi

    echo "[LAUNCH] $tag"
    (
        cd "$run_dir"
        conda run -n protac python "$SCRIPT" \
            Protac_params.txt \
            --n_dock  "$N_DOCK" \
            --n_workers "$N_WORKERS"
    ) > "$log_file" 2>&1 &

    local pid=$!
    echo "$pid" > "$pid_file"
    LAUNCHED_PIDS+=("$pid")
    echo "         PID $pid → $log_file"
}

# ── Launch all 13 POI runs ────────────────────────────────────────────────────
echo "--- high_affinity_degraders ---"
launch DDR2   high_affinity_degraders  DDR2
launch EPHB2  high_affinity_degraders  EPHB2
launch MAPK14 high_affinity_degraders  MAPK14_high
launch RIPK2  high_affinity_degraders  RIPK2

echo ""
echo "--- high_affinity_no_degradation ---"
launch ABL1   high_affinity_no_degradation  ABL1
launch AXL    high_affinity_no_degradation  AXL
launch EPHA2  high_affinity_no_degradation  EPHA2
launch MAP4K5 high_affinity_no_degradation  MAP4K5
launch SLK    high_affinity_no_degradation  SLK

echo ""
echo "--- low_affinity_degraders ---"
launch CDK17   low_affinity_degraders  CDK17
launch MAPK14  low_affinity_degraders  MAPK14_low
launch PIP4K2C low_affinity_degraders  PIP4K2C
launch TNIK    low_affinity_degraders  TNIK

echo ""
echo "All runs launched at $(date)."
echo ""
echo "Monitor live progress:"
echo "  tail -f $LOG_DIR/DDR2.log"
echo "  watch -n120 'for f in $LOG_DIR/*.log; do"
echo "    echo \"\$(basename \$f .log): \$(grep -c \"Generating cluster\" \$f 2>/dev/null || echo 0) clusters / \$(grep -c \"patchdock\" \$f 2>/dev/null || echo 0) docks\""
echo "  done'"
echo ""
echo "Results count per POI:"
echo "  for d in $PROTAC_DIR/pRosettaC/runs/*/*; do"
echo "    echo \"\$(basename \$d): \$(ls \$d/Results/*.pdb 2>/dev/null | wc -l) PDBs\""
echo "  done"
echo ""
echo "=== Waiting for all runs to complete (est. 13-14 hours) ==="

# ── Wait for all launched processes ──────────────────────────────────────────
FAILED=0
for pid in "${LAUNCHED_PIDS[@]}"; do
    if wait "$pid"; then
        : # success
    else
        echo "  WARNING: process $pid exited with non-zero status"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== All pRosettaC processes finished at $(date) ==="
if [[ $FAILED -gt 0 ]]; then
    echo "  $FAILED run(s) exited with errors — check logs before proceeding"
else
    echo "  All runs completed cleanly."
fi
echo ""

# ── Count results ─────────────────────────────────────────────────────────────
echo "=== Results summary ==="
for run_dir in "$PROTAC_DIR"/pRosettaC/runs/*/*; do
    poi=$(basename "$run_dir")
    cat=$(basename "$(dirname "$run_dir")")
    n=$(ls "$run_dir/Results/"*.pdb 2>/dev/null | wc -l)
    printf "  %-12s %-32s %d PDBs\n" "$poi" "$cat" "$n"
done
echo ""

if [[ $FAILED -gt 0 ]]; then
    echo "Fix errors above, then run analysis manually:"
    echo "  conda run -n protac python $ANALYSIS_SCRIPT --all"
    echo "  conda run -n protac python $AGGREGATE_SCRIPT"
    exit 1
fi

# ── Post-run analysis ─────────────────────────────────────────────────────────
echo "=== Running post-run analysis ==="
echo "[1/2] Per-POI cluster feature extraction..."
conda run -n protac python "$ANALYSIS_SCRIPT" --all 2>&1 | tee "$LOG_DIR/analysis.log"

echo ""
echo "[2/2] Cross-POI feature aggregation..."
conda run -n protac python "$AGGREGATE_SCRIPT" 2>&1 | tee -a "$LOG_DIR/analysis.log"

echo ""
echo "=== Pipeline complete at $(date) ==="
echo "  ML feature tables: $PROTAC_DIR/analysis/degradation_prediction/"
echo "  Per-POI reports:   each run directory (premd_report.txt, cluster_features.csv)"
echo ""
echo "Copy results back to your local machine:"
echo "  rsync -avz s1638308@ernie.bio.ed.ac.uk:~/protac/pRosettaC/runs/ \\"
echo "      \"$HOME/PROTAC degradation Prediction/pRosettaC/runs_ernie/\""
echo "  rsync -avz s1638308@ernie.bio.ed.ac.uk:~/protac/analysis/ \\"
echo "      \"$HOME/PROTAC degradation Prediction/analysis/\""
