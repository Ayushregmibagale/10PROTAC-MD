#!/bin/bash
# Correct per-POI status check for pRosettaC runs.
#
# `ls Results/*.pdb | wc -l` ALWAYS returns 0, even for finished runs --
# clustering.py writes Results/tmpN/ subdirectories (one per cluster), not flat
# .pdb files at the top level. The real completion marker is RUN_DIR/result_summary.txt,
# written one level above Results/. See the "Monitoring pitfall" note in
# docs/05_bugs_and_fixes.md.
#
# This script also distinguishes a live-but-stuck process (e.g. DDR2, Bug 10)
# from a dead/crashed one (e.g. CDK17, Bug 9) by checking the PID file written
# by 06_launch_parallel.sh.
#
# Usage:
#   bash check_status.sh [PROTAC_DIR]
#
# PROTAC_DIR defaults to ~/protac (the symlink used both on ernie and locally).

set -uo pipefail

PROTAC_DIR="${1:-$HOME/protac}"
RUNS_DIR="$PROTAC_DIR/pRosettaC/runs"
LOG_DIR="$PROTAC_DIR/ernie_logs"

# MAPK14 appears under two categories with separate launch tags
# (high_affinity_degraders/MAPK14 and low_affinity_degraders/MAPK14 -- the
# latter is a duplicate, see Bug 12). Map (poi, category) -> launch tag so we
# can find the right .pid file.
get_tag() {
    local poi="$1" category="$2"
    if [[ "$poi" == "MAPK14" ]]; then
        case "$category" in
            high_affinity_degraders) echo "MAPK14_high" ;;
            low_affinity_degraders)  echo "MAPK14_low"  ;;
            *) echo "$poi" ;;
        esac
    else
        echo "$poi"
    fi
}

is_alive() {
    local tag="$1"
    local pidfile="$LOG_DIR/${tag}.pid"
    [[ -f "$pidfile" ]] || return 1
    local pid
    pid=$(cat "$pidfile")
    kill -0 "$pid" 2>/dev/null
}

printf "%-10s %-32s %-12s %s\n" "POI" "CATEGORY" "STATUS" "DETAIL"
printf "%-10s %-32s %-12s %s\n" "---" "--------" "------" "------"

for run_dir in "$RUNS_DIR"/*/*; do
    [[ -d "$run_dir" ]] || continue
    poi=$(basename "$run_dir")
    category=$(basename "$(dirname "$run_dir")")
    pd="$run_dir/Patchdock_Results"
    summary="$run_dir/result_summary.txt"

    if [[ -f "$summary" ]]; then
        mtime=$(date -r "$summary" '+%Y-%m-%d %H:%M')
        printf "%-10s %-32s %-12s %s\n" "$poi" "$category" "DONE" "result_summary.txt @ $mtime"
        continue
    fi

    tag=$(get_tag "$poi" "$category")

    if [[ -d "$pd" ]]; then
        n_docked=$(ls "$pd"/pd.*_docking_????.pdb 2>/dev/null | wc -l)
        n_combined=$(ls "$pd"/combined_*_0001.pdb 2>/dev/null | wc -l)

        if [[ $n_docked -gt 0 && $n_combined -gt 0 ]]; then
            pct=$((100 * n_combined / n_docked))
            status="RUNNING"
            is_alive "$tag" || status="STALLED"
            printf "%-10s %-32s %-12s %s\n" "$poi" "$category" "$status" "constraint_gen: $n_combined/$n_docked combined (~${pct}%)"
        elif [[ $n_docked -gt 0 ]]; then
            status="RUNNING"
            is_alive "$tag" || status="CRASHED"
            printf "%-10s %-32s %-12s %s\n" "$poi" "$category" "$status" "docking done ($n_docked PDBs); constraint_gen: 0 combined"
        else
            status="RUNNING"
            is_alive "$tag" || status="CRASHED"
            printf "%-10s %-32s %-12s %s\n" "$poi" "$category" "$status" "docking in progress (0 docked PDBs yet)"
        fi
        continue
    fi

    printf "%-10s %-32s %-12s %s\n" "$poi" "$category" "NOT-STARTED" "no Patchdock_Results/"
done
