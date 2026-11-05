#!/bin/bash
# Run ON ERNIE to re-run pRosettaC for 4 artifact POIs with chain-D kinase fix.
# Run inside screen/tmux to survive SSH disconnect:
#   screen -S chainD
#   bash ~/protac/scripts/ernie_setup/08_ernie_launch_chain_d.sh
#   Ctrl-A D

set -uo pipefail

PROTAC_DIR="$HOME/protac"
LOG_DIR="$PROTAC_DIR/ernie_logs/chain_d_reruns"
SCRIPT="$PROTAC_DIR/scripts/rosetta_utils/run_local_prosettac.py"
N_WORKERS=7    # 4 POIs × 7 = 28 CPUs total
N_DOCK=200

mkdir -p "$LOG_DIR"

echo "=== Chain-D pRosettaC Re-runs — $(date) ==="
echo "POIs: TNIK (low_affinity_degraders)"
echo "      MAPK14 (low_affinity_degraders)"
echo "      MAP4K5 (high_affinity_no_degradation)"
echo "      SLK (high_affinity_no_degradation)"
echo ""

launch() {
    local poi="$1"
    local category="$2"
    local run_dir="$PROTAC_DIR/pRosettaC/runs/$category/$poi"
    local log_file="$LOG_DIR/${poi}.log"

    if [[ ! -f "$run_dir/Protac_params.txt" ]]; then
        echo "[ERROR] $poi: Protac_params.txt missing at $run_dir"; return 1
    fi

    # Safety check: Chains must be ABC D (not ABC A)
    chains_line=$(grep "^Chains:" "$run_dir/Protac_params.txt")
    if [[ "$chains_line" != "Chains: ABC D" ]]; then
        echo "[ERROR] $poi: Protac_params.txt has '$chains_line' — expected 'Chains: ABC D'"
        echo "        Run 07_prep_chain_d_reruns.py first, then re-transfer."
        return 1
    fi

    # Safety check: Results should be empty (prep script cleared them)
    n_res=$(ls "$run_dir/Results/"*.pdb 2>/dev/null | wc -l)
    if [[ $n_res -gt 50 ]]; then
        echo "[SKIP] $poi: Results/ still has $n_res PDBs — delete Results/*.pdb to re-run"
        return 0
    fi

    echo "[LAUNCH] $poi (${category}) → $log_file"
    (
        cd "$run_dir"
        conda run -n protac python "$SCRIPT" \
            Protac_params.txt \
            --n_dock    "$N_DOCK" \
            --n_workers "$N_WORKERS"
    ) > "$log_file" 2>&1 &
    echo "         PID $! → $log_file"
    echo $! >> "$LOG_DIR/pids.txt"
}

> "$LOG_DIR/pids.txt"

launch TNIK    low_affinity_degraders
launch MAPK14  low_affinity_degraders
launch MAP4K5  high_affinity_no_degradation
launch SLK     high_affinity_no_degradation

echo ""
echo "All 4 launched. PIDs: $(cat $LOG_DIR/pids.txt | tr '\n' ' ')"
echo ""
echo "Monitor:"
echo "  tail -f $LOG_DIR/TNIK.log"
echo "  tail -n5 $LOG_DIR/*.log | grep -E 'done|error|fail|cluster|complete'"
echo ""
echo "Check chain D in output (after run completes):"
echo "  python3 - <<'EOF'"
echo "  import subprocess, glob"
echo "  for poi in ['TNIK','MAPK14','MAP4K5','SLK']:"
echo "    cat = 'low_affinity_degraders' if poi in ['TNIK','MAPK14'] else 'high_affinity_no_degradation'"
echo "    pdbs = glob.glob(f'$PROTAC_DIR/pRosettaC/runs/{cat}/{poi}/Patchdock_Results/combined_*.pdb')[:1]"
echo "    if pdbs:"
echo "      chains = set(open(pdbs[0]).read()[21::80][:50])"
echo "      print(f'{poi}: sample chains in combined PDB = see file manually')"
echo "EOF"
echo ""
echo "=== Waiting for all 4 POIs to complete... ==="
wait
echo ""
echo "=== All done at $(date) ==="
echo ""

# Verify results
echo "=== Results count ==="
for poi_cat in "TNIK:low_affinity_degraders" "MAPK14:low_affinity_degraders" "MAP4K5:high_affinity_no_degradation" "SLK:high_affinity_no_degradation"; do
    poi="${poi_cat%%:*}"
    cat="${poi_cat##*:}"
    n=$(ls "$PROTAC_DIR/pRosettaC/runs/$cat/$poi/Results/"*.pdb 2>/dev/null | wc -l)
    echo "  $poi: $n Results PDBs"
done

echo ""
echo "If all ≥ 50 Results PDBs, rsync back to local:"
echo "  rsync -avz s1638308@ernie.bio.ed.ac.uk:~/protac/pRosettaC/runs/low_affinity_degraders/TNIK/ \\"
echo "    '/mnt/d/PROTAC degradation Prediction/pRosettaC/runs/low_affinity_degraders/TNIK/'"
echo "  rsync -avz s1638308@ernie.bio.ed.ac.uk:~/protac/pRosettaC/runs/low_affinity_degraders/MAPK14/ \\"
echo "    '/mnt/d/PROTAC degradation Prediction/pRosettaC/runs/low_affinity_degraders/MAPK14/'"
echo "  rsync -avz s1638308@ernie.bio.ed.ac.uk:~/protac/pRosettaC/runs/high_affinity_no_degradation/MAP4K5/ \\"
echo "    '/mnt/d/PROTAC degradation Prediction/pRosettaC/runs/high_affinity_no_degradation/MAP4K5/'"
echo "  rsync -avz s1638308@ernie.bio.ed.ac.uk:~/protac/pRosettaC/runs/high_affinity_no_degradation/SLK/ \\"
echo "    '/mnt/d/PROTAC degradation Prediction/pRosettaC/runs/high_affinity_no_degradation/SLK/'"
