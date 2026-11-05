#!/bin/bash
# Run this ON ERNIE to verify the full setup before launching runs.
# Usage: bash ~/protac/scripts/ernie_setup/05_verify_setup.sh

PASS=0; FAIL=0

check() {
    local desc="$1"; local cmd="$2"
    if eval "$cmd" &>/dev/null; then
        echo "  [OK]  $desc"; ((PASS++))
    else
        echo "  [FAIL] $desc"; ((FAIL++))
    fi
}

echo "=== pRosettaC Setup Verification ==="
echo ""

echo "-- Environment --"
check "conda env 'protac' exists" "conda env list | grep -q '^protac '"
check "Python 3.10 in protac env" "conda run -n protac python -c 'import sys; assert sys.version_info[:2]==(3,10)'"
check "RDKit importable" "conda run -n protac python -c 'from rdkit import Chem'"
check "BioPython importable" "conda run -n protac python -c 'from Bio.PDB import Superimposer'"
check "PyRosetta importable" "conda run -n protac python -c 'import pyrosetta'"
check "scikit-learn importable" "conda run -n protac python -c 'from sklearn.linear_model import LogisticRegression'"

echo ""
echo "-- Project files --"
check "~/protac directory exists" "test -d ~/protac"
check "scripts/rosetta_utils present" "test -d ~/protac/scripts/rosetta_utils"
check "scripts/analysis present" "test -d ~/protac/scripts/analysis"
check "MET Protac_params.txt present" "test -f ~/protac/pRosettaC/runs/high_affinity_degraders/MET/Protac_params.txt"
check "VHL structure present" "test -f ~/protac/data/structures/VHL/4W9H_VBC_clean.pdb"
check "MET protein_clean.pdb present" "test -f ~/protac/data/structures/POI/high_affinity_degraders/MET/MET_protein_clean.pdb"
check "SJF8240 SMILES present" "test -f ~/protac/data/structures/PROTAC/SJF8240/SJF8240.smi"
check "PRosettaC directory present" "test -d ~/PRosettaC"
check "PRosettaC main.py present" "test -f ~/PRosettaC/main.py"

echo ""
echo "-- Path patching --"
check "No old paths in scripts" "! grep -r '/home/s1638308' ~/protac/scripts/ 2>/dev/null | grep -qv '.pyc'"
check "No old paths in params files" "! grep -r '/home/s1638308' ~/protac/pRosettaC/runs/ 2>/dev/null | grep -q '.'"

echo ""
echo "-- Warhead SDFs present --"
for poi in MET DDR2 RIPK2 EPHB2 MAPK14 AXL SLK ABL1 EPHA2 MAP4K5; do
    f=$(find ~/protac/data/structures/POI -name "${poi}_warhead_placed_H.sdf" 2>/dev/null | head -1)
    if [[ -n "$f" ]]; then
        echo "  [OK]  ${poi}_warhead_placed_H.sdf"
        ((PASS++))
    else
        echo "  [FAIL] ${poi}_warhead_placed_H.sdf not found"
        ((FAIL++))
    fi
done

echo ""
echo "-- Quick pRosettaC import test --"
check "pRosettaC importable from Python" \
  "conda run -n protac python -c 'import sys; sys.path.insert(0,\"$HOME/PRosettaC\"); import utils'"

echo ""
echo "=============================="
echo "PASSED: $PASS   FAILED: $FAIL"
echo "=============================="
if [[ $FAIL -eq 0 ]]; then
    echo "All checks passed — ready to launch runs!"
    echo "Run: bash ~/protac/scripts/ernie_setup/06_launch_parallel.sh"
else
    echo "Fix failed checks before launching."
fi
