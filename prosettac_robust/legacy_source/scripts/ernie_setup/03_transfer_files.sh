#!/bin/bash
# Run this on YOUR LOCAL MACHINE (WSL) — NOT on ernie.
# Transfers the project and pRosettaC source to ernie.
#
# Usage: bash scripts/ernie_setup/03_transfer_files.sh
#
# What it transfers:
#   ~/protac/  (= ~/PROTAC degradation Prediction/) → ernie:~/protac/
#   ~/PRosettaC/                                    → ernie:~/PRosettaC/

set -e
ERNIE="s1638308@ernie.bio.ed.ac.uk"
LOCAL_PROJECT="$HOME/PROTAC degradation Prediction"
LOCAL_PROSETTAC="$HOME/PRosettaC"

echo "=== Transferring project to ernie:~/protac/ ==="
echo "This may take a few minutes (only large files: combined PDBs in Patchdock_Results)"
rsync -avz --progress \
  --exclude "Patchdock_Results/*.pdb" \
  --exclude "__pycache__" \
  --exclude "*.pyc" \
  --exclude ".git" \
  "$LOCAL_PROJECT/" \
  "$ERNIE:~/protac/"

echo ""
echo "=== Transferring pRosettaC source to ernie:~/PRosettaC/ ==="
rsync -avz --progress \
  --exclude "__pycache__" \
  --exclude "*.pyc" \
  --exclude ".git" \
  "$LOCAL_PROSETTAC/" \
  "$ERNIE:~/PRosettaC/"

echo ""
echo "=== Transfer complete ==="
echo "Next: run 04_fix_paths.sh ON ERNIE"
