#!/bin/bash
# Run this ON ERNIE after 03_transfer_files.sh completes.
# Replaces all hardcoded /home/s1638308 paths with /home/s1638308
# and rebuilds the symlink expected by PyRosetta.
#
# Usage: bash ~/protac/scripts/ernie_setup/04_fix_paths.sh

set -e
OLD_HOME="/home/s1638308"
NEW_HOME="/home/s1638308"
PROTAC_DIR="$HOME/protac"
PROSETTAC_DIR="$HOME/PRosettaC"

echo "=== Patching paths in Python scripts ==="
find "$PROTAC_DIR/scripts" -name "*.py" | while read f; do
    if grep -q "$OLD_HOME" "$f" 2>/dev/null; then
        sed -i "s|$OLD_HOME|$NEW_HOME|g" "$f"
        echo "  patched: $(basename $f)"
    fi
done

echo ""
echo "=== Patching paths in shell scripts ==="
find "$PROTAC_DIR/scripts" -name "*.sh" | while read f; do
    if grep -q "$OLD_HOME" "$f" 2>/dev/null; then
        sed -i "s|$OLD_HOME|$NEW_HOME|g" "$f"
        echo "  patched: $(basename $f)"
    fi
done

echo ""
echo "=== Patching Protac_params.txt files ==="
find "$PROTAC_DIR/pRosettaC/runs" -name "Protac_params.txt" | while read f; do
    if grep -q "$OLD_HOME" "$f" 2>/dev/null; then
        sed -i "s|$OLD_HOME|$NEW_HOME|g" "$f"
        echo "  patched: $f"
    fi
done

echo ""
echo "=== Patching VHL/POI structure paths in data files ==="
find "$PROTAC_DIR/data" -name "*.txt" -o -name "*.params" | while read f; do
    if grep -q "$OLD_HOME" "$f" 2>/dev/null; then
        sed -i "s|$OLD_HOME|$NEW_HOME|g" "$f"
        echo "  patched: $(basename $f)"
    fi
done

echo ""
echo "=== Creating ~/protac symlink (for PyRosetta path compatibility) ==="
# The scripts use ~/protac to avoid spaces; on ernie ~/protac IS the directory
# so no symlink needed — but we verify it exists
if [[ -d "$PROTAC_DIR" ]]; then
    echo "  ~/protac directory exists OK"
else
    echo "  ERROR: ~/protac not found — did 03_transfer_files.sh complete?"
    exit 1
fi

echo ""
echo "=== Creating PatchDock symlink (if needed) ==="
# pRosettaC expects PatchDock in ~/PRosettaC/PatchDock/
if [[ -d "$PROSETTAC_DIR/PatchDock" ]]; then
    echo "  PatchDock found at $PROSETTAC_DIR/PatchDock"
else
    echo "  WARNING: PatchDock not found in $PROSETTAC_DIR/PatchDock"
    echo "  You may need to recompile PatchDock on ernie:"
    echo "    cd ~/PRosettaC/PatchDock && make"
fi

echo ""
echo "=== Verifying patched scripts ==="
remaining=$(grep -rl "$OLD_HOME" "$PROTAC_DIR/scripts" "$PROTAC_DIR/pRosettaC/runs" 2>/dev/null | wc -l)
if [[ "$remaining" -eq 0 ]]; then
    echo "  All paths patched successfully ✓"
else
    echo "  WARNING: $remaining files still contain old paths:"
    grep -rl "$OLD_HOME" "$PROTAC_DIR/scripts" "$PROTAC_DIR/pRosettaC/runs" 2>/dev/null | head -10
fi

echo ""
echo "=== Path fix complete ==="
echo "Next: run 05_verify_setup.sh to confirm everything works"
