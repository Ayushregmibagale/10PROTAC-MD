#!/bin/bash
# Run this ON ERNIE after downloading the PyRosetta wheel.
#
# HOW TO GET THE WHEEL:
#   1. Go to https://www.pyrosetta.org/downloads#h.6vttn15ac69d
#   2. Log in with your university email (free academic license)
#   3. Download: PyRosetta4 Release for Python-3.10 Linux
#      File will be named something like:
#      PyRosetta4.Release.python310.linux.release-NNN.whl
#   4. Transfer to ernie:
#      From your LOCAL terminal:
#        scp ~/Downloads/PyRosetta4.Release.python310.linux.*.whl \
#            s1638308@ernie.bio.ed.ac.uk:~/
#
# Usage: bash 02_install_pyrosetta.sh PyRosetta4.Release.python310.linux.release-NNN.whl

set -e
WHEEL="${1:-}"
if [[ -z "$WHEEL" ]]; then
    echo "Usage: bash 02_install_pyrosetta.sh <path-to-wheel.whl>"
    echo "Looking for any .whl file in ~/ ..."
    WHEEL=$(ls ~/PyRosetta*.whl 2>/dev/null | head -1)
    if [[ -z "$WHEEL" ]]; then
        echo "No wheel found. Download from https://www.pyrosetta.org/downloads"
        exit 1
    fi
    echo "Found: $WHEEL"
fi

echo "=== Installing PyRosetta into protac env ==="
conda run -n protac pip install "$WHEEL"

echo "=== Verifying PyRosetta ==="
conda run -n protac python -c "
import pyrosetta
pyrosetta.init('-out:level 0')
print('PyRosetta OK:', pyrosetta.__version__)
"

echo "=== PyRosetta installed successfully ==="
