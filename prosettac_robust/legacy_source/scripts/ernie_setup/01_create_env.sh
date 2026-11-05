#!/bin/bash
# Run this ON ERNIE: creates the conda environment for pRosettaC
# Usage: bash 01_create_env.sh

set -e
echo "=== Creating protac conda environment (Python 3.10) ==="

conda create -y -n protac python=3.10

echo "=== Installing packages ==="
conda run -n protac conda install -y -c conda-forge \
    rdkit \
    biopython \
    numpy \
    pandas \
    scipy \
    scikit-learn

echo "=== Verifying ==="
conda run -n protac python -c "
from rdkit import Chem
from Bio.PDB import Superimposer
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
print('All packages OK')
print('  rdkit:', Chem.__version__ if hasattr(Chem,'__version__') else 'ok')
import Bio; print('  biopython:', Bio.__version__)
print('  numpy:', np.__version__)
print('  pandas:', pd.__version__)
"

echo ""
echo "=== DONE ==="
echo "Next step: install PyRosetta (see 02_install_pyrosetta.sh)"
