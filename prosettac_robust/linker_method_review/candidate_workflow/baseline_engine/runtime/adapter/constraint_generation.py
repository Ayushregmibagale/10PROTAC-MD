#!/usr/bin/env python3
"""
Wrapper for pRosettaC's constraint_generation.py that applies our PyRosetta patches
before importing the original script.  main.py builds subprocess commands pointing to
utils.SCRIPTS_FOL/constraint_generation.py; we redirect SCRIPTS_FOL to this directory
so this wrapper is invoked instead of the unpatched original.
"""
import sys
import os

# Apply all adapter patches (PyRosetta, OpenBabel-Python, sparse GenConstConf, etc.)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pyrosetta_adapter import apply_patches
apply_patches()

# Now run the original constraint_generation with a patched environment
sys.path.insert(0, str(os.environ['PROTAC_PROSETTAC']))
import importlib
cg = importlib.import_module('constraint_generation')
cg.main(sys.argv[0], sys.argv[1:])
