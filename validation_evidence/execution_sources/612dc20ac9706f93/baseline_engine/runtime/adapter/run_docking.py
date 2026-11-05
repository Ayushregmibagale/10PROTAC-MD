#!/usr/bin/env python3
"""
PyRosetta local docking script for pRosettaC.
Called as a subprocess by pyrosetta_adapter.local_docking().

Usage:
  python run_docking.py <struct.pdb> <chainsA> <chainsB> <paramsA> <paramsB> [nstruct]

Output files: <struct_prefix>_docking_NNNN.pdb  (Rosetta-style naming)

The docking protocol mirrors pRosettaC's docking.xml:
  1. DockingLowRes  (centroid global search)
  2. DockingHighRes (full-atom local refinement)
"""

import sys
import os
import pyrosetta
from pyrosetta_adapter import publication_seed
from pyrosetta.rosetta.protocols.docking import (
    DockingLowRes, DockMCMProtocol, setup_foldtree,
)
from pyrosetta.rosetta.utility import vector1_int


def main():
    if len(sys.argv) < 6:
        print('Usage: run_docking.py struct.pdb chainsA chainsB paramsA paramsB [nstruct]')
        sys.exit(1)

    struct   = sys.argv[1]
    chainsA  = sys.argv[2]  # e.g. 'ABCX'
    chainsB  = sys.argv[3]  # e.g. 'AY'
    paramsA  = sys.argv[4]
    paramsB  = sys.argv[5]
    nstruct  = int(sys.argv[6]) if len(sys.argv) > 6 else 50

    pyrosetta.init(
        f'-extra_res_fa {paramsA} '
        f'-extra_res_fa {paramsB} '
        f'-ignore_unrecognized_res false -constant_seed -jran {publication_seed("local_dock")} '
        '-out:level 100'
    )

    pose = pyrosetta.pose_from_file(struct)

    # Build partners string for foldtree setup, e.g. 'ABCX_AY'
    partners = f'{chainsA}_{chainsB}'
    jumps = vector1_int()
    jumps.append(1)
    setup_foldtree(pose, partners, jumps)

    scorefxn_low  = pyrosetta.create_score_function('interchain_cen')
    scorefxn_high = pyrosetta.create_score_function('ref2015')

    dock_low  = DockingLowRes(scorefxn_low, 1)
    dock_high = DockMCMProtocol(1, scorefxn_high, scorefxn_high)

    scorefxn = pyrosetta.create_score_function('ref2015')
    prefix = struct.replace('.pdb', '')

    for i in range(nstruct):
        wpose = pyrosetta.Pose(pose)
        dock_low.apply(wpose)
        dock_high.apply(wpose)
        out_pdb = f'{prefix}_docking_{i + 1:04d}.pdb'
        wpose.dump_pdb(out_pdb)
        score = scorefxn(wpose)
        desc = out_pdb.replace('.pdb', '')
        with open('local.fasc', 'a') as lf:
            lf.write(f'SCORE:  {score:.3f}  0  0  0  {score:.3f}  {desc}\n')


if __name__ == '__main__':
    main()
