#!/usr/bin/env python3
"""
PyRosetta local docking script for pRosettaC.
Called as a subprocess by pyrosetta_adapter.local_docking().

Usage:
  python run_docking.py <struct.pdb> <chainsA> <chainsB> <paramsA> <paramsB> [nstruct]

Output files: <struct_prefix>_docking_NNNN.pdb  (Rosetta-style naming)

SEPARATE UNQUALIFIED centroid+soft-anchor candidate; benchmark before adoption.
The docking stages are:
  1. DockingLowRes  (centroid global search)
  2. DockingHighRes (full-atom local refinement)
"""

import sys
import os
import subprocess, json, math
from pathlib import Path
import pyrosetta
from pyrosetta_adapter import publication_seed
from pyrosetta.rosetta.protocols.docking import (
    DockingLowRes, DockMCMProtocol, setup_foldtree,
)
from pyrosetta.rosetta.utility import vector1_int
from pyrosetta.rosetta.protocols.simple_moves import SwitchResidueTypeSetMover, ReturnSidechainMover
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring import atom_pair_constraint
from pyrosetta.rosetta.core.scoring.func import FlatHarmonicFunc
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint


def read_anchor_contract(raw):
    params = dict(line.split(':', 1) for line in (raw / 'Protac_params.txt').read_text().splitlines() if ':' in line)
    serials = [int(x) for x in params['Anchor atoms'].split()]
    if len(serials) != 2:
        raise ValueError('Expected exactly two anchor serials')
    names = []
    for i, serial in enumerate(serials):
        hits = [line for line in (raw / f'Init{i}.pdb').read_text().splitlines()
                if line.startswith('HETATM') and int(line[6:11]) == serial]
        if len(hits) != 1:
            raise ValueError('Ambiguous/missing anchor serial')
        names.append(hits[0][12:16].strip())
    lo, hi = json.loads((raw / 'anchor_distance_sampling.json').read_text())['range_A']
    if not (math.isfinite(lo) and math.isfinite(hi) and 0 <= lo < hi):
        raise ValueError('Invalid sampled anchor interval')
    return names, lo, hi


def anchor_atoms(pose, names):
    result = []
    for residue_name, atom_name in zip(('PT0', 'PT1'), names):
        indices = [i for i in range(1, pose.total_residue() + 1) if pose.residue(i).name3() == residue_name]
        if len(indices) != 1 or not pose.residue(indices[0]).has(atom_name):
            raise ValueError('Ambiguous/missing named head anchor')
        result.append(AtomID(pose.residue(indices[0]).atom_index(atom_name), indices[0]))
    return result


def add_anchor_restraint(pose, names, lo, hi):
    # Rebuild AtomIDs after residue-type conversion; do not assume atom-index identity.
    atoms = anchor_atoms(pose, names)
    pose.remove_constraints()
    pose.add_constraint(AtomPairConstraint(atoms[0], atoms[1],
        FlatHarmonicFunc((lo + hi) / 2, 1.0, (hi - lo) / 2)))


def anchor_distance(pose, names):
    atoms = anchor_atoms(pose, names)
    return float((pose.xyz(atoms[0]) - pose.xyz(atoms[1])).norm())


def restore_modified_residues(pose, original):
    restored = []
    for i in range(1, pose.total_residue() + 1):
        before = original.residue(i)
        after = pose.residue(i)
        expected = {before.atom_name(j).strip() for j in range(1, before.natoms() + 1)}
        observed = {after.atom_name(j).strip() for j in range(1, after.natoms() + 1)}
        if before.name() != after.name() or expected != observed:
            if not before.name().startswith('THR:phosphorylated'):
                raise ValueError('Unhandled residue restoration mismatch: ' + before.name())
            pose.replace_residue(i, before, True)
            restored.append(i)
            actual = {pose.residue(i).atom_name(j).strip() for j in range(1, pose.residue(i).natoms() + 1)}
            if actual != expected or not pose.residue(i).has('P'):
                raise ValueError('Failed phosphorylated Thr restoration')
    return restored


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

    # Candidate v1: explicit ligand centroid types and pose representation.
    # Full-atom parameters and all other docking/linker settings stay frozen.
    scratch = Path(Path(struct).stem + '_centroid_parameters').resolve()
    scratch.mkdir(exist_ok=False)
    centroid_params = []
    for head, parameter in [('PT0', paramsA), ('PT1', paramsB)]:
        sdf = Path(parameter).with_name(head + '_explicit_H.sdf')
        prefix = scratch / head
        command = [sys.executable, str(Path(__file__).parent / 'vendor/molfile_to_params.py'),
                   str(sdf), '-n', head, '-p', str(prefix), '--centroid', '--chain', 'X']
        result = subprocess.run(command, text=True, capture_output=True, timeout=120)
        prefix.with_suffix('.generation.log').write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError('Candidate centroid parameter generation failed')
        centroid_params.append(str(prefix) + '.cen.params')

    pyrosetta.init(
        f'-extra_res_fa {paramsA} '
        f'-extra_res_fa {paramsB} '
        f'-extra_res_cen {centroid_params[0]} -extra_res_cen {centroid_params[1]} '
        f'-ignore_unrecognized_res false -constant_seed -jran {publication_seed("local_dock")} '
        '-out:level 100'
    )

    names, lo, hi = read_anchor_contract(Path(paramsA).resolve().parent)
    pose = pyrosetta.pose_from_file(struct)

    # Build partners string for foldtree setup, e.g. 'ABCX_AY'
    partners = f'{chainsA}_{chainsB}'
    jumps = vector1_int()
    jumps.append(1)
    setup_foldtree(pose, partners, jumps)

    scorefxn_low  = pyrosetta.create_score_function('interchain_cen')
    scorefxn_high = pyrosetta.create_score_function('ref2015')
    scorefxn_low.set_weight(atom_pair_constraint, 1.0)
    scorefxn_high.set_weight(atom_pair_constraint, 1.0)

    dock_low  = DockingLowRes(scorefxn_low, 1)
    dock_high = DockMCMProtocol(1, scorefxn_high, scorefxn_high)

    scorefxn = pyrosetta.create_score_function('ref2015')
    prefix = struct.replace('.pdb', '')

    receipts = []
    for i in range(nstruct):
        wpose = pyrosetta.Pose(pose)
        original_fullatom = pyrosetta.Pose(wpose)
        restore = ReturnSidechainMover(original_fullatom)
        receipt = {'local_index': i + 1, 'before_A': anchor_distance(wpose, names)}
        SwitchResidueTypeSetMover('centroid').apply(wpose)
        if wpose.is_fullatom():
            raise RuntimeError('Candidate failed centroid conversion')
        add_anchor_restraint(wpose, names, lo, hi)
        dock_low.apply(wpose)
        receipt['after_low_A'] = anchor_distance(wpose, names)
        restore.apply(wpose)
        receipt['restored_modified_residues'] = restore_modified_residues(wpose, original_fullatom)
        if not wpose.is_fullatom():
            raise RuntimeError('Candidate failed full-atom restoration')
        add_anchor_restraint(wpose, names, lo, hi)
        dock_high.apply(wpose)
        receipt['after_high_A'] = anchor_distance(wpose, names)
        receipts.append(receipt)
        Path(prefix + '_anchor_restraint_receipt.json').write_text(json.dumps({'sampled_range_A': [lo, hi], 'soft_sd_A': 1.0, 'weight': 1.0, 'native_used': False, 'poses': receipts}, indent=2) + '\n')
        out_pdb = f'{prefix}_docking_{i + 1:04d}.pdb'
        wpose.dump_pdb(out_pdb)
        score = scorefxn(wpose)
        desc = out_pdb.replace('.pdb', '')
        with open('local.fasc', 'a') as lf:
            lf.write(f'SCORE:  {score:.3f}  0  0  0  {score:.3f}  {desc}\n')


if __name__ == '__main__':
    main()
