#!/usr/bin/env python3
"""Rebuild an explicitly modeled p38δ control using the published template pair.

Smith et al., doi:10.1038/s41467-018-08027-7, used foretinib from 5IA4
overlaid on p38δ 4EYJ. Here correspondence is the explicit KLIFS pocket map;
this is an inspectable starting model, not a reproduction of their OPLS3/Jaguar
optimization or experimental p38δ/foretinib coordinates.
"""
import argparse
import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from protac_qc import ValidationError, atomic_json, min_distance, pdb_atoms, sha256


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def pocket_pairs(source_map, target_map, source_ca, target_ca):
    """Use explicit observed KLIFS indices; reject duplicates, never guess gaps."""
    def index(rows, atoms):
        result = {}
        for row in rows:
            key = int(row['index'])
            if key in result:
                raise ValidationError('Duplicate KLIFS pocket index')
            position = str(row['Xray_position'])
            result[key] = atoms.get(position)
        return result
    source = index(source_map, source_ca)
    target = index(target_map, target_ca)
    keys = sorted(k for k in source.keys() & target.keys()
                  if source[k] is not None and target[k] is not None)
    if len(keys) < 68:  # >=80% of the 85-position pocket
        raise ValidationError('Insufficient explicitly mapped KLIFS pocket coverage')
    return keys, np.array([source[k] for k in keys]), np.array([target[k] for k in keys])


def build(root):
    root = root.resolve()
    here = Path(__file__).resolve().parent
    bound = module('delta_bound_reference', here / '00d_prepare_bound_references.py')
    fragment = module('delta_foretinib_fragment', here / '00f_prepare_foretinib_heads.py')
    placement = module('delta_control_placement', here / '01_place_warhead_batch.py')
    refs = root / 'reference_sources'
    dest = root / 'modeled_reference_candidates/SJF_delta_MAPK13_published_overlay'
    control = root / 'prepared_inputs/SJF_delta/MAPK13_O15264_published_overlay'
    if dest.exists() or control.exists():
        raise ValidationError('Control already exists; preserve it and use a fresh root')
    paths = [refs / p for p in ['5IA4.pdb', '4EYJ.pdb', '88Z.cif',
             'KLIFS_7189_pocket_residue_map.json', 'KLIFS_2949_pocket_residue_map.json',
             'SJF_delta_original_study.xml']]
    source = [a for a in pdb_atoms(paths[0]) if a['record'] == 'ATOM' and a['chain'] == 'A']
    target = [a for a in pdb_atoms(paths[1]) if a['record'] == 'ATOM' and a['chain'] == 'A']
    def ca_map(atoms):
        result = {}
        for a in atoms:
            if a['name'] != 'CA':
                continue
            key = str(a['resnum']) + a['icode'].strip()
            if key in result:
                raise ValidationError('Ambiguous alternate Cα coordinates')
            result[key] = a['xyz']
        return result
    keys, x, y = pocket_pairs(json.loads(paths[3].read_text()), json.loads(paths[4].read_text()),
                              ca_map(source), ca_map(target))
    rotation, translation = placement.kabsch(x, y)
    rmsd = float(np.sqrt(np.mean(np.sum((x @ rotation.T + translation - y)**2, axis=1))))
    if rmsd > 2.0:
        raise ValidationError('Published-template pocket RMSD exceeds existing 2 Å threshold')
    structure = PDBParser(QUIET=True).get_structure('EPHA2_5IA4', paths[0])
    parent, _, _ = bound.bound_molecule(structure, paths[2], '88Z', 'A')
    matches = parent.GetSubstructMatches(Chem.MolFromSmiles(fragment.CORE), useChirality=True)
    if len(matches) != 1:
        raise ValidationError('Foretinib head must map uniquely')
    head, parent_ids = fragment.fragment(parent, matches[0])
    supplier = Chem.MolFromSmiles((root / 'molecules/SJF_delta/full.smi').read_text().strip())
    if not supplier.HasSubstructMatch(head, useChirality=True):
        raise ValidationError('Template head differs from supplier PROTAC')
    xyz = head.GetConformer().GetPositions() @ rotation.T + translation
    distance = min_distance(xyz, np.array([a['xyz'] for a in target if a['element'] not in {'H','D'}]))
    if not 1.5 <= distance <= 4.5:
        raise ValidationError('Transferred head violates existing contact/clash thresholds')
    exits = [a.GetIdx() for a in head.GetAtoms() if a.GetSymbol() == 'O'
             and a.GetDegree() == 1 and a.GetNeighbors()[0].GetIsAromatic()]
    if len(exits) != 1:
        raise ValidationError('Ambiguous linker exit atom')
    for i, position in enumerate(xyz):
        head.GetConformer().SetAtomPosition(i, position)
    dest.mkdir(parents=True)
    with Chem.SDWriter(str(dest / 'poi_head.sdf')) as writer:
        writer.write(head)
    original = root / 'modeled_reference_candidates/SJF_delta'
    for filename in ['proteins.pdb', 'vhl_head.sdf']:
        shutil.copy2(original / filename, dest / filename)
    (dest / 'protein.pdb').write_text('\n'.join(a['line'][:21]+'D'+a['line'][22:] for a in target)+'\nTER\nEND\n')
    evidence = {'status':'MODELED_STARTING_CONTROL_NOT_EXPERIMENTAL_COMPLEX',
        'paper_doi':'10.1038/s41467-018-08027-7', 'foretinib_template':'EPHA2/5IA4 chain A, component 88Z',
        'target_template':'MAPK13/4EYJ chain A', 'method':'KLIFS 85-pocket Cα least-squares overlay',
        'fit_positions':keys, 'pocket_coverage':len(keys)/85, 'pocket_rmsd_A':rmsd,
        'head_min_contact_A':distance, 'rotation':rotation.tolist(), 'translation':translation.tolist(),
        'source_sha256':{str(p.relative_to(root)):sha256(p) for p in paths},
        'limitations':['No experimental p38δ/foretinib binding coordinates',
          'No reproduction of published OPLS3/Jaguar torsion correction',
          'VHL head retains separately documented modeled conjugation',
          'No affinity, degradation or selectivity inference from preparation']}
    atomic_json(dest / 'model_evidence.json', evidence)
    heads = {'production_ready':False, 'method':'published_template_pair_with_explicit_KLIFS_overlay',
        'modeled_role':'poi_and_vhl', 'model_evidence':dict(evidence,
            file=str((dest/'model_evidence.json').relative_to(root)),sha256=sha256(dest/'model_evidence.json')),
        'poi':{'file':str((dest/'poi_head.sdf').relative_to(root)), 'sha256':sha256(dest/'poi_head.sdf'),
            'anchor_one_based':exits[0]+1, 'parent_indices_zero_based':parent_ids,
            'reference_pdb':'5IA4', 'target_pdb':'4EYJ', 'component':'88Z'},
        'vhl':dict(json.loads((original/'heads.json').read_text())['vhl'],
                   file=str((dest/'vhl_head.sdf').relative_to(root)),sha256=sha256(dest/'vhl_head.sdf'))}
    atomic_json(dest/'heads.json',heads)
    meta = placement.prepare(root, {'protac':'SJF_delta','gene':'MAPK13','uniprot_id':'O15264',
        'run_name':control.name,'structure_path':'reference_sources/4EYJ.pdb','chain':'A',
        'res_start':min(a['resnum'] for a in target),'res_end':max(a['resnum'] for a in target),
        'reference_bundle':str(dest.relative_to(root)), 'alignment_method':'sequence',
        'source':'published_template_pair_control'}, True)
    meta['source_path']='reference_sources/4EYJ.pdb'
    atomic_json(control/'preparation.json',meta)
    atomic_json(root/'audits/SJF_delta_MAPK13_control.json', {
        'status':'PREPARATION_CANDIDATE','production_ready':False,
        'preparation':str((control/'preparation.json').relative_to(root)),
        'preparation_sha256':sha256(control/'preparation.json'), 'evidence':evidence,
        'screen_effect':'Supplementary known-target control; existing 112 screen candidates unchanged'})
    print(json.dumps({'control':str(control),'pocket_rmsd_A':rmsd,'head_min_contact_A':distance,
                      'production_ready':False},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    build(parser.parse_args().root)
