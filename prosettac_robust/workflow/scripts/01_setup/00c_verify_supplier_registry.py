#!/usr/bin/env python3
"""Validate Tocris structures and export molecular graphs, without inventing poses."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError, atomic_json
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

def verify(root):
    registry = json.loads((root/'config/molecules.json').read_text())
    results = {}
    for name, entry in registry['compounds'].items():
        if not entry.get('smiles'):
            results[name] = {'status':'BLOCKED', 'reason':entry['reason']}; continue
        mol = Chem.MolFromSmiles(entry['smiles'])
        if mol is None: raise ValidationError('Invalid supplier SMILES: '+name)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        key = Chem.MolToInchiKey(mol)
        if formula != entry['formula'] or key != entry['inchikey']:
            raise ValidationError(f'{name}: supplier identity mismatch: {formula}, {key}')
        undefined = [idx for idx,label in Chem.FindMolChiralCenters(mol,includeUnassigned=True) if label == '?']
        if undefined: raise ValidationError(f'{name}: unspecified stereocenters {undefined}')
        dest = root/'molecules'/name
        dest.mkdir(parents=True,exist_ok=True)
        (dest/'full.smi').write_text(Chem.MolToSmiles(mol,isomericSmiles=True)+'\n')
        results[name] = {'status':'IDENTITY_VERIFIED', 'formula':formula, 'inchikey':key,
                         'molecular_weight':Descriptors.MolWt(mol), 'catalog':entry['catalog'],
                         'source':entry['source'], 'pose_status':entry['status']}
    atomic_json(root/'molecules/verification.json',results)
    print(json.dumps(results,indent=2))
    return results

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    verify(p.parse_args().root)
