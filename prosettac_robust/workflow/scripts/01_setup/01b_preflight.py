#!/usr/bin/env python3
"""Audit every selected input. Never connect to Ernie or submit simulations."""
import argparse
import csv
import json
import sys
from pathlib import Path
from rdkit import Chem
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError, read_params, pdb_atoms, sha256, min_distance, atomic_json
import numpy as np

def validate_pair(folder, registry):
    folder=Path(folder); issues=[]
    molecule=registry.get(folder.parent.name,{})
    if not molecule.get('smiles'): issues.append('UNRESOLVED_CHEMICAL_IDENTITY')
    try:
        params=read_params(folder/'Protac_params.txt')
        fullpath=Path(params['Protac'])
        fullpath=fullpath if fullpath.is_absolute() else folder/fullpath
        full=Chem.MolFromSmiles(fullpath.read_text().split()[0])
        if full is None: raise ValidationError('Invalid full-molecule SMILES')
        if molecule.get('inchikey') != Chem.MolToInchiKey(full): issues.append('WRONG_FULL_MOLECULE')
        for struct,chains,head,anchor in zip(params['Structures'].split(),params['Chains'].split(),
                                            params['Heads'].split(),map(int,params['Anchor atoms'].split())):
            structpath=Path(struct); headpath=Path(head)
            structpath=structpath if structpath.is_absolute() else folder/structpath
            headpath=headpath if headpath.is_absolute() else folder/headpath
            atoms=pdb_atoms(structpath)
            if set(chains)-{a['chain'] for a in atoms if a['record']=='ATOM'}:
                raise ValidationError('Requested protein chain absent')
            mol=Chem.MolFromMolFile(str(headpath),removeHs=False)
            if mol is None or mol.GetNumConformers()!=1 or not mol.GetConformer().Is3D():
                raise ValidationError('Missing bound 3D head')
            if anchor>mol.GetNumAtoms(): raise ValidationError('Head anchor index out of bounds')
            if not full.HasSubstructMatch(Chem.RemoveHs(mol),useChirality=True):
                raise ValidationError('Head not stereochemically compatible with full molecule')
            protein=np.array([a['xyz'] for a in atoms if a['chain'] in chains and a['record']=='ATOM' and a['element'] not in {'H','D'}])
            ligand=np.array([mol.GetConformer().GetAtomPosition(a.GetIdx()) for a in mol.GetAtoms() if a.GetAtomicNum()>1])
            distance=min_distance(ligand,protein)
            if distance<1.5 or distance>4.5: raise ValidationError(f'Implausible head contacts: {distance:.3f} Å')
    except (ValidationError,OSError,ValueError,IndexError) as exc:
        issues.append('INVALID_INPUT: '+str(exc))
    try:
        meta=json.loads((folder/'preparation.json').read_text())
        if not meta.get('files'): issues.append('MISSING_PREPARATION_HASHES')
        for filename,digest in meta.get('files',{}).items():
            if Path(filename).name!=filename or sha256(folder/filename)!=digest:
                issues.append('PREPARATION_HASH_MISMATCH: '+filename)
        if meta.get('molecule_inchikey')!=molecule.get('inchikey'): issues.append('PREPARATION_IDENTITY_MISMATCH')
        if not meta.get('production_ready'): issues.append('BENCHMARK_RELEASE_PENDING')
        if meta.get('contains_modeled_head'): issues.append('MODELED_HEAD_REVIEW_PENDING')
    except (OSError,ValueError) as exc: issues.append('MISSING_OR_INVALID_PREPARATION: '+str(exc))
    return {'protac':folder.parent.name,'run_name':folder.name,'status':'BLOCKED' if issues else 'READY',
            'issues':issues,'config_sha256':sha256(folder/'Protac_params.txt') if (folder/'Protac_params.txt').exists() else None}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--inputs',type=Path,help='Defaults to legacy prosettac_inputs, so bad inputs remain visible')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); registry=json.loads((a.root/'config/molecules.json').read_text())['compounds']
    inputs=a.inputs or a.root/'prosettac_inputs'
    folders=sorted(d for d in inputs.glob('*/*') if d.is_dir())
    results=[validate_pair(d,registry) for d in folders]
    atomic_json(a.output,{'schema_version':1,'input_root':str(inputs.resolve()),'pairs':results,
                         'ready':sum(r['status']=='READY' for r in results),'blocked':sum(r['status']=='BLOCKED' for r in results)})
    print(f'{len(results)} pairs: {sum(r["status"]=="READY" for r in results)} READY; {sum(r["status"]=="BLOCKED" for r in results)} BLOCKED')
    if not results or any(r['status']!='READY' for r in results): raise SystemExit(1)

if __name__=='__main__': main()
