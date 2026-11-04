#!/usr/bin/env python3
"""Prepare exact common foretinib and unmethylated VHL heads from bound parents.

No invented coordinates: the chosen fragments occur in both the supplier PROTAC
and the deposited bound ligand. Modified portions outside the head remain linker.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path
from rdkit import Chem
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError, atomic_json, sha256

CORE='COc1cc2c(Oc3ccc(NC(=O)C4(CC4)C(=O)Nc4ccc(F)cc4)cc3F)ccnc2cc1O'

def fragment(mol,indices):
    ids=sorted(indices);rw=Chem.RWMol(mol)
    for i in reversed(range(mol.GetNumAtoms())):
        if i not in ids:rw.RemoveAtom(i)
    result=rw.GetMol();Chem.SanitizeMol(result)
    return result,ids

def prepare(root):
    refs=root/'reference_structures/validated'
    parent=Chem.MolFromMolFile(str(refs/'foretinib_parent/bound_ligand.sdf'),removeHs=False)
    query=Chem.MolFromSmiles(CORE)
    matches=parent.GetSubstructMatches(query)
    if len(matches)!=1:raise ValidationError('Foretinib core mapping must be unique')
    poi,poi_ids=fragment(parent,matches[0])
    # Terminal phenoxy oxygen is the known solvent-exposed linker exit.
    exit_atoms=[a.GetIdx() for a in poi.GetAtoms() if a.GetSymbol()=='O' and a.GetDegree()==1
                and a.GetNeighbors()[0].GetIsAromatic()]
    if len(exit_atoms)!=1:raise ValidationError('Ambiguous foretinib exit oxygen')
    vhl_parent=Chem.MolFromMolFile(str(refs/'VH032_parent/bound_ligand.sdf'),removeHs=False)
    acetyl=Chem.MolFromSmarts('[CH3][C](=[O])[N]')
    match=vhl_parent.GetSubstructMatches(acetyl)
    if len(match)!=1:raise ValidationError('Ambiguous VHL N-acetyl cap')
    rw=Chem.RWMol(vhl_parent);rw.RemoveBond(match[0][1],match[0][3]);split=rw.GetMol();Chem.SanitizeMol(split)
    vhl,vhl_ids=fragment(vhl_parent,max(Chem.GetMolFrags(split),key=len))
    if vhl.GetNumAtoms()!=30:raise ValidationError('Unexpected VHL head size')
    registry=json.loads((root/'config/molecules.json').read_text())['compounds']
    for compound in ['SJF8240','SJF_delta']:
        full=Chem.MolFromSmiles(registry[compound]['smiles'])
        if not full.HasSubstructMatch(poi,useChirality=True):raise ValidationError('Foretinib core not in '+compound)
    for compound in ['SJF8240','SJF1521']:
        full=Chem.MolFromSmiles(registry[compound]['smiles'])
        if not full.HasSubstructMatch(vhl,useChirality=True):raise ValidationError('VHL head not in '+compound)
    dest=refs/'SJF8240';dest.mkdir(exist_ok=True)
    heads={'production_ready':False,'method':'exact_bound_parent_fragments',
           'poi':{'file':str((dest/'poi_head.sdf').relative_to(root)),
                  'anchor_one_based':exit_atoms[0]+1,'parent_indices_zero_based':poi_ids,
                  'reference_pdb':'3LQ8','component':'88Z','parent_sha256':sha256(refs/'foretinib_parent/bound_ligand.sdf')},
           'vhl':{'file':str((dest/'vhl_head.sdf').relative_to(root)),
                  'anchor_one_based':vhl_ids.index(match[0][3])+1,'parent_indices_zero_based':vhl_ids,
                  'reference_pdb':'4W9H','component':'3JF','parent_sha256':sha256(refs/'VH032_parent/bound_ligand.sdf')}}
    for role,mol in [('poi',poi),('vhl',vhl)]:
        w=Chem.SDWriter(str(dest/(role+'_head.sdf')));w.write(mol);w.close()
        heads[role]['sha256']=sha256(dest/(role+'_head.sdf'))
    shutil.copy2(refs/'foretinib_parent/protein.pdb',dest/'protein.pdb')
    shutil.copy2(refs/'VH032_parent/proteins.pdb',dest/'proteins.pdb')
    atomic_json(dest/'heads.json',heads)
    print(json.dumps(heads,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    prepare(p.parse_args().root)
