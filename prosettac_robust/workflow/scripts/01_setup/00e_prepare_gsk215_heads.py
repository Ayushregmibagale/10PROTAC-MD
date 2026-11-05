#!/usr/bin/env python3
"""Derive GSK215 head fragments from its native bound CCD atom identities.

Cut N4-C20 and C21-N5, leaving the central CH2-C(=O) linker out of both heads.
The explicit CCD names, original indices and anchor indices are archived.
These fragments are preparation candidates until the exact engine is benchmarked.
"""
import argparse
import json
import sys
from pathlib import Path
from rdkit import Chem
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError, atomic_json, sha256

def prepare(root):
    dest=root/'reference_structures/validated/GSK215'
    source=dest/'bound_ligand.sdf'
    mol=Chem.MolFromMolFile(str(source),removeHs=False)
    mapping=json.loads((dest/'provenance.json').read_text())['atom_indices_zero_based']
    # CCD labels are used rather than relying on a supplier's atom ordering.
    cuts=[('N4','C20'),('C21','N5')]
    rw=Chem.RWMol(mol)
    for a,b in cuts:
        if rw.GetBondBetweenAtoms(mapping[a],mapping[b]) is None:
            raise ValidationError('CCD attachment bond missing: '+a+'-'+b)
        if rw.GetBondBetweenAtoms(mapping[a],mapping[b]).IsInRing():
            raise ValidationError('Attachment cut would break a ring')
        rw.RemoveBond(mapping[a],mapping[b])
    split=rw.GetMol(); Chem.SanitizeMol(split)
    fragments=Chem.GetMolFrags(split)
    if sorted(map(len,fragments)) != [3,31,36]:
        raise ValidationError('Unexpected GSK215 head/linker partition')
    info={'source_sha256':sha256(source),'cut_bonds_ccd':cuts,'production_ready':False}
    full=Chem.MolFromSmiles((root/'molecules/GSK215/full.smi').read_text().strip())
    for role,anchor in [('poi','N4'),('vhl','N5')]:
        ids=next(ids for ids in fragments if mapping[anchor] in ids)
        part=Chem.RWMol(mol)
        for i in reversed(range(mol.GetNumAtoms())):
            if i not in ids: part.RemoveAtom(i)
        head=part.GetMol(); Chem.SanitizeMol(head)
        if not full.HasSubstructMatch(head,useChirality=True):
            raise ValidationError('Head is not a stereochemically compatible full-molecule substructure')
        file=dest/(role+'_head.sdf')
        w=Chem.SDWriter(str(file)); w.write(head); w.close()
        info[role]={'file':str(file.relative_to(root)), 'anchor_one_based':list(ids).index(mapping[anchor])+1,
                    'anchor_ccd':anchor, 'full_ccd_indices_zero_based':list(ids), 'sha256':sha256(file)}
    atomic_json(dest/'heads.json',info); print(json.dumps(info,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    prepare(p.parse_args().root)
