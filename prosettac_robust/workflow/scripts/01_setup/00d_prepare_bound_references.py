#!/usr/bin/env python3
"""Extract deposited bound coordinates with CCD bonds and explicit chain mapping.

Exports preparation controls, not approved production inputs. Supplier molecular
identity and native ternary geometry are independent checks.
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError, atomic_json, sha256, pose_features

REFERENCES = {
    'GSK215': {'pdb':'7PI4','ccd':'7QB','poi_chain':'D','vhl_chain':'A', 'ligand_chain':'E',
               'chains':{'A':'A','B':'B','C':'C','D':'D'}},
    'SJF1521_parent': {'pdb':'1XKK','ccd':'FMM','poi_chain':'A','ligand_chain':'D','chains':{'A':'D'}},
    'JNJ1013_parent': {'pdb':'6BFN','ccd':'DL1','poi_chain':'A','ligand_chain':None,'chains':{'A':'D'}},
    'foretinib_parent': {'pdb':'3LQ8','ccd':'88Z','poi_chain':'A','ligand_chain':'B','chains':{'A':'D'}},
    'VH032_parent': {'pdb':'4W9H','ccd':'3JF','ligand_chain':'M',
                     'focus_chain':'A','chains':{'C':'A','A':'B','B':'C'}},
}

def bound_molecule(structure, ccd, component, ligand_chain=None):
    d=MMCIF2Dict(str(ccd)); atom_ids=d['_chem_comp_atom.atom_id']
    elements=d['_chem_comp_atom.type_symbol']
    # Select a single ligand residue and reject missing/ambiguous heavy atoms.
    residues=[r for c in structure[0] for r in c if r.resname==component
              and (ligand_chain is None or c.id==ligand_chain)]
    if not residues: raise ValidationError('No bound ligand '+component)
    residue=sorted(residues,key=lambda r:(r.parent.id,r.id))[0]
    present={a.name:a for a in residue}
    rw=Chem.RWMol(); mapping={}
    for i,(name,element) in enumerate(zip(atom_ids,elements)):
        if element in {'H','D'}: continue
        if name not in present: raise ValidationError('Missing ligand atom '+name)
        a=Chem.Atom(element.title())
        charge=d.get('_chem_comp_atom.charge',['0']*len(atom_ids))[i]
        a.SetFormalCharge(int(charge)); a.SetProp('ccd_atom_id',name)
        mapping[name]=rw.AddAtom(a)
    orders={'SING':Chem.BondType.SINGLE,'DOUB':Chem.BondType.DOUBLE,'TRIP':Chem.BondType.TRIPLE,'AROM':Chem.BondType.AROMATIC}
    for a,b,order in zip(d['_chem_comp_bond.atom_id_1'],d['_chem_comp_bond.atom_id_2'],d['_chem_comp_bond.value_order']):
        if a in mapping and b in mapping: rw.AddBond(mapping[a],mapping[b],orders[order])
    mol=rw.GetMol(); Chem.SanitizeMol(mol)
    conf=Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
    for name,idx in mapping.items(): conf.SetAtomPosition(idx,present[name].coord.astype(float))
    mol.AddConformer(conf); Chem.AssignStereochemistryFrom3D(mol)
    return mol, residue.parent.id, mapping

def prepare(root, sources):
    registry=json.loads((root/'config/molecules.json').read_text())['compounds']
    results={}
    for name,cfg in REFERENCES.items():
        cif=sources/(cfg['pdb']+'.cif'); ccd=sources/(cfg['ccd']+'.cif')
        structure=MMCIFParser(QUIET=True,auth_chains=False).get_structure(cfg['pdb'],str(cif))
        mol,ligand_chain,mapping=bound_molecule(structure,ccd,cfg['ccd'],cfg['ligand_chain'])
        key=Chem.MolToInchiKey(mol)
        if name=='GSK215' and key != registry[name]['inchikey']:
            raise ValidationError('7PI4 ligand differs from Tocris GSK215: '+key)
        dest=root/'reference_structures/validated'/name
        dest.mkdir(parents=True,exist_ok=True)
        w=Chem.SDWriter(str(dest/'bound_ligand.sdf')); w.write(mol); w.close()
        class ProteinOnly(Select):
            def accept_model(self,model): return model.id==0
            def accept_chain(self,chain): return chain.id in cfg['chains'].values()
            def accept_residue(self,residue): return residue.id[0]==' '
            def accept_atom(self,atom): return atom.altloc in (' ','A')
        # Normalize using label chain IDs, not multi-character author IDs in 7PI4.
        for chain in list(structure[0]):
            if chain.id not in cfg['chains']: structure[0].detach_child(chain.id)
        renamed={chain.id:chain for chain in structure[0]}
        for old,chain in renamed.items(): chain.id='temporary_'+old
        for old,chain in renamed.items(): chain.id=cfg['chains'][old]
        io=PDBIO(); io.set_structure(structure); io.save(str(dest/'proteins.pdb'),ProteinOnly())
        # Individual POI file uses canonical chain D, for explicit later mapping.
        proteins=(dest/'proteins.pdb').read_text().splitlines()
        focus=cfg.get('focus_chain','D')
        (dest/'protein.pdb').write_text('\n'.join(l for l in proteins if l.startswith(('ATOM  ','TER')) and l[21:22]==focus)+'\nEND\n')
        meta={'kind': 'native_ternary_preparation_control' if name=='GSK215' else 'parent_ligand_only',
              'pdb':cfg['pdb'],'component':cfg['ccd'],'source_sha256':sha256(cif),
              'ccd_sha256':sha256(ccd),'bound_inchikey':key,'ligand_label_chain':ligand_chain,
              'atom_indices_zero_based':mapping,'production_ready':False}
        if name=='GSK215':
            selection={'vhl':{'chains':['A'],'record':'ATOM'},'poi':{'chains':['D'],'record':'ATOM'},
                       'protac':{'chains':['X'],'record':'HETATM','resname':'PTC'}}
            ligand_lines=[]
            for i,a in enumerate(mol.GetAtoms(),1):
                x,y,z=mol.GetConformer().GetAtomPosition(i-1)
                ligand_lines.append(f'HETATM{i:5d} {a.GetSymbol()+str(i):>4s} PTC X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a.GetSymbol():>2s}  ')
            (dest/'native_complex.pdb').write_text('\n'.join(l for l in proteins if l.startswith('ATOM  '))+'\n'+ '\n'.join(ligand_lines)+'\nEND\n')
            meta['native_geometry']=pose_features(dest/'native_complex.pdb',selection)
            atomic_json(dest/'selection.json',selection)
        atomic_json(dest/'provenance.json',meta); results[name]=meta
    print(json.dumps(results,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--sources',type=Path,required=True)
    a=p.parse_args(); prepare(a.root,a.sources)
