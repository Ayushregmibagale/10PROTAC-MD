#!/usr/bin/env python3
"""Rebuild explicit internal sequence gaps and missing heavy atoms before MD.

Uses cached UniProt sequence identity and PDBFixer; retains observed termini and
the original ligand frame. Modeled residues/atoms and geometry checks are recorded.
Never use pdb2gmx -missing to conceal an incomplete protein model.
"""
import argparse
import io
import json
import os
import sys
from pathlib import Path
import numpy as np
from Bio.SeqUtils import seq3
from scipy.spatial import cKDTree
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import pdb_atoms,ValidationError,atomic_json,sha256
import importlib.util
spec=importlib.util.spec_from_file_location('domain_inventory',Path(__file__).with_name('00i_build_domain_inventory.py'))
inventory=importlib.util.module_from_spec(spec);spec.loader.exec_module(inventory)

def complete(root,pdb,poi_uid,output,seed):
    import openmm
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer
    os.environ['OPENMM_CPU_THREADS']='1'
    data=inventory.entries(root/'structure_provenance')
    chain_uid={'A':'P40337','B':'Q15370','C':'Q15369','D':poi_uid}
    all_atoms=pdb_atoms(pdb)
    atoms=[a for a in all_atoms if a['record']=='ATOM' and a['chain'] in chain_uid and a['element'] not in {'H','D'}]
    chains={c:[a for a in atoms if a['chain']==c] for c in chain_uid}
    if not all(chains.values()):raise ValidationError('Missing an explicitly identified protein chain')
    mappings={};expected={};identities={}
    for chain,group in chains.items():
        sequence=data[chain_uid[chain]]['sequence']['value']
        mapping,identity=inventory.map_sequence(group,sequence)
        ca=[a for a in group if a['name']=='CA']
        if len(mapping)<0.95*len(ca):raise ValidationError('Too many unmapped protein residues in chain '+chain)
        mappings[chain]=mapping;identities[chain]=identity;expected[chain]=sequence
    text=[]
    for chain,group in chains.items():
        # Intermediate gap termini were introduced deliberately for docking;
        # merge them for sequence-informed repair, retaining only true termini.
        final=max((a['resnum'],a['icode']) for a in group)
        text.extend(a['line'] for a in group if a['name']!='OXT' or (a['resnum'],a['icode'])==final)
        text.append('TER')
    fixer=PDBFixer(pdbfile=io.StringIO('\n'.join(text)+'\nEND\n'),platform=openmm.Platform.getPlatformByName('CPU'))
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:raise ValidationError('Nonstandard protein residues require an explicit chemical-state decision')
    fixer.missingResidues={};gaps=[]
    for chain in fixer.topology.chains():
        residues=list(chain.residues());mapping=mappings[chain.id];sequence=expected[chain.id]
        for position,(left,right) in enumerate(zip(residues,residues[1:]),1):
            a=mapping.get((int(left.id),left.insertionCode or ' '));b=mapping.get((int(right.id),right.insertionCode or ' '))
            if a is None or b is None:continue
            if b-a>1:
                missing=[seq3(letter).upper() for letter in sequence[a:b-1]]
                if len(missing)>30:raise ValidationError('Long/disordered gap requires a dedicated structural model: '+chain.id)
                fixer.missingResidues[(chain.index,position)]=missing
                gaps.append({'chain':chain.id,'left_uniprot':a,'right_uniprot':b,'modeled_residues':missing})
    fixer.findMissingAtoms()
    missing_atoms=[{'chain':r.chain.id,'residue':r.id,'resname':r.name,'atoms':[a.name for a in names]}
                   for r,names in fixer.missingAtoms.items()]
    fixer.addMissingAtoms(seed=seed)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w') as f:PDBFile.writeFile(fixer.topology,fixer.positions,f,keepIds=True)
    repaired=pdb_atoms(output)
    # Consecutive residues within each explicit chain must have a physical peptide bond.
    residues={}
    for a in repaired:residues.setdefault((a['chain'],a['resnum'],a['icode']),{})[a['name']]=a
    for chain in chains:
        groups=[group for key,group in residues.items() if key[0]==chain]
        for left,right in zip(groups,groups[1:]):
            if 'C' not in left or 'N' not in right:raise ValidationError('Incomplete repaired peptide backbone')
            d=float(np.linalg.norm(left['C']['xyz']-right['N']['xyz']))
            if not 1.0<d<1.7:raise ValidationError(f'Unclosed peptide bond after repair: {chain}, {d:.3f} A')
    before={(a['chain'],a['resnum'],a['icode'],a['name']):a['xyz'] for a in atoms}
    shifts=[float(np.linalg.norm(a['xyz']-before[(a['chain'],a['resnum'],a['icode'],a['name'])]))
            for a in repaired if (a['chain'],a['resnum'],a['icode'],a['name']) in before]
    if not shifts or max(shifts)>0.2:raise ValidationError('Repair moved observed protein atoms beyond the allowed frame tolerance')
    ligand=np.array([a['xyz'] for a in all_atoms if a['chain']=='X' and a['element'] not in {'H','D'}])
    if len(ligand):
        contact=float(cKDTree(np.array([a['xyz'] for a in repaired if a['element'] not in {'H','D'}])).query(ligand)[0].min())
        if contact<1.0:raise ValidationError('Protein repair created a severe ligand clash')
    else:contact=None
    report={'status':'PROTEIN_MODEL_COMPLETED','source_sha256':sha256(pdb),'output_sha256':sha256(output),
            'seed':seed,'chain_uniprot':chain_uid,'sequence_identity':identities,'internal_gaps':gaps,
            'missing_heavy_atoms':missing_atoms,'maximum_observed_atom_shift_A':max(shifts),
            'ligand_min_distance_A':contact,'new_coordinates_are_modeled':bool(gaps or missing_atoms),
            'openmm_version':openmm.__version__,'production_equilibration_completed':False}
    atomic_json(output.with_suffix('.repair.json'),report);print(json.dumps(report,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--pdb',type=Path,required=True);p.add_argument('--poi-uniprot',required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=int,required=True)
    a=p.parse_args();complete(a.root,a.pdb,a.poi_uniprot,a.output,a.seed)
