"""Prepare graph-mapped MS143 heads; modified AKT head is explicitly modeled.

Requires downloaded 4GV1 and CCD 0XZ. No native ternary pose is invented.
"""
import argparse, hashlib, importlib.util, json, itertools, sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from Bio.PDB import MMCIFParser,PDBIO,Select
from Bio.PDB import PDBParser
from scipy.spatial.distance import cdist

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_sdf(p,m):
    with Chem.SDWriter(str(p)) as w:w.write(m)

def overlap_review(mol,protein):
    atoms=[x for x in PDBParser(QUIET=True).get_structure('p',str(protein)).get_atoms() if x.element not in {'H','D'}]
    pt=Chem.GetPeriodicTable();xyz=np.array([x.coord for x in atoms])
    radii=np.array([pt.GetRvdw(pt.GetAtomicNumber(x.element.title())) for x in atoms])
    distance=cdist(mol.GetConformer().GetPositions(),xyz)
    mradii=np.array([pt.GetRvdw(x.GetAtomicNum()) for x in mol.GetAtoms()])
    overlap=mradii[:,None]+radii[None,:]-distance
    i,j=np.unravel_index(np.argmax(overlap),overlap.shape)
    return {'max_vdw_overlap_A':float(overlap[i,j]),'distance_at_max_overlap_A':float(distance[i,j]),
            'ligand_atom_zero_based':int(i),'protein_atom':str(atoms[j].get_full_id()),'radii':'RDKit periodic table van der Waals radii'}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,required=True);p.add_argument('--sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-overlap',type=float,default=.6,help='Operational steric screen, assessed against the deposited parent; not a binding-affinity criterion')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    ref=load('bound_reference',a.project/'scripts/01_setup/00d_prepare_bound_references.py')
    analog=load('analog_head',a.project/'scripts/01_setup/00g_model_modified_heads.py')
    full=Chem.MolFromSmiles((a.project/'molecules/MS143/full.smi').read_text().strip())
    entry=json.loads((a.project/'config/molecules.json').read_text())['compounds']['MS143']
    if Chem.MolToInchiKey(full)!=entry['inchikey']:raise ValueError('Supplier identity mismatch')
    vhl=Chem.SDMolSupplier(str(a.project/'reference_structures/validated/SJF8240/vhl_head.sdf'),removeHs=True)[0]
    amides=full.GetSubstructMatches(Chem.MolFromSmarts('[N]-[C](=[O])'));choices=[]
    for bonds in itertools.combinations(amides,2):
        rw=Chem.RWMol(full)
        for n,c,_ in bonds:rw.RemoveBond(n,c)
        split=rw.GetMol();Chem.SanitizeMol(split)
        frags=Chem.GetMolFrags(split)
        if sorted(map(len,frags))!=[14,30,35]:continue
        vhids=next(x for x in frags if len(x)==30)
        vh=analog.fragment(full,vhids)
        if Chem.MolToInchiKey(vh)!=Chem.MolToInchiKey(vhl):continue
        choices.append((bonds,frags))
    if len(choices)!=1:raise ValueError('Ambiguous head/linker partition')
    bonds,frags=choices[0];poiids=next(x for x in frags if len(x)==35);vhids=next(x for x in frags if len(x)==30)
    poi=analog.fragment(full,poiids)
    if not full.HasSubstructMatch(poi,useChirality=True) or not full.HasSubstructMatch(vhl,useChirality=True):raise ValueError('Head stereo/graph mismatch')
    structure=MMCIFParser(QUIET=True,auth_chains=False).get_structure('4GV1',str(a.sources/'4GV1.cif'))
    parent,ligchain,atommap=ref.bound_molecule(structure,a.sources/'0XZ.cif','0XZ')
    # Select the protein chain adjacent to the actually extracted ligand.
    candidates=[]
    for chain in structure[0]:
        ats=[x.coord for r in chain if r.id[0]==' ' for x in r if x.element not in {'H','D'}]
        if ats:
            from scipy.spatial import cKDTree
            candidates.append((float(cKDTree(ats).query(parent.GetConformer().GetPositions())[0].min()),chain.id))
    candidates.sort()
    if not candidates or candidates[0][0]>4.5:raise ValueError('Bound protein chain not identified')
    selected=candidates[0][1]
    for chain in list(structure[0]):
        if chain.id!=selected:structure[0].detach_child(chain.id)
    structure[0][selected].id='D'
    modified=[{'name':r.resname,'number':r.id[1]} for r in structure[0]['D'] if r.resname=='TPO']
    class Protein(Select):
        def accept_model(self,m):return m.id==0
        def accept_residue(self,r):return r.id[0]==' ' or r.resname=='TPO'
        def accept_atom(self,x):return x.altloc in (' ','A')
    io=PDBIO();io.set_structure(structure);io.save(str(a.output/'AKT1_parent_protein.pdb'),Protein())
    write_sdf(a.output/'AZD5363_bound_parent.sdf',parent)
    write_sdf(a.output/'vhl_head.sdf',vhl)
    result={'compound':'MS143','catalog':'8081','supplier_inchikey':entry['inchikey'],
      'source_structure':'4GV1','component':'0XZ','ligand_label_chain':ligchain,'selected_original_protein_label_chain':selected,
      'source_sha256':{x:sha(a.sources/x) for x in ('4GV1.cif','0XZ.cif')},
      'full_atom_index_basis':'molecules/MS143/full.smi canonical graph; zero based',
      'cut_bonds':[[n,c] for n,c,_ in bonds],'poi_full_indices':list(poiids),'vhl_full_indices':list(vhids),
      'head_heavy_atoms':{'poi':35,'vhl':30,'linker_excluded_from_heads':14},
      'retained_modified_residues':modified,'production_ready':False,
      'required_reviews':['modeled piperazine binding geometry','reference phosphorylation/state','bound-head anchor correspondence','docking/QC integration'],
      'vhl_reference':'exact graph/stereo match to retained 4W9H-derived head; independent protein frame, not a ternary pose'}
    def anchor(ids):
        endpoints=[x for n,c,_ in bonds for x in (n,c) if x in ids]
        if len(endpoints)!=1:raise ValueError('Ambiguous anchor')
        return list(ids).index(endpoints[0])+1
    result['poi_anchor_one_based']=anchor(poiids)
    # Existing VHL coordinates/order may differ from supplier fragment ordering.
    vhlmatches=full.GetSubstructMatches(vhl,useChirality=True)
    vhlendpoint=next(x for n,c,_ in bonds for x in (n,c) if x in vhids)
    mapped_anchors={match.index(vhlendpoint)+1 for match in vhlmatches if vhlendpoint in match}
    if len(mapped_anchors)!=1:raise ValueError('Ambiguous VHL anchor')
    result['vhl_anchor_one_based']=mapped_anchors.pop()
    parent_check=overlap_review(parent,a.output/'AKT1_parent_protein.pdb')
    result['parent_steric_control']=parent_check;result['max_overlap_cutoff_A']=a.max_overlap
    if parent_check['max_vdw_overlap_A']>a.max_overlap:raise ValueError('Deposited parent fails chosen operational steric screen; review reference instead of tuning to outcomes')
    candidates=a.output/'candidates';candidates.mkdir(exist_ok=True)
    trials=[];passing=[]
    for seed in range(1,13):
        try:
            chosen,meta=analog.model(poi,parent,a.output/'AKT1_parent_protein.pdb',[seed])
            review=overlap_review(chosen,a.output/'AKT1_parent_protein.pdb')
            write_sdf(candidates/f'seed{seed}_head.sdf',chosen)
            good=review['max_vdw_overlap_A']<=a.max_overlap
            trials.append({'seed':seed,'coarse_model_pass':True,'steric_pass':good,**review,'UFF_energy_internal_only':meta['UFF_energy_internal_only']})
            if good:passing.append((meta['UFF_energy_internal_only'],seed,chosen,meta,review))
        except (ValueError,RuntimeError) as exc:trials.append({'seed':seed,'coarse_model_pass':False,'steric_pass':False,'reason':str(exc)})
    result['trials']=trials;result['strict_steric_pass_count']=len(passing)
    if passing:
        _,seed,chosen,meta,review=min(passing,key=lambda x:(x[0],x[1]))
        write_sdf(a.output/'poi_modeled_head.sdf',chosen)
        result.update(status='MODELED_HEAD_CANDIDATE_REQUIRES_REVIEW',poi_model=meta,selected_steric_review=review)
    else:
        result.update(status='BLOCKED_NO_STERICALLY_ACCEPTABLE_MODIFIED_HEAD',reason='All12fixed-scaffold trials failed convergence/coarse geometry or the operational vdW overlap screen; do not launch this placement')
    (a.output/'PREPARATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
