"""Inventory actual domain/pair coverage and independently review prepared heads.

This is read-only with respect to the source project. Existing candidate labels
are preserved; steric flags indicate review needs, not biological negatives.
"""
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
from rdkit import Chem
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write_table(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sterics(head_file,protein_file,chains=None):
    mol=Chem.SDMolSupplier(str(head_file),removeHs=True)[0]
    if mol is None or mol.GetNumConformers()!=1:raise ValueError('Missing valid head coordinates')
    atoms=[a for a in PDBParser(QUIET=True).get_structure('p',str(protein_file)).get_atoms()
           if a.element not in ('H','D') and a.get_parent().resname not in ('HOH','WAT') and (chains is None or a.get_parent().get_parent().id in chains)]
    if not atoms:raise ValueError('Empty protein selection')
    pt=Chem.GetPeriodicTable();xyz=np.array([a.coord for a in atoms]);tree=cKDTree(xyz)
    pr=np.array([pt.GetRvdw(pt.GetAtomicNumber(a.element.title())) for a in atoms])
    hx=mol.GetConformer().GetPositions();closest=tree.query(hx)[0]
    best=-float('inf');pair=None;backbone=-float('inf');sidechain=-float('inf');atom_name=None
    for i,near in enumerate(tree.query_ball_point(hx,4.5)):
        for j in near:
            distance=float(np.linalg.norm(hx[i]-xyz[j]));overlap=pt.GetRvdw(mol.GetAtomWithIdx(i).GetAtomicNum())+pr[j]-distance
            if atoms[j].name in {'N','CA','C','O','OXT'}:backbone=max(backbone,float(overlap))
            else:sidechain=max(sidechain,float(overlap))
            if overlap>best:best=float(overlap);pair=[i,str(atoms[j].get_full_id())];atom_name=atoms[j].name
    if pair is None:raise ValueError('No head-protein contacts within4.5A; review placement/frame instead of assigning a zero feature')
    return {'head_heavy_atoms':mol.GetNumAtoms(),'min_distance_A':float(closest.min()),
            'max_vdw_overlap_A':best if pair else None,'max_overlap_pair':pair,
            'head_atoms_within4A':int((closest<=4.).sum()),'max_overlap_protein_atom_name':atom_name,
            'max_backbone_vdw_overlap_A':backbone if np.isfinite(backbone) else None,
            'max_sidechain_vdw_overlap_A':sidechain if np.isfinite(sidechain) else None}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    root=a.project.resolve();a.output.mkdir(parents=True,exist_ok=True)
    registry=json.loads((root/'config/molecules.json').read_text())['compounds']
    verified=json.loads((root/'molecules/verification.json').read_text())
    invpath=root/'structure_provenance/domains/domain_inventory.tsv';domains=table(invpath)
    merged=table(root/'structures/merged_manifest.tsv')
    baselines={}
    for name in registry:
        if name=='MS143':continue
        folder=root/('reference_structures/validated' if name in ('GSK215','SJF8240') else 'modeled_reference_candidates')/name
        baselines[name]={'poi':sterics(folder/'poi_head.sdf',folder/'protein.pdb'),
                         'vhl':sterics(folder/'vhl_head.sdf',folder/'proteins.pdb',chains={'A','B','C'}),
                         'reference_folder':str(folder)}
    candidate_rows=[];method_summary=[];sources={str(invpath):sha(invpath),
       str(root/'config/molecules.json'):sha(root/'config/molecules.json'),
       str(root/'structures/merged_manifest.tsv'):sha(root/'structures/merged_manifest.tsv')}
    for baseline in baselines.values():
        for filename in ('poi_head.sdf','protein.pdb','vhl_head.sdf','proteins.pdb'):
            source=Path(baseline['reference_folder'])/filename;sources[str(source)]=sha(source)
    for method in ('prepared_screen','prepared_pocket_screen'):
        source=root/method/'screen_preparation.json';manifest=json.loads(source.read_text());sources[str(source)]=sha(source)
        if manifest.get('inventory_sha256',sha(invpath))!=sha(invpath):raise ValueError('Preparation inventory drift: '+method)
        for name in registry:
            selected=[r for r in manifest['pairs'] if r['protac']==name]
            candidates=[r for r in selected if r['status']=='CANDIDATE']
            method_summary.append({'method':method,'compound':name,'attempted_structure_pairs':len(selected),
               'candidate_structures':len(candidates),'candidate_unique_proteins':len({r['uniprot_id'] for r in candidates}),
               'rejected':sum(r['status']=='REJECTED' for r in selected),
               'correspondence_unavailable':sum(r['status']=='CORRESPONDENCE_UNAVAILABLE' for r in selected),
               'production_ready':False})
            for row in candidates:
                run=root/method/name/row['run_name'];run.resolve().relative_to(root)
                info=run/'preparation.json';integrity=[]
                if sha(info)!=row['preparation_sha256']:integrity.append('PREPARATION_HASH_MISMATCH')
                meta=json.loads(info.read_text())
                for filename,digest in meta['files'].items():
                    local=run/filename
                    if not local.is_file() or sha(local)!=digest:integrity.append('INPUT_HASH_MISMATCH:'+filename)
                if meta['molecule_inchikey']!=registry[name]['inchikey']:integrity.append('REGISTRY_IDENTITY_MISMATCH')
                actual=Chem.MolFromSmiles((run/'full.smi').read_text().strip())
                if actual is None or Chem.MolToInchiKey(actual)!=registry[name]['inchikey']:integrity.append('FULL_SMILES_IDENTITY_MISMATCH')
                entry={'method':method,'compound':name,'gene':row['gene'],'uniprot_id':row['uniprot_id'],
                  'domain_id':row['domain_id'],'run_name':row['run_name'],'original_status':row['status'],
                  'contains_modeled_head':meta['contains_modeled_head'],'pocket_rmsd_A':row['pocket_rmsd_A'],
                  'hash_identity_integrity':'PASS' if not integrity else ';'.join(integrity),
                  'preparation_sha256':sha(info),'poi_min_distance_A':None,'poi_max_vdw_overlap_A':None,
                  'poi_head_atoms_within4A':None,'vhl_min_distance_A':None,'vhl_max_vdw_overlap_A':None,
                  'vdw_review_flag':'NOT_EVALUATED','poi_reference_max_overlap_A':baselines[name]['poi']['max_vdw_overlap_A'],
                  'vhl_reference_max_overlap_A':baselines[name]['vhl']['max_vdw_overlap_A'],
                  'poi_overlap_change_from_reference_A':None,'vhl_overlap_change_from_reference_A':None,
                  'poi_worst_protein_atom':None,'poi_worst_contact_class':None,
                  'poi_max_backbone_overlap_A':None,'poi_backbone_overlap_change_A':None,
                  'poi_max_sidechain_overlap_A':None,'poi_sidechain_overlap_change_A':None,
                  'production_ready':False}
                try:
                    poi=sterics(run/'poi_head.sdf',run/'protein.pdb');vhl=sterics(run/'vhl_head.sdf',run/'vbc.pdb')
                    entry.update(poi_min_distance_A=poi['min_distance_A'],poi_max_vdw_overlap_A=poi['max_vdw_overlap_A'],
                         poi_head_atoms_within4A=poi['head_atoms_within4A'],vhl_min_distance_A=vhl['min_distance_A'],
                         vhl_max_vdw_overlap_A=vhl['max_vdw_overlap_A'],
                         vdw_review_flag='OVERLAP_REVIEW_NEEDED' if any(x['max_vdw_overlap_A'] is not None and x['max_vdw_overlap_A']>.6 for x in (poi,vhl)) else 'WITHIN_OPERATIONAL_SCREEN',
                         poi_overlap_change_from_reference_A=poi['max_vdw_overlap_A']-baselines[name]['poi']['max_vdw_overlap_A'],
                         vhl_overlap_change_from_reference_A=vhl['max_vdw_overlap_A']-baselines[name]['vhl']['max_vdw_overlap_A'],
                         poi_worst_protein_atom=poi['max_overlap_protein_atom_name'],
                         poi_worst_contact_class='BACKBONE' if poi['max_overlap_protein_atom_name'] in {'N','CA','C','O','OXT'} else 'SIDECHAIN',
                         poi_max_backbone_overlap_A=poi['max_backbone_vdw_overlap_A'],
                         poi_backbone_overlap_change_A=poi['max_backbone_vdw_overlap_A']-baselines[name]['poi']['max_backbone_vdw_overlap_A'] if poi['max_backbone_vdw_overlap_A'] is not None and baselines[name]['poi']['max_backbone_vdw_overlap_A'] is not None else None,
                         poi_max_sidechain_overlap_A=poi['max_sidechain_vdw_overlap_A'],
                         poi_sidechain_overlap_change_A=poi['max_sidechain_vdw_overlap_A']-baselines[name]['poi']['max_sidechain_vdw_overlap_A'] if poi['max_sidechain_vdw_overlap_A'] is not None and baselines[name]['poi']['max_sidechain_vdw_overlap_A'] is not None else None)
                except (ValueError,RuntimeError,OSError) as exc:entry['vdw_review_flag']='ERROR:'+str(exc)
                candidate_rows.append(entry)
    write_table(a.output/'COVERAGE_BY_COMPOUND.tsv',method_summary)
    write_table(a.output/'PREPARED_HEAD_REVIEW.tsv',candidate_rows)
    summary={'scope':'Inventory and static preparation review only; no new docking/MD; no biological labels inferred',
       'registered_compounds':list(registry),'supplier_verified':sum(v['status']=='IDENTITY_VERIFIED' for v in verified.values()),
       'merged_manifest_identifiers':len({r['uniprot_id'] for r in merged}),
       'defined_domain_structure_rows':len(domains),'proteins_with_defined_domains':len({r['uniprot_id'] for r in domains}),
       'missing_defined_domain_files':sum(not (root/r['structure_path']).is_file() for r in domains),
       'methods':{},'source_head_steric_baselines':baselines,'sources_sha256':sources,
       'limitations':['Candidate counts include structure variants; methods are alternative preparation cohorts, not additive biological samples',
         'An ATP-pocket fit is not evidence of binding across the kinome','Overlap0.6A is an operational review flag, not an experimentally calibrated binding threshold; native GSK215 VHL baseline itself exceeds it',
         'Continuous differences from the corresponding bound-head baseline distinguish inherited source contacts from transfer-induced overlaps',
         'Protein repacking/alternative conformations may resolve some overlaps; this review does not change or relax the source inputs',
         'No MS143 cross-kinase placements exist in these two preparation cohorts','Full campaign release pending reference recovery, preparation review, state coverage, resources and explicit user approval']}
    for method in ('prepared_screen','prepared_pocket_screen'):
        rows=[x for x in candidate_rows if x['method']==method]
        summary['methods'][method]={'candidate_structures':len(rows),'unique_compound_protein_pairs':len({(x['compound'],x['uniprot_id']) for x in rows}),
          'hash_identity_integrity':dict(Counter(x['hash_identity_integrity'] for x in rows)),
          'steric_review':dict(Counter(x['vdw_review_flag'] for x in rows)),
          'poi_overlap_increase_over_reference_counts':{str(margin):sum(x['poi_overlap_change_from_reference_A'] is not None and x['poi_overlap_change_from_reference_A']>margin for x in rows) for margin in (0.,.2,.5)},
          'vhl_overlap_increase_over_reference_gt0p2A':sum(x['vhl_overlap_change_from_reference_A'] is not None and x['vhl_overlap_change_from_reference_A']>.2 for x in rows),
          'worst_poi_contact_class':dict(Counter(x['poi_worst_contact_class'] for x in rows)),
          'poi_backbone_overlap_increase_over_reference_gt0p5A':sum(x['poi_backbone_overlap_change_A'] is not None and x['poi_backbone_overlap_change_A']>.5 for x in rows),
          'poi_sidechain_overlap_increase_over_reference_gt0p5A':sum(x['poi_sidechain_overlap_change_A'] is not None and x['poi_sidechain_overlap_change_A']>.5 for x in rows)}
    (a.output/'COVERAGE_AUDIT.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='sources_sha256'},indent=2))

if __name__=='__main__':main()
