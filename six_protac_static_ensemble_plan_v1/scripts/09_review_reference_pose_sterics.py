"""Independent atom-type-aware steric readouts for the two reference ensembles.

Native baselines and continuous differences are reported. No threshold is
tuned to degradation labels, and no structural prediction is called binding.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial import cKDTree

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def overlap(x,rx,y,ry):
    if not len(x) or not len(y):raise ValueError('Empty atom selection')
    tree=cKDTree(y);near=tree.query_ball_point(x,4.5)
    ii=np.repeat(np.arange(len(x)),[len(v) for v in near])
    if len(ii):
        jj=np.concatenate([np.asarray(v,dtype=int) for v in near if v])
        distances=np.linalg.norm(x[ii]-y[jj],axis=1);values=rx[ii]+ry[jj]-distances
        maximum=float(values.max())
    else:maximum=None
    return {'min_heavy_distance_A':float(tree.query(x)[0].min()),'max_vdw_overlap_within4p5A':maximum}

def measure(path,full):
    lines=[x for x in path.read_text().splitlines() if x.startswith(('ATOM  ','HETATM')) and x[16:17] in (' ','A')]
    if any(x[17:20].strip() in {'TPO','SEP','PTR','MSE'} for x in lines):raise ValueError('Modified residue requires an explicit extended role schema; this reference-only review does not silently exclude it')
    protein=[x for x in lines if x.startswith('ATOM  ') and x[76:78].strip().upper() not in {'H','D'}]
    if {x[21] for x in protein}!={'A','B','C','D'}:raise ValueError('Unexpected canonical reference protein chains')
    ligand=[x for x in lines if x.startswith('HETATM') and x[21]=='X' and x[17:20]=='PTC' and x[76:78].strip().upper() not in {'H','D'}]
    if not ligand:raise ValueError('Expected PTC chainX absent')
    mol=Chem.MolFromPDBBlock('\n'.join(ligand)+'\nEND\n',sanitize=False,removeHs=True,proximityBonding=True)
    if mol is None or mol.GetNumAtoms()!=full.GetNumAtoms():raise ValueError('Full ligand heavy-atom count mismatch')
    mol=AllChem.AssignBondOrdersFromTemplate(full,mol);Chem.AssignStereochemistryFrom3D(mol,replaceExistingTags=True)
    if Chem.MolToInchiKey(mol)!=Chem.MolToInchiKey(full):raise ValueError('Full ligand graph/stereochemistry mismatch')
    xyz=np.array([[float(x[30:38]),float(x[38:46]),float(x[46:54])] for x in protein])
    pt=Chem.GetPeriodicTable();pr=np.array([pt.GetRvdw(pt.GetAtomicNumber(x[76:78].strip().title())) for x in protein])
    ligand_xyz=mol.GetConformer().GetPositions();lr=np.array([pt.GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    lp=overlap(ligand_xyz,lr,xyz,pr);poi=np.array([x[21]=='D' for x in protein])
    pp=overlap(xyz[poi],pr[poi],xyz[~poi],pr[~poi])
    topology=Chem.GetDistanceMatrix(mol);mask=np.triu(topology>3,k=1)
    d=np.linalg.norm(ligand_xyz[:,None,:]-ligand_xyz[None,:,:],axis=2)
    intra=float((lr[:,None]+lr[None,:]-d)[mask].max()) if mask.any() else None
    return {'ligand_stereo_graph_match':True,'ligand_protein_min_distance_A':lp['min_heavy_distance_A'],
            'ligand_protein_max_vdw_overlap_A':lp['max_vdw_overlap_within4p5A'],
            'E3_POI_min_distance_A':pp['min_heavy_distance_A'],'E3_POI_max_vdw_overlap_A':pp['max_vdw_overlap_within4p5A'],
            'intraligand_max_vdw_overlap_topological_distance_gt3_A':intra}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--references',type=Path,required=True);p.add_argument('--recovery-table',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    with a.recovery_table.open() as f:recovery=list(csv.DictReader(f,delimiter='\t'))
    rows=[];baselines={};delta_keys=('ligand_protein_max_vdw_overlap_A','E3_POI_max_vdw_overlap_A','intraligand_max_vdw_overlap_topological_distance_gt3_A')
    for ref in ('5T35_MZ1','7PI4_GSK215'):
        run=a.runs/f'global_{ref}_seed42';full=Chem.MolFromSmiles((run/'full.smi').read_text().strip())
        native=a.references/ref/'native_complex.pdb';baselines[ref]=measure(native,full);baselines[ref]['native_sha256']=sha(native)
        for row in [r for r in recovery if r['reference']==ref]:
            pose=run/'Patchdock_Results'/(row['pose']+'.pdb')
            out={**row,'independent_steric_status':'NOT_EVALUATED','ligand_stereo_graph_match':None,
                 'ligand_protein_min_distance_A':None,'ligand_protein_max_vdw_overlap_A':None,
                 'E3_POI_min_distance_A':None,'E3_POI_max_vdw_overlap_A':None,
                 'intraligand_max_vdw_overlap_topological_distance_gt3_A':None,
                 **{'change_from_native_'+k:None for k in delta_keys},'error':''}
            try:
                out.update(measure(pose,full));out['independent_steric_status']='MEASURED_NOT_A_RELEASE_CERTIFICATE'
                for key in delta_keys:
                    out['change_from_native_'+key]=out[key]-baselines[ref][key] if out[key] is not None and baselines[ref][key] is not None else None
            except (ValueError,RuntimeError,OSError) as exc:out.update(independent_steric_status='FAILED_REVIEW',error=str(exc))
            rows.append(out)
    with (a.output/'REFERENCE_STERIC_REVIEW.tsv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    summary={}
    for ref in baselines:
        rr=[r for r in rows if r['reference']==ref]
        summary[ref]={'poses':len(rr),'chemistry_or_review_errors':sum(r['independent_steric_status']=='FAILED_REVIEW' for r in rr),
            'ligand_protein_overlap_increase_gt0p5A':sum(r['change_from_native_ligand_protein_max_vdw_overlap_A'] is not None and r['change_from_native_ligand_protein_max_vdw_overlap_A']>.5 for r in rr),
            'E3_POI_overlap_increase_gt0p5A':sum(r['change_from_native_E3_POI_max_vdw_overlap_A'] is not None and r['change_from_native_E3_POI_max_vdw_overlap_A']>.5 for r in rr),
            'top_ranked':min(rr,key=lambda r:int(r['rank_by_raw_total_score'])) if rr else None}
    report={'status':'INDEPENDENT_CONTINUOUS_STERIC_REVIEW_COMPLETE','native_baselines':baselines,'summary':summary,
       'limitations':['RDKit element-based van der Waals radii; not a force-field energy or binding affinity',
          '0.5A excess counts are exploratory diagnostic summaries, not a validated acceptance cutoff',
          'Protein interaction checks include ELOB/ELOC as well as VHL','Intraligand overlap excludes bonded, angle and1-4pairs (topological distance<=3)',
          'No modification/minimization or selective deletion of poses','This reference-only schema explicitly rejects unhandled modified residues'],
       'source_script_sha256':sha(Path(__file__)),'recovery_table_sha256':sha(a.recovery_table)}
    (a.output/'REFERENCE_STERIC_REVIEW.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'baselines':baselines,'summary':{k:{x:y for x,y in v.items() if x!='top_ranked'} for k,v in summary.items()}},indent=2))

if __name__=='__main__':main()
