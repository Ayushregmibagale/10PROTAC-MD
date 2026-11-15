"""Evaluate retained docking outputs with the previously frozen reference gate.

No docking or input coordinates are modified. Native self/rigid-body controls
check this measurement only, not predictive recovery.
"""
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def residues(path, chain):
    model=PDBParser(QUIET=True).get_structure('s',str(path))[0]
    if chain not in model: raise ValueError(f'Missing chain {chain}: {path}')
    rs=[r for r in model[chain] if r.id[0]==' ' and 'CA' in r]
    keys=[(r.id[1],r.id[2],r.resname) for r in rs]
    if len(keys)!=len(set(keys)): raise ValueError('Ambiguous residue identities')
    return keys,np.array([r['CA'].coord for r in rs],dtype=float),[
        np.array([a.coord for a in r if a.element not in {'H','D'}],dtype=float) for r in rs]

def transform(moving,fixed):
    mc=moving.mean(0);fc=fixed.mean(0)
    u,_,vt=np.linalg.svd((moving-mc).T@(fixed-fc))
    d=np.eye(3);d[-1,-1]=np.linalg.det(u@vt)
    rotation=u@d@vt
    return rotation,fc-mc@rotation

def contact_pairs(a,b,cutoff=5.):
    xyz=np.concatenate(b);owners=np.concatenate([np.full(len(x),i,dtype=int) for i,x in enumerate(b)])
    tree=cKDTree(xyz);pairs=set()
    for i,atoms in enumerate(a):
        for near in tree.query_ball_point(atoms,cutoff):
            pairs.update((i,int(owners[j])) for j in near)
    return pairs

def measure(native,model,cutoff=5.):
    for role in ('vhl','poi'):
        if native[role][0]!=model[role][0]:raise ValueError('Residue identity/order differs: '+role)
    r,t=transform(model['vhl'][1],native['vhl'][1])
    rms=lambda x,y:float(np.sqrt(np.mean(np.sum((x-y)**2,axis=1))))
    nc=contact_pairs(native['vhl'][2],native['poi'][2],cutoff)
    if not nc:raise ValueError('No native contacts')
    pc=contact_pairs(model['vhl'][2],model['poi'][2],cutoff)
    return {'VHL_fit_CA_RMSD_A':rms(model['vhl'][1]@r+t,native['vhl'][1]),
      'POI_CA_RMSD_after_VHL_fit_A':rms(model['poi'][1]@r+t,native['poi'][1]),
      'native_residue_contact_fraction':len(nc&pc)/len(nc),'native_contact_pairs':len(nc),'model_contact_pairs':len(pc)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--references',type=Path,help='Separate evaluation-only <reference>/native_complex.pdb files; never required as generation inputs')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    rows=[];summaries=[]
    for ref in ('5T35_MZ1','7PI4_GSK215'):
        run=a.runs/f'global_{ref}_seed42';state=json.loads((run/'engine_status.json').read_text())
        if state['status'] not in ('COMPLETE','COMPLETE_NO_VALID_POSES'):raise ValueError('Not a complete retained-output run')
        selection=json.loads((run/'selection.json').read_text())
        chains={role:selection[role]['chains'][0] for role in ('vhl','poi')}
        if any(len(selection[x]['chains'])!=1 for x in chains):raise ValueError('Expected one chain per role')
        native_path=a.references/ref/'native_complex.pdb' if a.references else run/'native_complex.pdb'
        native={role:residues(native_path,chain) for role,chain in chains.items()}
        selfcheck=measure(native,native)
        if selfcheck['POI_CA_RMSD_after_VHL_fit_A']>1e-8 or selfcheck['native_residue_contact_fraction']!=1:raise ValueError('Native measurement control failed')
        # Common rigid transformation must not change recovery.
        q=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
        translated={k:(v[0],v[1]@q+np.array([13.,-7.,4.]),[x@q+np.array([13.,-7.,4.]) for x in v[2]]) for k,v in native.items()}
        rc=measure(native,translated)
        if rc['POI_CA_RMSD_after_VHL_fit_A']>1e-7 or rc['native_residue_contact_fraction']!=1:raise ValueError('Rigid measurement control failed')
        incompatible={k:(v[0],v[1].copy(),[x.copy() for x in v[2]]) for k,v in native.items()}
        v=incompatible['poi'];incompatible['poi']=(v[0],v[1]+1000,[x+1000 for x in v[2]])
        neg=measure(native,incompatible)
        if neg['native_residue_contact_fraction']!=0 or neg['POI_CA_RMSD_after_VHL_fit_A']<100:raise ValueError('Incompatible measurement control failed')
        score=run/'Patchdock_Results/score.sc';items=[]
        for line in score.read_text().splitlines():
            s=line.split()
            if len(s)==3 and s[0]=='SCORE:' and s[1]!='total_score':items.append((float(s[1]),s[2]))
        if len(items)!=state['stages']['pose_qc']['succeeded']:raise ValueError('Retained score/QC count mismatch')
        current=[]
        for rank,(energy,name) in enumerate(sorted(items),1):
            pose=run/'Patchdock_Results'/f'{name}.pdb'
            m={role:residues(pose,chain) for role,chain in chains.items()}
            metrics=measure(native,m)
            row={'reference':ref,'seed':42,'rank_by_raw_total_score':rank,'pose':name,'raw_total_score':energy,
                 **metrics,'frozen_gate_pass':metrics['POI_CA_RMSD_after_VHL_fit_A']<=5 and metrics['native_residue_contact_fraction']>=.5,
                 'pose_sha256':sha(pose)}
            rows.append(row);current.append(row)
        summaries.append({'reference':ref,'retained_QC_poses':len(current),'raw_scored_poses':state['stages']['scoring']['succeeded'],
          'gate':'POI CA RMSD after VHL fit <=5A AND native residue contact fraction >=0.5; retained engine chemical/gross-clash QC',
          'gate_pass_count':sum(x['frozen_gate_pass'] for x in current),'top1_pass':bool(current and current[0]['frozen_gate_pass']),
          'top10_any_pass':any(x['frozen_gate_pass'] for x in current[:10]),
          'best_POI_RMSD_A':min((x['POI_CA_RMSD_after_VHL_fit_A'] for x in current),default=None),
          'best_native_contact_fraction':max((x['native_residue_contact_fraction'] for x in current),default=None),
          'run_seconds':state['finished_unix']-state['started_unix'],
          'measurement_controls':'native self, common rigid transform, incompatible target placement: passed',
          'engine_status_sha256':sha(run/'engine_status.json'),'score_sha256':sha(score),'native_sha256':sha(native_path),
          'biological_validation':False})
    with (a.output/'REFERENCE_POSE_RECOVERY.tsv').open('w') as f:
        fields=['reference','seed','rank_by_raw_total_score','pose','raw_total_score','VHL_fit_CA_RMSD_A',
                'POI_CA_RMSD_after_VHL_fit_A','native_residue_contact_fraction','native_contact_pairs',
                'model_contact_pairs','frozen_gate_pass','pose_sha256']
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    report={'status':'EVALUATED_FROZEN_GATE','summary':summaries,
      'limitations':['One seed at a small budget; not a convergence study','Engine chemistry/gross-clash QC inherited; no new independent chemistry assay',
                    'This is a structural benchmark, not a degradation predictor','Native coordinates not used for generation in these global runs']}
    (a.output/'REFERENCE_RECOVERY.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
