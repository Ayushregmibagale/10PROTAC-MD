#!/usr/bin/env python3
"""Build one verified pose-ensemble graph; no valid pose remains missing data.

Nine node features: role one-hot(3), within-role normalized position(1),
C/N/O/F/S one-hot(5). Three edge features: contact frequency, mean contact
 distance(nm), minimum contact distance(nm). Cross-role contacts use 8 Å.
This new version preserves model dimensions, not historical weight compatibility.
"""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,verify_receipt,pdb_atoms,sha256,atomic_json

def build(run,poi,label,output):
    if label not in (-1,0,1):raise ValidationError('Use an experimental binary label or -1 unknown')
    if output.exists():raise ValidationError('Output already exists')
    for old in ['/home/ayushregmibagale/protac-md','/mnt/d/PROTAC degradation Prediction']:
        if output.resolve().is_relative_to(Path(old)):raise ValidationError('Output points into original work')
    from rdkit import Chem
    mol=Chem.MolFromSmiles((run/'full.smi').read_text().split()[0])
    if mol is None or Chem.MolToInchiKey(mol)!='SCUZPTNDZHEHLB-ZWUDBOSPSA-N':raise ValidationError('Input is not verified SJF8240 chemistry')
    source_digest=sha256(Path(__file__));qc_path=Path(__file__).resolve().parents[1]/'lib/protac_qc.py';qc_digest=sha256(qc_path)
    selection=json.loads((run/'selection.json').read_text());rows=verify_receipt(run,selection)
    if not rows:
        output.mkdir(parents=True);atomic_json(output/'receipt.json',{'status':'NO_VALID_POSES_MISSING_GRAPH','poi':poi,'label':label,'biological_negative_inferred':False,'completion_sha256':sha256(run/'completion.json')});return
    reference=None;sums=None;counts=None;mins=None;features=None;roles=None
    for row in rows:
        p=run/'Patchdock_Results'/(row['description']+'.pdb');atoms=pdb_atoms(p);groups=[]
        for role in ['vhl','poi']:
            groups.append(sorted([a for a in atoms if a['record']=='ATOM' and a['chain'] in selection[role]['chains'] and a['name']=='CA'],key=lambda a:(a['chain'],a['resnum'],a['icode'])))
        groups.append(sorted([a for a in atoms if a['record']=='HETATM' and a['chain'] in selection['protac']['chains'] and a['resname']==selection['protac']['resname'] and a['element'] not in {'H','D'}],key=lambda a:a['name']))
        if any(not g for g in groups) or len(groups[2])!=79:raise ValidationError('Incomplete mapped SJF8240 graph nodes')
        keys=[(r,a['chain'],a['resnum'],a['icode'],a['name'],a['element']) for r,g in enumerate(groups) for a in g]
        if len(set(keys))!=len(keys):raise ValidationError('Nonunique graph atom identifiers')
        if reference is None:
            reference=keys;n=len(keys);roles=np.array([x[0] for x in keys]);counts=np.zeros((n,n),dtype=np.uint64);sums=np.zeros((n,n));mins=np.full((n,n),np.inf);features=np.zeros((n,9),dtype=np.float32)
            elements={'C':0,'N':1,'O':2,'F':3,'S':4};start=0
            for r,g in enumerate(groups):
                nums=np.array([a['resnum'] for a in g]) if r<2 else np.arange(len(g));norm=(nums-nums.min())/max(int(nums.max()-nums.min()),1)
                for j,a in enumerate(g):
                    if a['element'] not in elements:raise ValidationError('Unexpected node element')
                    features[start+j,r]=1;features[start+j,3]=norm[j];features[start+j,4+elements[a['element']]]=1
                start+=len(g)
        elif keys!=reference:raise ValidationError('Graph node correspondence changed across poses')
        xyz=np.array([a['xyz'] for g in groups for a in g]);distance=cdist(xyz,xyz);contact=(distance<8.)&(roles[:,None]!=roles[None,:]);counts+=contact;sums+=np.where(contact,distance,0.);mins=np.minimum(mins,np.where(contact,distance,np.inf))
    if sha256(Path(__file__))!=source_digest or sha256(qc_path)!=qc_digest or verify_receipt(run,selection)!=rows:raise ValidationError('Static graph source or inputs changed during extraction')
    freq=counts/len(rows);src,dst=np.nonzero(freq>=.005);edge=np.column_stack([freq[src,dst],sums[src,dst]/counts[src,dst]/10.,mins[src,dst]/10.]).astype(np.float32)
    output.mkdir(parents=True);np.savez_compressed(output/'graph.npz',x=features,edge_index=np.array([src,dst],dtype=np.int64),edge_attr=edge,label=np.array([label]),node_identifiers=np.array([json.dumps(x) for x in reference]),poi=np.array([poi]))
    atomic_json(output/'receipt.json',{'status':'GRAPH_PACKET_REQUIRES_REVIEW','feature_schema':'SJF8240-periodic-contacts-v1','component':'static','poi':poi,'label':label,'n_poses':len(rows),'nodes':len(reference),'directed_edges':len(src),'contact_cutoff_A':8.,'edge_frequency_threshold':.005,'mean_distance':'Conditional on contact, nm','normalization':'Within-role residue range; ligand sorted unique atom names','source_sha256':sha256(Path(__file__)),'completion_sha256':sha256(run/'completion.json'),'selection_sha256':sha256(run/'selection.json'),'graph_sha256':sha256(output/'graph.npz'),'automatic_training_release':False})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--poi',required=True);p.add_argument('--label',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();build(a.run,a.poi,a.label,a.output)
