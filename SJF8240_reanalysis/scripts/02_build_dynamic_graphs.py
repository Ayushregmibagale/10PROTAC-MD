"""Regenerate paired periodic joint-contact MD graphs from retained trajectories."""
from pathlib import Path
import json,sys,time,argparse,hashlib
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts/lib'))
from read_only_trajectory import ReadOnlyXTCReader
from protac_qc import ValidationError
from gnn_contract import periodic_minima
OLD=Path('/home/ayushregmibagale/protac-md');AUDIT=Path('/home/ayushregmibagale/10PROTAC-MD/prosettac_robust')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(2**20):h.update(b)
 return h.hexdigest()
def construct(x,joint,left,right):
 i,j=np.nonzero(joint>=.01);n=len(left);src=np.concatenate([i,n+j]);dst=np.concatenate([n+j,i]);a=np.column_stack([joint[i,j],left[i]/10,right[j]/10]).astype(np.float32)
 return x,np.array([src,dst],dtype=np.int64),np.concatenate([a,a])
def case(item):
 poi,pose=item;out=ROOT/'analysis/dynamic/poses'/poi/pose
 if (out/'receipt.json').exists():
  record=json.loads((out/'receipt.json').read_text())
  if any(sha(out/name)!=digest for name,digest in record['outputs_sha256'].items()):raise ValidationError('Cached dynamic output drift')
  return record
 if out.exists():raise ValidationError('Partial dynamic output requires explicit review')
 tpr=OLD/'random_screen'/poi/pose/'prod/prod.tpr';kind='kinome_md' if poi in ['SIK3','RPS6KB2'] else 'random_screen';xtc=(OLD/'random_screen'/poi/pose/'prod/prod.xtc') if kind=='kinome_md' else (OLD/'analysis'/kind/poi/pose/'prod_nojump.xtc')
 if not xtc.exists():raise ValidationError('Retained trajectory not found: '+str(xtc))
 paths=[tpr,xtc];hashes={str(p):sha(p) for p in paths};u=mda.Universe(str(tpr),str(xtc),format=ReadOnlyXTCReader)
 groups={r:u.select_atoms('moltype '+name).select_atoms('name CA') for r,name in [('elob','Protein_chain_A'),('vhl','Protein_chain_C'),('poi','Protein_chain_D')]}
 ligand=u.select_atoms('resname UNL PTC');ligand=ligand[ligand.masses>1.5]
 if len(ligand)!=79:raise ValidationError('Historical ligand heavy-atom count is not 79')
 cfg=json.loads((ROOT/'config/role_reference.json').read_text());uid=cfg['targets'].get(poi,cfg['targets'].get(poi+'_full'))
 # Bound exact per-case sequence assignments to the completed independent topology audit.
 core=json.loads((AUDIT/'existing_result_audit/protein_roles_core_complete/per_pose_roles.json').read_text())
 pred=json.loads((AUDIT/'existing_result_audit/prediction_compiled_roles_v3/per_pose_roles.json').read_text())
 evidence=next((x for x in core+pred if x['poi']==poi and x['pose_id']==pose),None)
 if evidence is None:raise ValidationError('No independent compiled role evidence for case')
 proteins=evidence['proteins'];proteins=list(proteins.values()) if isinstance(proteins,dict) else proteins
 for role,expected in [('elob','Q15370'),('vhl','P40337'),('poi',uid)]:
  protein=next(x for x in proteins if x['uniprot']==expected)
  if protein['moltype']!=str(groups[role].moltypes[0]):raise ValidationError('Protein molecule differs from independent sequence audit')
  if 'CA_atoms' in protein and protein['CA_atoms']!=len(groups[role]):raise ValidationError('CA count differs from sequence audit')
 if evidence.get('source_sha256',{}).get('tpr',hashes[str(tpr)])!=hashes[str(tpr)]:raise ValidationError('Compiled topology changed since identity audit')
 arrays={};npoi=len(groups['poi'])
 for role in ['elob','vhl']:arrays[role]={'joint':np.zeros((len(groups[role]),npoi)),'sum_left':np.zeros(len(groups[role])),'sum_right':np.zeros(npoi),'freq_left':np.zeros(len(groups[role])),'freq_right':np.zeros(npoi)}
 times=[];last_coordinates=None;duplicates=0;heavy_bridge=[]
 for ts in u.trajectory:
  t=float(ts.time)/1000
  if t>20.+1e-7:break
  if not np.isfinite(t) or t< -1e-7 or (times and t<times[-1]-1e-8):raise ValidationError('Invalid trajectory time order')
  if times and abs(t-times[-1])<1e-8:
   if not np.array_equal(ts.positions,last_coordinates):raise ValidationError('Repeated timestamp with different coordinates')
   duplicates+=1;continue
  right=periodic_minima(ligand.positions,groups['poi'].positions,ts.dimensions)
  for role,a in arrays.items():
   left=periodic_minima(ligand.positions,groups[role].positions,ts.dimensions)
   a['joint']+=np.outer(left<5.,right<5.);a['sum_left']+=left;a['sum_right']+=right;a['freq_left']+=left<5.;a['freq_right']+=right<5.
  times.append(t);last_coordinates=ts.positions.copy()
 if len(times)!=201 or not np.allclose(times,np.linspace(0,20,201),atol=5.1e-6,rtol=0):raise ValidationError('Explicit 201-frame 0–20 ns time grid unavailable')
 if any(sha(p)!=hashes[str(p)] for p in paths):raise ValidationError('Original dynamic source drift')
 out.mkdir(parents=True,exist_ok=False)
 for role,a in arrays.items():
  a={k:v/len(times) for k,v in a.items()};a.update(left_resids=np.array(groups[role].resids),right_resids=np.array(groups['poi'].resids),times_ns=np.array(times))
  np.savez_compressed(out/(role+'.npz'),**a)
 record={'status':'DYNAMIC_CASE_COMPLETE','poi':poi,'pose_id':pose,'frames':len(times),'duplicate_identical_frames_removed':duplicates,'source_sha256':hashes,'sequence_identity_evidence':str((AUDIT/'existing_result_audit/protein_roles_core_complete/per_pose_roles.json') if kind=='random_screen' else (AUDIT/'existing_result_audit/prediction_compiled_roles_v3/per_pose_roles.json')),'outputs_sha256':{p.name:sha(p) for p in out.glob('*.npz')},'source_script_sha256':sha(__file__),'new_MD':False}
 (out/'receipt.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'poi':poi,'pose':pose,'frames':len(times)}),flush=True);return record
def aggregate(records):
 for poi in sorted({r['poi'] for r in records}):
  name=poi+'_full' if poi in ['SIK3','RPS6KB2'] else poi
  for role in ['elob','vhl']:
   sums=None;ref=None;n=0
   for row in records:
    if row['poi']!=poi:continue
    path=ROOT/'analysis/dynamic/poses'/poi/row['pose_id']/(role+'.npz')
    if sha(path)!=row['outputs_sha256'][path.name]:raise ValidationError('Derived pose drift')
    with np.load(path) as z:
     ids=(z['left_resids'],z['right_resids'])
     if ref is None:ref=tuple(x.copy() for x in ids)
     elif any(not np.array_equal(a,b) for a,b in zip(ref,ids)):raise ValidationError('Residue correspondence differs across poses')
     if sums is None:sums={k:z[k].copy() for k in ['joint','sum_left','sum_right']}
     else:
      for k in sums:sums[k]+=z[k]
     n+=1
   for k in sums:sums[k]/=n
   x=np.zeros((len(ref[0])+len(ref[1]),3),dtype=np.float32);start=0
   for ri,ids in enumerate(ref):
    x[start:start+len(ids),ri]=1;x[start:start+len(ids),2]=(ids-ids.min())/max(int(ids.max()-ids.min()),1);start+=len(ids)
   x,ei,ea=construct(x,sums['joint'],sums['sum_left'],sums['sum_right'])
   dest=ROOT/'analysis/dynamic/targets'/name/role;dest.mkdir(parents=True,exist_ok=False)
   labels={r['poi']:r['historical_label'] for r in json.loads((ROOT/'config/source_inventory.json').read_text())['records']}
   label=labels.get(name);np.savez_compressed(dest/'graph.npz',x=x,edge_index=ei,edge_attr=ea,poi=np.array([name]),label=np.array([label if label is not None else -1]))
   (dest/'receipt.json').write_text(json.dumps({'status':'DYNAMIC_TARGET_COMPLETE','poi':name,'poses':n,'frames_per_pose':201,'role':role,'equal_pose_weight':True,'joint_frequency_not_product_of_marginals':True,'graph_sha256':sha(dest/'graph.npz'),'nodes':len(x),'directed_edges':len(ea)},indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=2);p.add_argument('--core-only',action='store_true');a=p.parse_args()
 roster=json.loads((AUDIT/'sjf8240_review/md_graph_pose_lineage.json').read_text())
 print('roster type',type(roster).__name__,flush=True)
 # Actual roster layout is validated below; never infer missing cases.
 if isinstance(roster,dict):roster=roster.get('poses',roster.get('rows',roster.get('records',[])))
 items=[(r['poi'],r['pose_id']) for r in roster]
 if a.core_only:items=[x for x in items if x[0] not in ['SIK3','RPS6KB2']]
 if len(set(items))!=len(items) or not items:raise ValidationError('Nonempty unique source roster required')
 from concurrent.futures import ProcessPoolExecutor
 with ProcessPoolExecutor(max_workers=a.workers) as pool:records=list(pool.map(case,items))
 aggregate(records)
 (ROOT/'analysis/dynamic/completion.json').write_text(json.dumps({'status':'DYNAMIC_COMPLETE','cases':len(records),'records':records,'new_MD':False},indent=2)+'\n')
