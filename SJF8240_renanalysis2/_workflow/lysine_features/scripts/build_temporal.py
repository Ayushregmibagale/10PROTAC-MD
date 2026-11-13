"""Predeclared 29 node features and four equal-pose time-window graphs."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from protac_qc import atomic_json
AA='ACDEFGHIKLMNPQRSTVWY';ALIASES={'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H','CYX':'C','MSE':'M'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 labels=json.loads((ROOT/'config/protocol.json').read_text())['labels'];completion=json.loads((ROOT/'analysis/lysines/completion.json').read_text());assert completion['cases']==141
 target_rows=[];lys_rows=[];graphs=[]
 for target,label in labels.items():
  records=[r for r in completion['records'] if r['poi']==target];packets=[]
  for r in records:
   p=ROOT/'analysis/lysines/poses'/target/r['pose']/'packet.npz';assert sha(p)==r['packet_sha256'];packets.append({k:v.copy() for k,v in np.load(p).items()})
  first=packets[0];lyscount=len(first['lys_resids'])
  for z in packets:
   for k in ['lys_resids','lys_node_indices','resids_A','resids_C','resids_D','resnames_A','resnames_C','resnames_D']:assert np.array_equal(z[k],first[k]),'Pose correspondence changed'
  # Equal pose aggregation: every pose has 41 sampled frames.
  perkeys={'complex_sasa_A2':'complex_sasa_A2','burial_A2':'burial_A2','min_e2_A':'min_e2_A','reach80_fraction':'reach80_fraction'}
  desc={k:np.mean([z[v].mean(0) for z in packets],axis=0) for k,v in perkeys.items()}
  desc['exposed_fraction']=np.mean([(z['complex_sasa_A2']>5).mean(0) for z in packets],axis=0)
  desc['joint_exposed_reachable_bridge_fraction']=np.mean([((z['complex_sasa_A2']>5)&(z['min_e2_A']<80)&z['bridge'][:,None]).mean(0) for z in packets],axis=0)
  for j,resid in enumerate(first['lys_resids']):lys_rows.append({'target':target,'label':label,'resid':int(resid),'poses':len(packets),**{k:float(v[j]) for k,v in desc.items()}})
  target_rows.append({'target':target,'label':label,'poses':len(packets),'lysines':lyscount,'mean_NZ_sasa_A2':float(desc['complex_sasa_A2'].mean()),'mean_NZ_burial_A2':float(desc['burial_A2'].mean()),'mean_exposed_fraction':float(desc['exposed_fraction'].mean()),'mean_joint_fraction':float(desc['joint_exposed_reachable_bridge_fraction'].mean()),'best_lysine_joint_fraction':float(desc['joint_exposed_reachable_bridge_fraction'].max()),'mean_min_e2_A':float(desc['min_e2_A'].mean())})
  for role,chain in [('elob','A'),('vhl','C')]:
   nl=len(first['resids_'+chain]);nr=len(first['resids_D']);base=np.zeros((nl+nr,29),dtype=np.float32)
   for offset,c,ri in [(0,chain,0),(nl,'D',1)]:
    ids=first['resids_'+c];base[offset:offset+len(ids),ri]=1;base[offset:offset+len(ids),2]=(ids-ids.min())/max(int(ids.max()-ids.min()),1)
    for j,name in enumerate(first['resnames_'+c]):
     a=seq1(str(name),custom_map=ALIASES,undef_code='X');assert a in AA,'Unknown amino acid '+str(name);base[offset+j,3+AA.index(a)]=1
   arr={}
   for w in range(4):
    mask=(first['times']>=5*w)&((first['times']<5*(w+1)) if w<3 else (first['times']<=20));assert mask.sum() in [10,11]
    feats=[];joints=[];lefts=[];rights=[]
    for z in packets:
     complex=z['complex_sasa_A2'][mask];burial=z['burial_A2'][mask];dist=z['min_e2_A'][mask];reach=z['reach80_fraction'][mask];bridge=z['bridge'][mask,None]
     feats.append(np.stack([np.ones(lyscount),complex.mean(0)/150,burial.mean(0)/150,dist.mean(0)/100,reach.mean(0),((complex>5)&(dist<80)&bridge).mean(0)],axis=1))
     l=z[role+'_left'][mask];r=z[role+'_right'][mask];joints.append(np.einsum('ti,tj->ij',(l<5).astype(float),(r<5).astype(float))/len(l));lefts.append(l.mean(0));rights.append(r.mean(0))
    x=base.copy();x[nl+first['lys_node_indices'],23:]=np.mean(feats,axis=0)
    joint=np.mean(joints,axis=0);left=np.mean(lefts,axis=0);right=np.mean(rights,axis=0);i,j=np.nonzero(joint>=.01);ei=np.array([np.r_[i,nl+j],np.r_[nl+j,i]],dtype=np.int64);ea=np.column_stack([joint[i,j],left[i]/10,right[j]/10]).astype(np.float32);ea=np.concatenate([ea,ea]);assert np.isfinite(x).all() and np.isfinite(ea).all()
    arr.update({f'x{w}':x,f'edge_index{w}':ei,f'edge_attr{w}':ea})
   dest=ROOT/'inputs/temporal'/role;dest.mkdir(parents=True,exist_ok=True);p=dest/(target+'.npz');np.savez_compressed(p,**arr);graphs.append({'role':role,'target':target,'sha256':sha(p),'path':str(p.relative_to(ROOT))})
 import pandas as pd
 dest=ROOT/'results';dest.mkdir(exist_ok=True);pd.DataFrame(target_rows).to_csv(dest/'lysine_target_summary.tsv',sep='\t',index=False);pd.DataFrame(lys_rows).to_csv(dest/'per_lysine_summary.tsv',sep='\t',index=False)
 atomic_json(ROOT/'config/temporal_input_manifest.json',{'graphs':graphs,'node_features':['role_E3','role_POI','normalized_resid']+list(AA)+['is_POI_lysine','NZ_complex_SASA/150','NZ_burial/150','mean_min_E2_distance/100','reach80_point_fraction','joint_exposed5A2_reachable80A_trueVHLbridge'],'windows_ns':[[0,5],[5,10],[10,15],[15,20]],'last_window_includes_20':True,'equal_pose_weight':True,'label_based_feature_selection':False,'missing_E2Ub_catalytic_geometry':True})
 print(pd.DataFrame(target_rows).to_string(index=False))
if __name__=='__main__':main()
