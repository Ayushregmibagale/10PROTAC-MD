"""Observed true-VHL MD bridging; reuse audited core and compute33 extension cases."""
from pathlib import Path
import os
os.environ['OMP_NUM_THREADS']='1'
import sys,json,concurrent.futures,hashlib
import numpy as np,pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
ROOT=Path(__file__).resolve().parents[1];CAN=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis');AUDIT=CAN.parent/'prosettac_robust/existing_result_audit';sys.path.insert(0,str(CAN/'full_frame_e2_delta50_v1/lib'))
from read_only_trajectory import ReadOnlyXTCReader
from protac_qc import sha256,atomic_json

def require(ok,msg):
 if not ok:raise ValueError(msg)
def summarize(item,distances,source):
 times=np.asarray(distances['time_ns']);require(len(times)==201 and np.allclose(times,np.arange(201)/10,atol=5.1e-5,rtol=0),'Incomplete frame grid');v=np.asarray(distances['true_VHL_heavy_A']);k=np.asarray(distances['POI_heavy_A']);b=(v<5)&(k<5);runs=[];start=None
 for i,state in enumerate(np.r_[b,False]):
  if state and start is None:start=i
  elif not state and start is not None:runs.append({'frames':i-start,'observed_span_ns':float(times[i-1]-times[start]),'left_censored':start==0,'right_censored':i==len(b)});start=None
 return {'poi':item['poi'],'pose':item['pose_id'],'frames':201,'bridge_frames':int(b.sum()),'bridge_fraction':float(b.mean()),'q90_worst_contact_A':float(np.quantile(np.maximum(v,k),.9)),'mean_VHL_contact_A':float(v.mean()),'mean_POI_contact_A':float(k.mean()),'longest_observed_bridge_span_ns':max((r['observed_span_ns'] for r in runs),default=0.),'bridge_runs':len(runs),'starts_bridged':bool(b[0]),'ends_bridged':bool(b[-1]),'source':source,'source_TPR_XTC_sha256':item['source_sha256']}
def extension(item):
 poi,pose=item['poi'],item['pose_id'];paths=item['source_sha256'];tpr=Path(next(x for x in paths if x.endswith('.tpr')));xtc=Path(next(x for x in paths if x.endswith('.xtc')))
 require(all(sha256(Path(p))==h for p,h in paths.items()),'MD source drift');evidence=json.loads((AUDIT/'prediction_compiled_roles_v3/per_pose_roles.json').read_text());row=next(r for r in evidence if r['poi']==poi and r['pose_id']==pose);require(row['tpr_sha256']==sha256(tpr),'Compiled identity evidence drift');expected={'SIK3':'Q9Y2K2','RPS6KB2':'Q9UBS0'}[poi]
 for uid,mol in [('P40337','Protein_chain_C'),(expected,'Protein_chain_D')]:require(any(r['uniprot']==uid and r['moltype']==mol for r in row['proteins']),'Wrong target identity')
 u=mda.Universe(str(tpr),str(xtc),format=ReadOnlyXTCReader);v=u.select_atoms('moltype Protein_chain_C');p=u.select_atoms('moltype Protein_chain_D');lig=u.select_atoms('resname UNL PTC');v=v[v.masses>1.5];p=p[p.masses>1.5];lig=lig[lig.masses>1.5];require(len(lig)==79 and len(u.select_atoms('moltype Protein_chain_C and name CA'))==140,'Incomplete ligand/VHL construct');rows=[]
 for ts in u.trajectory:
  t=float(ts.time)/1000
  if t>20+1e-6:break
  rows.append({'time_ns':t,'true_VHL_heavy_A':float(distance_array(lig.positions,v.positions,box=ts.dimensions).min()),'POI_heavy_A':float(distance_array(lig.positions,p.positions,box=ts.dimensions).min())})
 u.trajectory.close();frame=pd.DataFrame(rows);folder=ROOT/'results/md_frames'/poi/pose;folder.mkdir(parents=True,exist_ok=True);frame.to_csv(folder/'contacts.tsv',sep='\t',index=False);r=summarize(item,frame,'fresh_periodic_heavy_contact_replay');r['contact_table_sha256']=sha256(folder/'contacts.tsv');r['true_VHL_CA']=140;return r

def main():
 completion=json.loads((CAN/'analysis/dynamic/completion.json').read_text());manifest=pd.read_csv('/home/ayushregmibagale/protac-md/random_screen/manifest.tsv',sep='\t');metadata=json.loads((CAN/'config/role_reference.json').read_text());labels=json.loads((CAN/'config/current_md_labels.json').read_text())['labels'];records=completion['records'];require(len(records)==174,'Wrong completed roster');rows=[];ext=[]
 for item in records:
  if item['poi'] in ['SIK3','RPS6KB2']:ext.append(item);continue
  table=AUDIT/'contact_core_v1/poses'/item['poi']/item['pose_id']/'contacts_A.tsv';d=pd.read_csv(table,sep='\t');r=summarize(item,d,'previously_audited_periodic_core_contacts');r['contact_table_sha256']=sha256(table);r['true_VHL_CA']=140;rows.append(r)
 require(len(ext)==33 and len(rows)==141,'Incomplete case coverage')
 with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
  for r in pool.map(extension,ext):rows.append(r);print(r['poi'],r['pose'],r['bridge_fraction'],flush=True)
 for r in rows:
  target=r['poi']+'_full' if r['poi'] in ['SIK3','RPS6KB2'] else r['poi'];r['target']=target;r['uid']=metadata['targets'][target];r['label']=labels.get(target);r['label_status']='unknown' if target not in labels else 'confirmed_user';matched=manifest[(manifest.poi==r['poi'])&(manifest.pose_id==r['pose'])];require(len(matched)==1,'Missing source-pose lineage');r['source_pdb']=str(matched.source_pdb.iloc[0])
 atomic_json(ROOT/'validation/MD_case_receipts.json',{'status':'PASS','cases':174,'core_reused':141,'extension_recomputed':33,'frames':34974,'records':rows});table=pd.DataFrame([{k:v for k,v in r.items() if k!='source_TPR_XTC_sha256'} for r in rows]);table=table.sort_values(['bridge_fraction','q90_worst_contact_A','uid','pose'],ascending=[False,True,True,True]);table.to_csv(ROOT/'results/MD_ranked_poses.tsv',sep='\t',index=False);prospective=table[table.label_status=='unknown'].copy();prospective.to_csv(ROOT/'results/MD_prospective_ranked_poses.tsv',sep='\t',index=False);prospective.head(10).to_csv(ROOT/'results/top10_MD_prospective_poses.tsv',sep='\t',index=False);table.groupby(['target','uid','label_status'],as_index=False).agg(completed_poses=('pose','size'),mean_bridge_fraction=('bridge_fraction','mean'),maximum_bridge_fraction=('bridge_fraction','max')).to_csv(ROOT/'results/MD_target_summary.tsv',sep='\t',index=False)
 print('MD ranking complete;174 retained poses,33 prospective poses on2 kinases.')
if __name__=='__main__':main()
