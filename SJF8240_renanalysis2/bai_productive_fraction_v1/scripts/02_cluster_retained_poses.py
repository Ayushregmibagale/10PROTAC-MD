"""3A clustering of retained production-frame0 poses; not full Rosetta ensemble."""
from pathlib import Path
import sys,json,importlib.util,concurrent.futures
import numpy as np,pandas as pd
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import pdist
ROOT=Path(__file__).resolve().parents[1];SOURCE=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2/crl2_ubiquitination_v1');spec=importlib.util.spec_from_file_location('md_extract',SOURCE/'scripts/02_extract_md_compatibility.py');e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def coordinates(case):
 poi,pose=case;receipt=json.loads((SOURCE/'results/poses'/poi/pose/'receipt.json').read_text());tpr=e.OLD/'random_screen'/poi/pose/'prod/prod.tpr';xtc=e.OLD/'analysis/random_screen'/poi/pose/'prod_nojump.xtc'
 for p in [tpr,xtc]:e.require(e.sha256(p)==receipt['input_sha256'][str(p)],'Source changed')
 u=e.mda.Universe(str(tpr),str(xtc),format=e.ReadOnlyXTCReader);u.trajectory[0];e.require(abs(u.trajectory.ts.time)<.051,'Frame0 is not0ns');e.WholeProteinTree(u.select_atoms('protein'))(u.trajectory.ts);ligall=u.select_atoms('resname UNL');e.WholeProteinTree(ligall)(u.trajectory.ts);lig=ligall[ligall.masses>1.5];v=u.select_atoms('moltype Protein_chain_C');p=u.select_atoms('moltype Protein_chain_D');vh=v[v.masses>1.5];ph=p[p.masses>1.5];vc=v.select_atoms('name CA');pc=p.select_atoms('name CA');lc=e.ligand_connected_images(vh.positions,lig.positions,ph.positions,np.zeros(3),u.trajectory.ts.dimensions);result={'poi':poi,'pose':pose,'vhl_CA':vc.positions.astype(float).copy(),'poi_CA':pc.positions.astype(float).copy()+lc['poi_shift_A'],'vhl_sequence':''.join(e.seq1(r) for r in vc.resnames),'poi_sequence':''.join(e.seq1(r) for r in pc.resnames),'input_sha256':{str(tpr):e.sha256(tpr),str(xtc):receipt['input_sha256'][str(xtc)]}};u.trajectory.close();return result

def main():
 flags=pd.read_csv(ROOT/'results/frame0_pose_flags.tsv',sep='\t');cases=sorted(set(zip(flags.poi,flags.pose)));records=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
  for r in executor.map(coordinates,cases):records.append(r)
 membership=[];summaries=[];sources=[]
 for poi in sorted({r['poi'] for r in records}):
  subset=[r for r in records if r['poi']==poi];ref=subset[0];aligned=[]
  for r in subset:
   assert r['vhl_sequence']==ref['vhl_sequence'] and r['poi_sequence']==ref['poi_sequence'],'Unmatched coordinate correspondence';xyz,rmsd,rotation=e.fit_reference_points(r['vhl_CA'],ref['vhl_CA'],r['poi_CA']);aligned.append(xyz);sources.append({'poi':poi,'pose':r['pose'],'input_sha256':r['input_sha256']})
  arr=np.array(aligned);nca=arr.shape[1];dist=pdist(arr.reshape(len(arr),-1))/np.sqrt(nca);groups=fcluster(linkage(dist,method='complete'),t=3.,criterion='distance');square=np.zeros((len(arr),len(arr)));square[np.triu_indices(len(arr),1)]=dist;square+=square.T
  for cluster in sorted(set(groups)):
   indices=np.flatnonzero(groups==cluster);medoid=min(indices,key=lambda i:(float(square[i,indices].sum()),subset[i]['pose']));representative=subset[medoid]['pose']
   for i in indices:membership.append({'poi':poi,'pose':subset[i]['pose'],'cluster':int(cluster),'representative_pose':representative,'cluster_size':len(indices),'representative_method':'RMSD medoid; Rosetta interface energies unavailable for this retained MD ensemble'})
  for scope in ['closed_only','closed_plus_modeled_open']:
   f=flags[(flags.poi==poi)&(flags.scaffold_scope==scope)].set_index('pose');weighted=0;equal=[];naive=[]
   for cluster in sorted(set(groups)):
    indices=np.flatnonzero(groups==cluster);medoid=min(indices,key=lambda i:(float(square[i,indices].sum()),subset[i]['pose']));productive=bool(f.loc[subset[medoid]['pose'],'frame0_productive']);weighted+=len(indices)*productive;equal.append(productive)
   summaries.append({'poi':poi,'label':int(f.label.iloc[0]),'scaffold_scope':scope,'retained_poses':len(subset),'clusters_3A':len(equal),'productive_clusters':int(sum(equal)),'cluster_representative_productive_fraction':float(np.mean(equal)),'member_weighted_productive_fraction':weighted/len(subset),'unclustered_frame0_fraction':float(f.frame0_productive.mean())})
 pd.DataFrame(membership).to_csv(ROOT/'results/cluster_membership.tsv',sep='\t',index=False);summary=pd.DataFrame(summaries);summary.to_csv(ROOT/'results/cluster_target_scores.tsv',sep='\t',index=False);(ROOT/'validation/frame0_coordinate_sources.json').write_text(json.dumps(sources,indent=2)+'\n');from sklearn.metrics import roc_auc_score
 metrics=[]
 for scope,g in summary.groupby('scaffold_scope'):
  for score in ['cluster_representative_productive_fraction','member_weighted_productive_fraction','unclustered_frame0_fraction']:metrics.append({'scaffold_scope':scope,'score':score,'AUROC':float(roc_auc_score(g.label,g[score])),'independent_targets':8})
 pd.DataFrame(metrics).to_csv(ROOT/'results/cluster_metrics.tsv',sep='\t',index=False);(ROOT/'validation/clustering_method.json').write_text(json.dumps({'RMSD_threshold_A':3.,'alignment':'whole trueVHL C-alpha fit; matching target C-alpha sequence/order required','clustering':'complete-linkage,max pairwise RMSD within cluster<=3A','representative':'coordinate medoid,not minimum Rosetta interface energy','model_scope':'141 retained production-frame0 models,not originalfull Rosetta decoy pools','members_weighted':'Classify representative; propagate to original retained cluster members','not_literal_MMTSB_reproduction':True,'cases':len(cases),'script_sha256':e.sha256(__file__)},indent=2)+'\n');print(summary.to_string(index=False));print(pd.DataFrame(metrics).to_string(index=False))
if __name__=='__main__':main()
