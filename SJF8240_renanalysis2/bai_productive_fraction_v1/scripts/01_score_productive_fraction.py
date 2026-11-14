"""Bai-style16A exposed-Lys/clash rule applied to retained VHL ensembles.
Not a literal reproduction of CRBN/CRL4A or the complete Rosetta ensemble.
"""
from pathlib import Path
import json,hashlib,itertools
import pandas as pd,numpy as np
from sklearn.metrics import roc_auc_score,balanced_accuracy_score,recall_score
ROOT=Path(__file__).resolve().parents[1];SOURCE=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2/crl2_ubiquitination_v1')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 protocol=json.loads((ROOT/'config/protocol.json').read_text());complete=json.loads((SOURCE/'completion.json').read_text());assert complete['cases']==141 and complete['frames']==28341;labels=json.loads((SOURCE/'config/labels.json').read_text())['labels'];rows=[];poseflags=[]
 for case in complete['records']:
  path=SOURCE/'results/poses'/case['poi']/case['pose']/'packet.npz';assert sha(path)==case['packet_sha256'];p=np.load(path);names=p['spec_names'].tolist();closed=names.index('closed_observed__VHL');opened=names.index('open_RBX1_transferred_module__VHL');ex=p['nz_exposed_fraction'];dist=p['donor_distance_A'];clash=p['target_scaffold_max_overlap_A'];fit=p['anchor_RMSD_A'];
  for tail_index,tail in [(0,'tail_fit74'),(1,'tail_fit75')]:
   for exposure in [.05,.1,.2]:
    for cutoff in [16.,20.]:
     # <=16A follows paper's stated inclusive boundary; previous20A implementation used<.
     productive=((dist[:,:,tail_index,:]<=cutoff)&(ex[None,:,:]>=exposure)).any(axis=2)&(clash<=.8);never_noclash=clash<=.8
     for scope,indices in [('closed_only',[closed]),('closed_plus_modeled_open',[closed,opened])]:
      binary=productive[indices].any(axis=0);quality=(fit[indices]<=3.).any(axis=0);strict=(productive[indices]&(fit[indices]<=3.)).any(axis=0)
      rows.append({'poi':case['poi'],'pose':case['pose'],'label':labels[case['poi']],'tail':tail,'exposure_fraction':exposure,'distance_A':cutoff,'scaffold_scope':scope,'productive_frame_fraction':float(binary.mean()),'productive_frame_fraction_fit3_sensitivity':float(strict.mean()),'anchor_fit_any_fraction':float(quality.mean()),'frame0_productive':bool(binary[0]),'frame0_fit3_productive':bool(strict[0]),'bridge_fraction_reference_only':float(p['bridge'].mean()),'orientation_filter_applied':False,'bridge_filter_applied':False})
      if tail=='tail_fit75' and exposure==.1 and cutoff==16.:
       poseflags.append({'poi':case['poi'],'pose':case['pose'],'label':labels[case['poi']],'scaffold_scope':scope,'frame0_productive':bool(binary[0]),'frame0_productive_fit3':bool(strict[0]),'MD_productive_fraction':float(binary.mean()),'MD_productive_fraction_fit3':float(strict.mean())})
 d=pd.DataFrame(rows);flags=pd.DataFrame(poseflags);target=d.groupby(['poi','label','tail','exposure_fraction','distance_A','scaffold_scope'],as_index=False).agg(productive_frame_fraction=('productive_frame_fraction','mean'),productive_frame_fraction_fit3_sensitivity=('productive_frame_fraction_fit3_sensitivity','mean'),frame0_productive_fraction=('frame0_productive','mean'),frame0_productive_fraction_fit3=('frame0_fit3_productive','mean'),poses=('pose','nunique'));primary=target[(target['tail']=='tail_fit75')&(target.exposure_fraction==.1)&(target.distance_A==16)&(target.scaffold_scope=='closed_plus_modeled_open')].copy();assert len(primary)==8;metrics=[]
 for keys,g in target.groupby(['tail','exposure_fraction','distance_A','scaffold_scope']):
  for score in ['productive_frame_fraction','productive_frame_fraction_fit3_sensitivity','frame0_productive_fraction','frame0_productive_fraction_fit3']:
   y=g.label.to_numpy();v=g[score].to_numpy();auc=float(roc_auc_score(y,v));pred=v>0;permutations=[]
   for inds in itertools.combinations(range(8),2):
    yy=np.ones(8,int);yy[list(inds)]=0;permutations.append(roc_auc_score(yy,v))
   metrics.append(dict(zip(['tail','exposure_fraction','distance_A','scaffold_scope'],keys))|{'score':score,'AUROC':auc,'balanced_accuracy_at_nonzero':float(balanced_accuracy_score(y,pred)),'sensitivity_at_nonzero':float(recall_score(y,pred,zero_division=0)),'specificity_at_nonzero':float(recall_score(1-y,1-pred,zero_division=0)),'exact_two_sided_p':float(np.mean(abs(np.array(permutations)-.5)>=abs(auc-.5)-1e-12)),'classifier_threshold_validated':False})
 for name,frame in [('pose_scores.tsv',d),('frame0_pose_flags.tsv',flags),('target_scores_all_variants.tsv',target),('primary_target_scores.tsv',primary),('metrics.tsv',pd.DataFrame(metrics))]:frame.to_csv(ROOT/'results'/name,sep='\t',index=False)
 receipt={'source_cohort_sha256':sha(SOURCE/'completion.json'),'source_extraction_manifest_sha256':sha(SOURCE/'config/execution_manifest.json'),'cases':141,'frames':28341,'targets':8,'script_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'config/protocol.json'),'outputs_sha256':{p.name:sha(p) for p in (ROOT/'results').glob('*.tsv')},'literal_paper_reproduction':False};(ROOT/'validation/scoring_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(primary[['poi','label','productive_frame_fraction','frame0_productive_fraction']].to_string(index=False));print(pd.DataFrame(metrics).query("tail=='tail_fit75' and exposure_fraction==0.1 and distance_A==16 and scaffold_scope=='closed_plus_modeled_open'").to_string(index=False))
if __name__=='__main__':main()
