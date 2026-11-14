"""Frozen three-reference retrospective target comparison; no classifier fitting."""
from pathlib import Path
import json,hashlib,itertools
import numpy as np,pandas as pd
from sklearn.metrics import roc_auc_score,balanced_accuracy_score,recall_score
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 complete=json.loads((ROOT/'completion.json').read_text());assert complete['cases']==141 and complete['frames']==28341;labels=json.loads((ROOT/'config/labels.json').read_text())['labels'];rows=[];sens=[]
 for rr in complete['records']:
  path=ROOT/'results/poses'/rr['poi']/rr['pose']/'packet.npz';assert sha(path)==rr['packet_sha256'];loaded=np.load(path);p={key:loaded[key] for key in loaded.files};loaded.close();ex=p['nz_exposed_fraction'];bridge=p['bridge'].astype(bool)
  for si,name in enumerate(p['spec_names']):
   ref,anchor=str(name).rsplit('__',1);dist=p['donor_distance_A'][si,:,0,:];cos=p['sidechain_outward_cosine'][si,:,0,:];clear=p['target_scaffold_max_overlap_A'][si]<=.8;fit=p['anchor_RMSD_A'][si]<=3.;near=dist<=20;exposed=near&(ex>=.1);prod=exposed.any(axis=1)&clear;oriented=(exposed&(cos>=0)).any(axis=1)&clear;row={'poi':rr['poi'],'pose':rr['pose'],'label':labels[rr['poi']],'reference':ref,'anchor':anchor,'bridge_fraction':float(bridge.mean()),'distance_only_fraction':float(near.any(axis=1).mean()),'scaffold_compatible_fraction':float(clear.mean()),'anchor_fit3_fraction':float(fit.mean()),'productive_fraction':float(prod.mean()),'joint_bridge_productive_fraction':float((prod&bridge).mean()),'productive_conditional_on_bridge':float(prod[bridge].mean()) if bridge.any() else np.nan,'oriented_productive_fraction':float(oriented.mean()),'fit3_productive_fraction':float((prod&fit).mean()),'maximum_anchor_RMSD_A':float(p['anchor_RMSD_A'][si].max())};rows.append(row)
   for cutoff,exposure,angle in itertools.product([12,16,20,24,30],[.05,.1,.2],['none','outward']):
    good=(dist<=cutoff)&(ex>=exposure)
    if angle=='outward':good&=cos>=0
    score=good.any(axis=1)&clear;sens.append({'poi':rr['poi'],'reference':ref,'anchor':anchor,'distance_A':cutoff,'exposure_fraction':exposure,'orientation':angle,'score':float(score.mean())})
 pose=pd.DataFrame(rows);cols=[x for x in pose.select_dtypes('number').columns if x!='label'];target=pose.groupby(['poi','label','reference','anchor'],as_index=False)[cols].mean();target['poses']=target.poi.map(pose.groupby('poi').pose.nunique());primary=target[(target.reference=='8RX0_relaxed_G76')&(target.anchor=='VHL')];assert len(primary)==8;metrics=[]
 for (ref,anchor),g in target.groupby(['reference','anchor']):
  for score in ['bridge_fraction','distance_only_fraction','scaffold_compatible_fraction','productive_fraction','joint_bridge_productive_fraction','productive_conditional_on_bridge','oriented_productive_fraction','fit3_productive_fraction']:
   v=g[score].to_numpy();y=g.label.to_numpy();good=np.isfinite(v);auc=float(roc_auc_score(y[good],v[good]));pred=v>0;null=[]
   if good.all():
    for neg in itertools.combinations(range(8),2):
     yy=np.ones(8,int);yy[list(neg)]=0;null.append(roc_auc_score(yy,v))
   metrics.append({'reference':ref,'anchor':anchor,'score':score,'AUROC':auc,'balanced_accuracy_nonzero':float(balanced_accuracy_score(y[good],pred[good])),'sensitivity_nonzero':float(recall_score(y[good],pred[good],zero_division=0)),'specificity_nonzero':float(recall_score(1-y[good],1-pred[good],zero_division=0)),'exact_two_sided_AUROC_p':float(np.mean(abs(np.array(null)-.5)>=abs(auc-.5)-1e-12)) if null else np.nan,'trained_classifier':False})
 st=pd.DataFrame(sens).groupby(['poi','reference','anchor','distance_A','exposure_fraction','orientation'],as_index=False).score.mean();st['label']=st.poi.map(labels);sm=[]
 for key,g in st.groupby(['reference','anchor','distance_A','exposure_fraction','orientation']):sm.append(dict(zip(['reference','anchor','distance_A','exposure_fraction','orientation'],key))|{'AUROC':float(roc_auc_score(g.label,g.score))})
 for name,df in [('pose_native_scores.tsv',pose),('target_native_scores.tsv',target),('primary_target_scores.tsv',primary),('native_metrics.tsv',pd.DataFrame(metrics)),('sensitivity_target_scores.tsv',st),('sensitivity_metrics.tsv',pd.DataFrame(sm))]:df.to_csv(ROOT/'results'/name,sep='\t',index=False)
 receipt={'cases':141,'frames':28341,'targets':8,'secondary':'chemically modeledG76,not replacement ofparent primary;20A exposed/no clash','script_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'config/protocol.json'),'outputs_sha256':{p.name:sha(p) for p in (ROOT/'results').glob('*.tsv')},'no_training':True};(ROOT/'validation/comparison_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(primary[['poi','label','productive_fraction','joint_bridge_productive_fraction','scaffold_compatible_fraction']].to_string(index=False));print(pd.DataFrame(metrics).query("reference=='8RX0_relaxed_G76' and anchor=='VHL'").to_string(index=False))
if __name__=='__main__':main()
