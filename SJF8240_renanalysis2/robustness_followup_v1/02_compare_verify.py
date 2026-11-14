"""Fresh-instance CPU replay, fixed baseline comparison and held-out metrics."""
import importlib.util,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('runner',ROOT/'01_run.py');r=importlib.util.module_from_spec(sp);sp.loader.exec_module(r)
import numpy as np,pandas as pd,torch
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score

def main():
 r.freeze();torch.set_num_threads(1);rows=json.loads((ROOT/'completion.json').read_text());pred=[];replays=[]
 for row in rows:
  seed,v,pair=row['base_seed'],row['variant'],row['pair'];path=ROOT/'folds'/str(seed)/v/(pair+'.pt');assert r.sha(path)==row['weight_sha256'];m=r.TemporalEncoder();m.load_state_dict(torch.load(path,map_location='cpu',weights_only=True)['weights']);scores=r.score(m,r.graphs(v));err=max(abs(scores[t]-row['scores'][t]) for t in scores);assert err<=1e-7;replays.append({'seed':seed,'variant':v,'pair':pair,'maximum_error':err,'checkpoint_sha256':r.sha(path)})
  other=next(x for x in rows if x['base_seed']==seed and x['pair']==pair and x['variant']!=v);assert row['initial_digest']==other['initial_digest']
  pred.append({'model':v,'seed':seed,'pair':pair,'positive':row['positive'],'negative':row['negative'],'positive_score':scores[row['positive']],'negative_score':scores[row['negative']],'threshold':row['training_only_threshold']})
 # Five catalytic node features, pooled over POI lysines and time windows; fixed, label-independent definition.
 X={};geometry={}
 for t in r.P['cohort']:
  z=np.load(r.OLD/'inputs/native_features'/(t+'.npz'));a=np.stack([z[f'x{w}'][z[f'x{w}'][:,23]>0,29:34].mean(0) for w in range(4)]);X[t]=np.r_[a.mean(0),a.std(0),a[-1]-a[0]];geometry[t]=float(a[:,4].mean())
 pd.DataFrame([{'target':t,**{f'feature{i}':x for i,x in enumerate(X[t])}} for t in X]).to_csv(ROOT/'baseline_features.tsv',sep='\t',index=False)
 for f in r.P['folds']:
  tr=f['training'];assert not set(tr)&{f['held_out_positive'],f['held_out_negative']};model=make_pipeline(StandardScaler(),LogisticRegression(C=1.,solver='lbfgs',max_iter=1000,random_state=42));model.fit(np.stack([X[t] for t in tr]),[r.P['labels'][t] for t in tr]);scores=dict(zip(r.P['cohort'],model.predict_proba(np.stack([X[t] for t in r.P['cohort']]))[:,1]));np.savez(ROOT/('baseline_'+f['pair']+'.npz'),coef=model[-1].coef_,intercept=model[-1].intercept_,mean=model[0].mean_,scale=model[0].scale_)
  for name,sc in [('logistic_catalytic',scores),('geometry_joint_fraction',geometry)]:
   pred.append({'model':name,'seed':42,'pair':f['pair'],'positive':f['held_out_positive'],'negative':f['held_out_negative'],'positive_score':sc[f['held_out_positive']],'negative_score':sc[f['held_out_negative']],'threshold':r.threshold(sc,tr)})
 df=pd.DataFrame(pred);df.to_csv(ROOT/'heldout_predictions.tsv',sep='\t',index=False);summary=[]
 for (model,seed),g in df.groupby(['model','seed']):
  for policy in ['fixed_0.5','training_only']:
   cuts=.5 if policy=='fixed_0.5' else g.threshold;sen=float((g.positive_score>=cuts).mean());spe=float((g.negative_score<cuts).mean());summary.append({'model':model,'seed':int(seed),'policy':policy,'matched_pair_AUROC':float(((g.positive_score>g.negative_score)+.5*(g.positive_score==g.negative_score)).mean()),'sensitivity':sen,'specificity':spe,'balanced_accuracy':(sen+spe)/2})
 metrics=pd.DataFrame(summary);metrics.to_csv(ROOT/'metrics.tsv',sep='\t',index=False);r.save(ROOT/'CPU_VERIFICATION.json',{'status':'PASS','models':len(replays),'scores_replayed':len(replays)*8,'records':replays})
 # Unknown inputs are not interchangeable with temporal catalytic graphs: availability reported separately.
 lines=['# Five-seed robustness, baselines and threshold assessment','', 'EPHA2 remains label 1. SIK3 and RPS6KB2 are unknown and excluded from fitting/evaluation. Historical chemistry retained; no supplier equivalence claimed.','', '## Results', '',metrics.to_markdown(index=False),'','## Interpretation and method','', 'Five predeclared seeds: 42, 7, 123, 2026, 31415. Same 12 outer pair splits, 300 epochs, unweighted BCE, architecture and optimizer as prior ablation. No seed, hyperparameter or feature selection based on outcomes. All checkpoints saved and independently replayed on CPU.','', 'Logistic regression: 15 fixed catalytic summaries (five lysine means across windows, temporal standard deviations and end-minus-start); StandardScaler fitted only on outer training targets, L2 C=1. Geometry score is mean simultaneous exposed/reachable/scaffold-compatible/bridging fraction (feature 5); higher is fixed to degrader-like. The geometry score is not a calibrated probability, so its 0.5 threshold is only a numerical comparator.','', 'Threshold maximizes balanced accuracy on the six training proteins; ties choose closest to 0.5 then smaller threshold. Outer held-out labels never select thresholds. This is training-resubstitution calibration, not nested CV: outer training has only one negative, so a non-degenerate inner negative holdout is unavailable. Training fit can be optimistic and score calibration may shift to unseen proteins.','', 'Eight independent labelled proteins and two negatives; pairs, seeds and repeated predictions are correlated, not new biological replicates. Report seed variability descriptively without a naive significance test.']
 (ROOT/'FINDINGS.md').write_text('\n'.join(lines)+'\n');print(metrics.to_string(index=False));print('CPU replay PASS',len(replays))
if __name__=='__main__':main()
