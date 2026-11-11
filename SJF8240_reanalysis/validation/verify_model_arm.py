"""Read-only artifact audit and independent CPU inference replay; no fitting."""
from pathlib import Path
import argparse,json,hashlib,sys,importlib.util,datetime
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,accuracy_score,confusion_matrix
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('frozen_training',ROOT/'scripts/04_train_dual_lopo.py')
training=importlib.util.module_from_spec(spec);spec.loader.exec_module(training)
torch=training.torch;torch.set_num_threads(1)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def require(condition,message):
 if not condition:raise ValueError(message)
def main(role):
 protocol=json.loads((ROOT/'config/protocol.json').read_text());cohort=sorted(protocol['primary_cohort']);labels=protocol['historical_labels']
 run=ROOT/'analysis/gnn'/role/'seed_42';completion=json.loads((run/'completion.json').read_text())
 require(completion['status']=='NEW_EXPLORATORY_LOPO_COMPLETE' and completion['cases']==16 and completion['role']==role,'Incomplete model arm')
 settings=completion['settings'];require(settings['cohort']==cohort and settings['labels']==labels,'Frozen cohort/labels mismatch')
 for path,digest in settings['inputs_sha256'].items():require(sha(path)==digest,'Input/source drift: '+path)
 for name,digest in completion['output_sha256'].items():
  require(Path(name).name==name and sha(run/name)==digest,'Model output drift: '+name)
 df=pd.read_csv(run/'predictions.tsv',sep='\t');extra=pd.read_csv(run/'prediction_only_fold_scores.tsv',sep='\t')
 require(len(df)==16 and not df.poi.duplicated().any() and set(df.poi)==set(cohort),'Held-out cohort not unique/complete')
 require(len(extra)==32 and set(extra.poi)==set(protocol['prediction_only']) and not extra.duplicated(['held_out','poi']).any(),'Prediction-only fold coverage mismatch')
 graphs={};dynamic={}
 for poi in cohort+protocol['prediction_only']:
  label=labels.get(poi,-1);graphs[poi]=training.load_graph(ROOT/'analysis/static'/poi/role/'graph.npz',poi,label,9)
  path=ROOT/'analysis/dynamic/targets'/poi/role/'graph.npz'
  if path.exists():dynamic[poi]=training.load_graph(path,poi,label,3)
 replay=[];max_error=0.
 for poi in cohort:
  receipt=json.loads((run/f'fold_{poi}.json').read_text());checkpoint=run/f'fold_{poi}.pt'
  require(sha(checkpoint)==receipt['checkpoint_sha256'],'Weight hash mismatch')
  state=torch.load(checkpoint,map_location='cpu',weights_only=True)
  require(state['held_out']==poi and state['training_targets']==[p for p in cohort if p!=poi],'Fold leakage or target mismatch')
  require(state['settings']==settings and state['seed']==training.fold_seed(42,poi),'Fold settings/seed drift')
  row=df.loc[df.poi==poi].iloc[0]
  require(row.label==labels[poi] and bool(row.has_md)==(poi in dynamic),'Prediction label/MD availability mismatch')
  require(int(row.pred)==int(row.score>=.5) and np.isfinite(row.score) and 0<=row.score<=1,'Invalid prediction score/class')
  require(receipt['prediction']==state['prediction'],'Checkpoint versus receipt prediction mismatch')
  require(abs(row.score-receipt['prediction']['score'])<=1e-12,'Prediction table versus fold mismatch')
  require(abs(receipt['replayed_score']-row.score)<=1e-7,'Original saved-GPU-state replay failed')
  model=training.DualEncoderGNN(hidden=64,heads=4,dropout=.3);model.load_state_dict(state['model_state_dict'])
  score=training.predict(model,graphs[poi],dynamic.get(poi));error=abs(score-row.score);max_error=max(max_error,error)
  require(error<=1e-5,'Independent CPU replay error exceeds explicitly separate 1e-5 tolerance')
  for target in protocol['prediction_only']:
   erow=extra.loc[(extra.held_out==poi)&(extra.poi==target)].iloc[0]
   predicted=training.predict(model,graphs[target],dynamic.get(target));difference=abs(predicted-erow.score);max_error=max(max_error,difference)
   require(difference<=1e-5 and pd.isna(erow.experimental_label),'Prediction-only replay/label mismatch')
  replay.append({'held_out':poi,'training_targets':state['training_targets'],'saved_GPU_score':float(row.score),'independent_CPU_score':score,'absolute_difference':error,'checkpoint_sha256':sha(checkpoint)})
 metrics={'AUROC':float(roc_auc_score(df.label,df.score)),'accuracy':float(accuracy_score(df.label,df.pred)),'confusion_matrix_labels_0_1':confusion_matrix(df.label,df.pred,labels=[0,1]).tolist()}
 for key in ['AUROC','accuracy']:require(abs(metrics[key]-completion[key])<=1e-12,'Reported '+key+' mismatch')
 subgroups=training.subgroup_aurocs(df.label,df.score,{'md':df.has_md.to_numpy(dtype=bool),'nonmd':~df.has_md.to_numpy(dtype=bool)})
 require(subgroups==completion['subgroups'],'Subgroup metric mismatch')
 report={'status':'INDEPENDENT_MODEL_ARTIFACT_AND_CPU_REPLAY_PASS','role':role,'checked_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'folds':16,'prediction_only_inferences':32,'source_inputs_hashed':len(settings['inputs_sha256']),'completion_sha256':sha(run/'completion.json'),'validator_sha256':sha(__file__),'CPU_torch':str(torch.__version__),'CPU_replay_absolute_tolerance':1e-5,'original_GPU_replay_tolerance':1e-7,'maximum_CPU_score_difference_including_prediction_only':max_error,'metrics':metrics,'subgroups':subgroups,'per_fold_replay':replay,'new_training_executed':False,'independent_biological_validation':False}
 path=ROOT/'validation'/f'verification_{role}.json'
 require(not path.exists(),'Preserve existing validation receipt; do not overwrite')
 path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='per_fold_replay'},indent=2))
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--role',choices=['elob','vhl'],required=True);main(parser.parse_args().role)
