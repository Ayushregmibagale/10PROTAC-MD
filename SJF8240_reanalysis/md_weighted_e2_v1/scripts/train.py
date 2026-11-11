"""Frozen target-LOPO MD-only GNN; corrected user labels, no new simulations."""
from pathlib import Path
import sys,argparse,json,random,time,os,fcntl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from common import *
import torch.nn.functional as F

def execute(role,device):
 frozen=verify_frozen();p,graphs=load_inputs(role,device);torch.set_num_threads(1);torch.use_deterministic_algorithms(True)
 require(device!='cuda' or torch.cuda.is_available(),'Requested GPU unavailable')
 run=ROOT/'analysis'/role;run.mkdir(parents=True,exist_ok=True);lock=(run/'training.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 settings={'role':role,'device':device,'protocol_sha256':sha(ROOT/'config/protocol.json'),'execution_manifest_sha256':frozen,'input_manifest_sha256':sha(ROOT/'config/input_manifest.json'),'label_map_sha256':sha(ROOT/'config/user_labels.json'),'torch':str(torch.__version__),'cohort':p['cohort'],'labels':p['labels'],'epochs':p['epochs'],'seed':p['seed'],'architecture':'MD_GAT_backbone_plus_two_graph_E2_scalars','loss':p['loss'],'all_evaluation_targets_have_MD':True,'deterministic_algorithms':True,'disabled_optional_backends':DISABLED}
 intent=run/'intent.json'
 if intent.exists():require(json.loads(intent.read_text())==settings,'Run intent/source drift')
 else:atomic_json(intent,settings)
 rows=[];external=[]
 for held in p['cohort']:
  checkpoint=run/f'fold_{held}.pt';receipt=run/f'fold_{held}.json';training=[x for x in p['cohort'] if x!=held];seed=fold_seed(p['seed'],held);class_weights=balanced_class_weights(p['labels'],training);e2_scaler=apply_e2_features(graphs,p,training)
  if receipt.exists():
   record=json.loads(receipt.read_text());require(record['status']=='FOLD_COMPLETE' and sha(checkpoint)==record['checkpoint_sha256'],'Cached fold receipt/weight drift')
   state=torch.load(checkpoint,map_location=device,weights_only=True);require(state['settings']==settings and state['held_out']==held and state['training_targets']==training and state['seed']==seed and state['e2_normalization']==e2_scaler,'Cached fold leakage/settings mismatch')
   model=MDEncoder().to(device);model.load_state_dict(state['model_state_dict']);require(state['prediction']==record['prediction'] and abs(predict(model,graphs[held])-record['prediction']['score'])<=1e-7,'Cached fold replay mismatch')
   require(record['prediction']['label']==p['labels'][held] and record['prediction']['poi']==held,'Cached fold label/identity mismatch')
   for item in record['prediction_only']:
    require(item['poi'] in p['prediction_only'] and item['held_out']==held and abs(predict(model,graphs[item['poi']])-item['score'])<=1e-7,'Cached prediction-only replay mismatch')
   rows.append(record['prediction']);external.extend(record['prediction_only']);continue
  random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
  if device=='cuda':torch.cuda.manual_seed_all(seed)
  model=MDEncoder().to(device);opt=torch.optim.Adam(model.parameters(),lr=p['learning_rate'],weight_decay=p['weight_decay']);scheduler=torch.optim.lr_scheduler.CosineAnnealingLR(opt,p['epochs']);order=list(training);start=time.monotonic()
  require(held not in training and len(set(training))==7,'Held-out leakage')
  for epoch in range(p['epochs']):
   np.random.shuffle(order)
   for poi in order:
    model.train();opt.zero_grad();loss=balanced_bce(model(graphs[poi]),p['labels'][poi],class_weights);require(bool(torch.isfinite(loss)),'Nonfinite training loss');loss.backward();opt.step()
   scheduler.step()
   if (epoch+1)%50==0:print(json.dumps({'role':role,'held_out':held,'epoch':epoch+1,'epochs':p['epochs'],'elapsed_s':round(time.monotonic()-start,2)}),flush=True)
  score=predict(model,graphs[held]);prediction={'poi':held,'label':p['labels'][held],'score':score,'pred':int(score>=p['decision_threshold']),'fold_seed':seed,'has_MD':True}
  extra=[{'held_out':held,'poi':poi,'score':predict(model,graphs[poi]),'experimental_label':None} for poi in p['prediction_only']]
  state={'model_state_dict':{k:v.detach().cpu().clone() for k,v in model.state_dict().items()},'settings':settings,'held_out':held,'training_targets':training,'class_weights':{str(k):v for k,v in class_weights.items()},'e2_normalization':e2_scaler,'seed':seed,'prediction':prediction}
  if checkpoint.exists():checkpoint.rename(checkpoint.with_suffix('.unreceipted_'+str(time.time_ns())+'.pt'))
  temp=checkpoint.with_suffix('.pending.pt');torch.save(state,temp);os.replace(temp,checkpoint)
  replay=MDEncoder().to(device);replay.load_state_dict(torch.load(checkpoint,map_location=device,weights_only=True)['model_state_dict']);replayed=predict(replay,graphs[held]);require(abs(score-replayed)<=1e-7,'Saved-weight replay mismatch')
  for item in extra:require(abs(predict(replay,graphs[item['poi']])-item['score'])<=1e-7,'Prediction-only saved-weight replay mismatch')
  record={'status':'FOLD_COMPLETE','prediction':prediction,'prediction_only':extra,'training_targets':training,'class_weights':{str(k):v for k,v in class_weights.items()},'e2_normalization':e2_scaler,'training_class_counts':{'positive':sum(p['labels'][x] for x in training),'negative':sum(1-p['labels'][x] for x in training)},'checkpoint_sha256':sha(checkpoint),'replayed_score':replayed,'elapsed_s':time.monotonic()-start}
  atomic_json(receipt,record);rows.append(prediction);external.extend(extra);print(json.dumps(record),flush=True)
 require(verify_frozen()==frozen,'Source changed during training');load_inputs(role,'cpu')
 import pandas as pd
 pd.DataFrame(rows).to_csv(run/'predictions.tsv',sep='\t',index=False);pd.DataFrame(external,columns=['held_out','poi','score','experimental_label']).to_csv(run/'prediction_only_fold_scores.tsv',sep='\t',index=False)
 names=['intent.json','predictions.tsv','prediction_only_fold_scores.tsv']+[f'fold_{poi}.{suffix}' for poi in p['cohort'] for suffix in ['pt','json']]
 summary={'status':'CORRECTED_MD_ONLY_LOPO_COMPLETE','role':role,'cases':8,'settings':settings,'metrics':metrics(rows),'all_evaluation_targets_have_MD':True,'EPHA2_degraded':True,'E2_features_used':True,'static_graphs_used':False,'new_MD_or_docking':False,'independent_biological_validation':False,'outputs_sha256':{name:sha(run/name) for name in names}}
 atomic_json(run/'completion.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--role',choices=['elob','vhl'],required=True);parser.add_argument('--device',choices=['cpu','cuda'],default='cuda');args=parser.parse_args();execute(args.role,args.device)
