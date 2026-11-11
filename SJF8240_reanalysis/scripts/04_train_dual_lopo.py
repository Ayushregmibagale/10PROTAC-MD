"""New paired 04d raw-graph LOPO; no tabular ML features or old weights."""
from pathlib import Path
import argparse,json,sys,hashlib,random,time,os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts/lib'))
from gnn_contract import configure_native_pyg,fold_seed,validate_feature_graph,subgroup_aurocs
DISABLED=configure_native_pyg()
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from dual_encoder_model import DualEncoderGNN
from sklearn.metrics import roc_auc_score,accuracy_score
from protac_qc import ValidationError,atomic_json
import pandas as pd
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_graph(path,poi,label,width):
 with np.load(path,allow_pickle=False) as d:
  if str(d['poi'][0])!=poi:raise ValidationError('Graph target mismatch')
  if d['label'].shape!=(1,) or not np.isfinite(d['label'][0]) or d['label'][0]!=label:raise ValidationError('Graph label mismatch')
  if d['edge_index'].dtype.kind not in 'iu':raise ValidationError('Raw edge indices must be integer')
  g=Data(x=torch.tensor(d['x'],dtype=torch.float32),edge_index=torch.tensor(d['edge_index'],dtype=torch.long),edge_attr=torch.tensor(d['edge_attr'],dtype=torch.float32),y=torch.tensor([label]),poi=poi,feature_schema='SJF8240-retained-paired-v1')
 validate_feature_graph(g,width,schema='SJF8240-retained-paired-v1');return g
@torch.no_grad()
def predict(model,sg,dg):
 model.eval();return float(torch.sigmoid(model(sg,dg)).item())
def run(args):
 if args.probe:
  print(json.dumps({'torch':torch.__version__,'cuda_available':torch.cuda.is_available(),'gpu':torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,'cuda_arch_list':torch.cuda.get_arch_list()}));return
 cfg=json.loads((ROOT/'config/protocol.json').read_text())
 if args.epochs!=cfg['training']['epochs'] or args.seed!=cfg['training']['seed']:raise ValidationError('Training budget or seed differs from frozen protocol')
 records=json.loads((ROOT/'config/source_inventory.json').read_text())['records'];cohort=sorted(r['poi'] for r in records if r['training_inclusion']);labels={r['poi']:r['historical_label'] for r in records if r['training_inclusion']}
 if cohort!=sorted(cfg['primary_cohort']) or labels!=cfg['historical_labels']:raise ValidationError('Frozen cohort/label drift')
 if len(cohort)!=16 or len(set(labels.values()))!=2 or 'CDK17' in cohort:raise ValidationError('Frozen identity-matched sixteen-target cohort required')
 if not (ROOT/'analysis/dynamic/completion.json').exists():raise ValidationError('Dynamic completion required')
 if any(not (ROOT/'analysis/static'/r['poi']/'receipt.json').exists() for r in records):raise ValidationError('All nineteen static targets must finish before training')
 # Completion markers are not sufficient: bind the actual graphs and provenance.
 for r in records:
  folder=ROOT/'analysis/static'/r['poi'];receipt=json.loads((folder/'receipt.json').read_text())
  if receipt.get('status')!='STATIC_PAIRED_COMPLETE' or receipt.get('poi')!=r['poi'] or sha(folder/'pose_provenance.json')!=receipt['pose_provenance_sha256']:raise ValidationError('Static receipt/provenance mismatch')
  if sha(ROOT/r['pose_list'])!=r['pose_list_sha256']:raise ValidationError('Pose list drift')
  for role in ['elob','vhl']:
   if receipt['outputs'][role]['poses']!=r['poses'] or sha(folder/role/'graph.npz')!=receipt['outputs'][role]['graph_sha256']:raise ValidationError('Static graph hash/count drift')
 dynamic=json.loads((ROOT/'analysis/dynamic/completion.json').read_text())
 if dynamic.get('status')!='DYNAMIC_COMPLETE' or dynamic.get('cases')!=174 or len(dynamic['records'])!=174:raise ValidationError('Dynamic completion mismatch')
 device=args.device
 if device=='cuda' and not torch.cuda.is_available():raise ValidationError('Requested CUDA unavailable; use explicit CPU or authorized host execution')
 torch.set_num_threads(1);torch.use_deterministic_algorithms(True)
 inputs={str(ROOT/'config/protocol.json'):sha(ROOT/'config/protocol.json'),str(ROOT/'scripts/lib/dual_encoder_model.py'):sha(ROOT/'scripts/lib/dual_encoder_model.py'),str(Path(__file__)):sha(__file__)};sg={};dg={}
 for relative in ['config/source_inventory.json','config/role_reference.json','scripts/lib/gnn_contract.py','scripts/lib/protac_qc.py','analysis/dynamic/completion.json']:
  path=ROOT/relative;inputs[str(path)]=sha(path)
 for r in records:
  path=ROOT/'analysis/static'/r['poi']/'receipt.json';inputs[str(path)]=sha(path)
 prediction_targets=['SIK3_full','RPS6KB2_full']
 for poi in cohort+prediction_targets:
  label=labels.get(poi,-1);path=ROOT/'analysis/static'/poi/args.role/'graph.npz';inputs[str(path)]=sha(path);sg[poi]=load_graph(path,poi,label,9).to(device)
  path=ROOT/'analysis/dynamic/targets'/poi/args.role/'graph.npz'
  if path.exists():
   receipt=path.parent/'receipt.json';metadata=json.loads(receipt.read_text())
   if metadata.get('status')!='DYNAMIC_TARGET_COMPLETE' or metadata.get('poi')!=poi or metadata.get('role')!=args.role or sha(path)!=metadata.get('graph_sha256'):raise ValidationError('Dynamic target graph receipt mismatch')
   inputs[str(receipt)]=sha(receipt);inputs[str(path)]=sha(path);dg[poi]=load_graph(path,poi,label,3).to(device)
 run=ROOT/'analysis/gnn'/args.role/f'seed_{args.seed}';run.mkdir(parents=True,exist_ok=True);intent=run/'intent.json'
 settings={'role':args.role,'seed':args.seed,'device':device,'epochs':args.epochs,'hidden':64,'heads':4,'dropout':.3,'lr':.001,'weight_decay':.0001,'cohort':cohort,'labels':labels,'inputs_sha256':inputs,'disabled_optional_backends':DISABLED,'torch':str(torch.__version__),'deterministic_algorithms':True,'purpose':'EXPLORATORY_NEW_REANALYSIS_NOT_HISTORICAL_REPLAY_OR_CALIBRATED_PROBABILITY'}
 if intent.exists() and json.loads(intent.read_text())!=settings:raise ValidationError('Training intent/input drift')
 atomic_json(intent,settings);rows=[];external=[]
 for held in cohort:
  checkpoint=run/f'fold_{held}.pt';receipt=run/f'fold_{held}.json'
  if receipt.exists():
   record=json.loads(receipt.read_text())
   if sha(checkpoint)!=record['checkpoint_sha256']:raise ValidationError('Fold checkpoint drift')
   state=torch.load(checkpoint,map_location=device,weights_only=True)
   if record.get('status')!='FOLD_COMPLETE' or state.get('settings')!=settings or state.get('held_out')!=held or state.get('training_targets')!=[p for p in cohort if p!=held] or state.get('seed')!=fold_seed(args.seed,held):raise ValidationError('Resumed fold provenance mismatch')
   replay=DualEncoderGNN(hidden=64,heads=4,dropout=.3).to(device);replay.load_state_dict(state['model_state_dict'])
   score=predict(replay,sg[held],dg.get(held))
   if state['prediction']!=record['prediction'] or record['prediction']['poi']!=held or record['prediction']['label']!=labels[held] or not np.isclose(score,record['prediction']['score'],rtol=0,atol=1e-7):raise ValidationError('Resumed fold replay mismatch')
   for prediction in record['external_predictions']:
    if prediction['held_out']!=held or prediction['poi'] not in prediction_targets or not np.isclose(predict(replay,sg[prediction['poi']],dg.get(prediction['poi'])),prediction['score'],rtol=0,atol=1e-7):raise ValidationError('Resumed external replay mismatch')
   rows.append(record['prediction']);external.extend(record['external_predictions']);continue
  seed=fold_seed(args.seed,held);random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
  if torch.cuda.is_available():torch.cuda.manual_seed_all(seed)
  model=DualEncoderGNN(hidden=64,heads=4,dropout=.3).to(device);optimizer=torch.optim.Adam(model.parameters(),lr=.001,weight_decay=.0001);scheduler=torch.optim.lr_scheduler.CosineAnnealingLR(optimizer,args.epochs)
  training=[p for p in cohort if p!=held];order=list(training);start=time.monotonic()
  for epoch in range(args.epochs):
   np.random.shuffle(order)
   for poi in order:
    model.train();optimizer.zero_grad();logit=model(sg[poi],dg.get(poi));loss=F.binary_cross_entropy_with_logits(logit,torch.tensor([float(labels[poi])],device=device));loss.backward();optimizer.step()
   scheduler.step()
   if (epoch+1)%50==0:print(json.dumps({'role':args.role,'held_out':held,'epoch':epoch+1,'epochs':args.epochs,'elapsed_s':round(time.monotonic()-start,1)}),flush=True)
  score=predict(model,sg[held],dg.get(held));prediction={'poi':held,'label':labels[held],'score':score,'pred':int(score>=.5),'has_md':held in dg,'fold_seed':seed}
  extra=[{'held_out':held,'poi':p,'score':predict(model,sg[p],dg.get(p)),'experimental_label':None,'purpose':'Prediction-only fold ensemble; not independent LOPO validation'} for p in prediction_targets]
  if held in training:raise ValidationError('Held-out leakage')
  state={'model_state_dict':{k:v.detach().cpu().clone() for k,v in model.state_dict().items()},'held_out':held,'training_targets':training,'seed':seed,'settings':settings,'prediction':prediction}
  if checkpoint.exists():
   archived=checkpoint.with_suffix('.incomplete_'+str(time.time_ns())+'.pt');os.replace(checkpoint,archived)
  temporary=checkpoint.with_suffix('.pending.pt');torch.save(state,temporary);os.replace(temporary,checkpoint)
  # Replay the saved weight state, not merely the in-memory score.
  replay=DualEncoderGNN(hidden=64,heads=4,dropout=.3).to(device);replay.load_state_dict(torch.load(checkpoint,map_location=device,weights_only=True)['model_state_dict']);replayed=predict(replay,sg[held],dg.get(held))
  if not np.isclose(score,replayed,rtol=0,atol=1e-7):raise ValidationError('Saved fold inference replay mismatch')
  record={'status':'FOLD_COMPLETE','prediction':prediction,'external_predictions':extra,'checkpoint_sha256':sha(checkpoint),'replayed_score':replayed,'elapsed_s':time.monotonic()-start}
  atomic_json(receipt,record);rows.append(prediction);external.extend(extra);print(json.dumps(record),flush=True)
 if any(sha(Path(p))!=v for p,v in inputs.items()):raise ValidationError('Source/graph drift during training')
 df=pd.DataFrame(rows);df.to_csv(run/'predictions.tsv',sep='\t',index=False);pd.DataFrame(external).to_csv(run/'prediction_only_fold_scores.tsv',sep='\t',index=False)
 group=subgroup_aurocs(df.label,df.score,{'md':df.has_md.to_numpy(dtype=bool),'nonmd':~df.has_md.to_numpy(dtype=bool)})
 summary={'status':'NEW_EXPLORATORY_LOPO_COMPLETE','role':args.role,'cases':len(df),'accuracy':float(accuracy_score(df.label,df.pred)),'AUROC':float(roc_auc_score(df.label,df.score)),'subgroups':group,'cohort':cohort,'historical_labels_unchanged':True,'CDK17_excluded_for_identity_mismatch':True,'new_MD_or_docking':False,'independent_validation':False,'settings':settings,'output_sha256':{p.name:sha(p) for p in run.iterdir() if p.is_file() and p.name!='completion.json'}}
 atomic_json(run/'completion.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--role',choices=['elob','vhl'],default='vhl');p.add_argument('--seed',type=int,default=42);p.add_argument('--epochs',type=int,default=400);p.add_argument('--device',choices=['cpu','cuda'],default='cpu');p.add_argument('--probe',action='store_true');a=p.parse_args()
 if a.epochs<1:raise ValidationError('Positive training budget required')
 run(a)
