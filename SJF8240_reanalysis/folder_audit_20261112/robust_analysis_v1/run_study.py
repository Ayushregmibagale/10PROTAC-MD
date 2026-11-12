"""Configurable grouped GNN study; atomic/resumable folds, separate CPU replay."""
from pathlib import Path
import os,sys,json,argparse,random,fcntl,time,hashlib,tempfile
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
from contracts import *
from runtime import configure_native_pyg
DISABLED=configure_native_pyg()
import torch
from torch.nn import functional as F
from torch_geometric.data import Data,Batch
from models import MDEncoder,TemporalEncoder

def model_for(config):return (MDEncoder if config['model']=='original_md' else TemporalEncoder)()
def load_graphs(config,paths,role,device):
 out={};width,windows=(3,1) if config['model']=='original_md' else (29,4)
 for (target,r),p in paths.items():
  if role!=r:continue
  arrays=graph_arrays(p,width,windows);gs=[Data(x=torch.tensor(x,dtype=torch.float32),edge_index=torch.tensor(ei,dtype=torch.long),edge_attr=torch.tensor(ea,dtype=torch.float32)) for x,ei,ea in arrays];out[target]=(gs[0] if windows==1 else Batch.from_data_list(gs)).to(device)
 return out
@torch.no_grad()
def predict(model,g):
 model.eval();value=float(torch.sigmoid(model(g)).item());require(math.isfinite(value) and 0<=value<=1,'Nonfinite prediction');return value

def run(study,output,role,device,verify=False):
 c,labels,paths=validate_study(study);require(role in c['roles'],'Undeclared role');require(not verify or device=='cpu','Independent verifier must use CPU');output=Path(output).resolve();source=Path(__file__).parent.resolve()
 for p in [Path(study).resolve(),*paths.values(),*source.glob('*.py')]:require(not p.is_relative_to(output),'Output would overwrite inputs/source')
 output.mkdir(parents=True,exist_ok=True)
 lock=(output/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 frozen={'study':sha(study),'source':{p.name:sha(p) for p in source.glob('*.py')},'graphs':{target+'|'+r:sha(p) for (target,r),p in paths.items()},'role':role,'model':c['model'],'compound_id':c['compound_id']}
 intent=output/'intent.json'
 if intent.exists():require(json.loads(intent.read_text())==frozen,'Run intent drift')
 else:require(not verify,'Cannot verify absent run');atomic_json(intent,frozen)
 environment=output/'training_environment.json';current={'torch':str(torch.__version__),'python':sys.version,'device':device,'cuda':torch.version.cuda,'disabled_optional_backends':DISABLED}
 if verify:require(environment.is_file(),'Training runtime receipt missing')
 elif environment.exists():require(json.loads(environment.read_text())==current,'Resume runtime changed; use original environment or a new study output')
 else:atomic_json(environment,current)
 require(device!='cuda' or torch.cuda.is_available(),'CUDA requested but unavailable');torch.set_num_threads(1);torch.use_deterministic_algorithms(True);graphs=load_graphs(c,paths,role,device);records=[];maxerr=0.;settings=c['training'];start=time.monotonic()
 for fold in folds(labels,c['evaluation']):
  key=fold['id'];weight=output/(key+'.pt');receipt=output/(key+'.json');seed=int.from_bytes(hashlib.sha256((str(settings['seed'])+'|'+json.dumps(fold['held_out'],separators=(',',':'))).encode()).digest()[:4],'big')%2147483647
  if receipt.exists():
   rec=json.loads(receipt.read_text());require(sha(weight)==rec['checkpoint_sha256'],'Weight drift');state=torch.load(weight,map_location=device,weights_only=True);require(state['frozen']==frozen and state['fold']==fold and state['seed']==seed and rec['seed']==seed and rec['fold']==fold and state['scores']==rec['scores'],'Cached fold drift');model=model_for(c).to(device);model.load_state_dict(state['weights'])
  else:
   require(not verify,'Missing fold '+fold['name']);random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
   if device=='cuda':torch.cuda.manual_seed_all(seed)
   model=model_for(c).to(device);opt=torch.optim.Adam(model.parameters(),lr=settings['lr'],weight_decay=settings['weight_decay']);schedule=torch.optim.lr_scheduler.CosineAnnealingLR(opt,settings['epochs']);weights=balanced_weights(labels,fold['training']) if settings['loss']=='balanced' else {x:1. for x in fold['training']}
   for epoch in range(settings['epochs']):
    model.train()
    for target in fold['training']:
     opt.zero_grad();loss=F.binary_cross_entropy_with_logits(model(graphs[target]),torch.tensor([float(labels[target])],device=device))*weights[target];require(bool(torch.isfinite(loss)),'Nonfinite loss');loss.backward();opt.step()
    schedule.step()
   # Unknowns are prediction-only; label never enters loss or split.
   targets=fold['held_out']+[r['target_id'] for r in c['targets'] if r['label'] is None];scores={x:predict(model,graphs[x]) for x in targets};state={'weights':{k:v.detach().cpu().clone() for k,v in model.state_dict().items()},'frozen':frozen,'fold':fold,'seed':seed,'scores':scores,'weights_per_training_target':weights,'runtime':{'torch':str(torch.__version__),'device':device}}
   if weight.exists():weight.replace(weight.with_name(weight.stem+'.unreceipted_'+str(time.time_ns())+'.pt'))
   with tempfile.NamedTemporaryFile(dir=output,delete=False) as f:name=f.name
   torch.save(state,name);os.replace(name,weight);model.load_state_dict(torch.load(weight,map_location=device,weights_only=True)['weights']);rec={'fold':fold,'seed':seed,'scores':scores,'checkpoint_sha256':sha(weight),'label_status':{r['target_id']:r['label_status'] for r in c['targets']}};atomic_json(receipt,rec)
  require(set(rec['scores'])==set(fold['held_out']+[r['target_id'] for r in c['targets'] if r['label'] is None]),'Prediction target coverage changed')
  for target,score in rec['scores'].items():
   err=abs(predict(model,graphs[target])-score);require(err<=(1e-5 if verify else 1e-7),'Prediction replay failure');maxerr=max(maxerr,err)
  records.append(rec);print(json.dumps({'fold':fold['name'],'complete':len(records),'verify':verify}),flush=True)
 # Re-read all bound inputs after execution, including role-evidence files.
 validate_study(study);require(sha(study)==frozen['study'],'Study changed during execution');require(all(sha(source/p)==h for p,h in frozen['source'].items()),'Source changed during execution')
 if c['evaluation']=='leave_pair_out':metrics=rank_metrics([{'positive_score':r['scores'][r['fold']['held_out'][0]],'negative_score':r['scores'][r['fold']['held_out'][1]]} for r in records])
 else:
  from sklearn.metrics import roc_auc_score,confusion_matrix
  held=[r['fold']['held_out'][0] for r in records];s=[r['scores'][x] for r,x in zip(records,held)];y=[labels[x] for x in held];tn,fp,fn,tp=confusion_matrix(y,np.array(s)>=.5,labels=[0,1]).ravel();sens=float(tp/(tp+fn));spec=float(tn/(tn+fp));metrics={'pooled_LOPO_AUROC':float(roc_auc_score(y,s)),'sensitivity':sens,'specificity':spec,'balanced_accuracy':(sens+spec)/2,'warning':'Pooled LOPO ranking can be affected by fold class-prior differences; do not equate with matched-pair AUROC.'}
 report={'status':'CPU_VERIFIED' if verify else 'TRAINING_COMPLETE_CPU_REPLAY_REQUIRED','folds':len(records),'metrics':metrics,'max_replay_error':maxerr,'independent_negative_targets':sum(v==0 for v in labels.values()),'elapsed_s':time.monotonic()-start,'frozen':frozen};atomic_json(output/('cpu_verification.json' if verify else 'completion.json'),report);print(json.dumps(report),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--role',required=True);p.add_argument('--device',choices=['cpu','cuda'],default='cpu');p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.study,a.output,a.role,a.device,a.verify)
