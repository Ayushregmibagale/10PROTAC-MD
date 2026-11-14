"""Deterministic24-fold CPU ablation;fixed300epochs,no outcome-based selection."""
from pathlib import Path
import os,sys,json,random,time,argparse,hashlib,concurrent.futures,multiprocessing
os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from contracts import *
from runtime import configure_native_pyg,fold_seed
DISABLED=configure_native_pyg()
import numpy as np,torch
from torch.nn import functional as F
from torch_geometric.data import Data,Batch
from temporal_model import TemporalEncoder

def state_digest(model):
 h=hashlib.sha256()
 for name,value in model.state_dict().items():h.update(name.encode());h.update(str(tuple(value.shape)).encode());h.update(value.detach().cpu().contiguous().numpy().tobytes())
 return h.hexdigest()
def graphs(variant,p):
 data={};manifest=json.loads((ROOT/'config/input_manifest.json').read_text())
 for row in manifest['entries']:
  if row['variant']!=variant:continue
  path=ROOT/row['path'];require(sha(path)==row['sha256'],'Graphdrift');z=np.load(path);windows=[Data(x=torch.tensor(z[f'x{w}'],dtype=torch.float32),edge_index=torch.tensor(z[f'edge_index{w}'],dtype=torch.long),edge_attr=torch.tensor(z[f'edge_attr{w}'],dtype=torch.float32)) for w in range(4)];data[row['target']]=Batch.from_data_list(windows)
 require(set(data)==set(p['cohort']),'Incompletecohort');return data
@torch.no_grad()
def predict(model,g):
 model.eval();result=float(torch.sigmoid(model(g)).item());require(np.isfinite(result) and 0<=result<=1,'Invalidscore');return result

def train(task):
 variant,fold=task;torch.set_num_threads(1);torch.use_deterministic_algorithms(True);p=protocol();f=frozen();out=ROOT/'results'/variant;out.mkdir(exist_ok=True,parents=True);pair=fold['pair'];receipt=out/(pair+'.json');weight=out/(pair+'.pt');seed=fold_seed(p['seed'],pair)
 if receipt.exists():
  r=json.loads(receipt.read_text());require(sha(weight)==r['checkpoint_sha256'] and r['source_manifest_sha256']==f,'Cachedfold drift');return r
 require(not weight.exists(),'Orphancheckpoint;inspect before overwriting');random.seed(seed);np.random.seed(seed);torch.manual_seed(seed);g=graphs(variant,p);model=TemporalEncoder();initial=state_digest(model);optimizer=torch.optim.Adam(model.parameters(),lr=p['learning_rate'],weight_decay=p['weight_decay']);scheduler=torch.optim.lr_scheduler.CosineAnnealingLR(optimizer,p['epochs']);start=time.monotonic();losses=[];training=fold['training'];held=[fold['held_out_positive'],fold['held_out_negative']];require(not set(training)&set(held),'Heldouttargetintraining')
 for epoch in range(p['epochs']):
  model.train();last=[]
  for target in training:
   optimizer.zero_grad();loss=F.binary_cross_entropy_with_logits(model(g[target]),torch.tensor([float(p['labels'][target])]));require(bool(torch.isfinite(loss)),'Nonfiniteloss');loss.backward();optimizer.step();last.append(float(loss.detach()))
  scheduler.step()
  if epoch==0 or (epoch+1)%25==0:
   losses.append({'epoch':epoch+1,'mean_training_loss':float(np.mean(last))});atomic_json(out/(pair+'.progress.json'),{'variant':variant,'pair':pair,'epoch':epoch+1,'epochs':p['epochs'],'elapsed_seconds':time.monotonic()-start,'mean_training_loss':float(np.mean(last))})
 scores={x:predict(model,g[x]) for x in held};state={'weights':{k:v.detach().cpu() for k,v in model.state_dict().items()},'variant':variant,'pair':held,'training':training,'seed':seed,'scores':scores,'initial_state_digest':initial,'source_manifest_sha256':f};temp=weight.with_suffix('.pending.pt');torch.save(state,temp);os.replace(temp,weight);r={'variant':variant,'pair':pair,'positive':held[0],'negative':held[1],'positive_score':scores[held[0]],'negative_score':scores[held[1]],'training':training,'training_label_counts':{'positive':5,'negative':1},'epochs':p['epochs'],'optimizer_steps':len(training)*p['epochs'],'seed':seed,'initial_state_digest':initial,'final_state_digest':state_digest(model),'checkpoint_sha256':sha(weight),'source_manifest_sha256':f,'seconds':time.monotonic()-start,'training_loss_trace':losses,'CPU_training':True};atomic_json(receipt,r)
 # Immediate new-instance replay in this worker;independentprocess replay is a separate step.
 replay=TemporalEncoder();replay.load_state_dict(torch.load(weight,map_location='cpu',weights_only=True)['weights']);errors=[abs(predict(replay,g[x])-scores[x]) for x in held];require(max(errors)<=1e-7,'Immediatecheckpoint replay failed');atomic_json(out/(pair+'.replay.json'),{'status':'PASS','maximum_absolute_score_error':max(errors),'checkpoint_sha256':sha(weight)});return r

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--workers',type=int,default=2);parser.add_argument('--limit',type=int);a=parser.parse_args();require(a.workers>0,'Positiveworkersrequired');require(json.loads((ROOT/'validation/preflight.json').read_text())['status']=='PASS','Preflightrequired');p=protocol();tasks=[(variant,fold) for fold in p['folds'] for variant in p['variants']];tasks=tasks[:a.limit] if a.limit else tasks;rows=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers,mp_context=multiprocessing.get_context('spawn')) as executor:
  for future in concurrent.futures.as_completed([executor.submit(train,task) for task in tasks]):r=future.result();rows.append(r);print(json.dumps({'completed':len(rows),'expected':len(tasks),'variant':r['variant'],'pair':r['pair'],'positive_score':r['positive_score'],'negative_score':r['negative_score'],'seconds':r['seconds']}),flush=True);atomic_json(ROOT/'STATUS.json',{'phase':'TRAINING','completed':len(rows),'expected':len(tasks),'all24_required':True})
 atomic_json(ROOT/('diagnostic_completion.json' if a.limit else 'training_completion.json'),{'folds':len(rows),'records':rows,'source_manifest_sha256':frozen(),'device':'cpu','epochs':300});print('Training complete:',len(rows),'folds',flush=True)
if __name__=='__main__':main()
