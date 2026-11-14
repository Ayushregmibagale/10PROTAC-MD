"""Five-seed matched GNN robustness and training-only threshold calibration."""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='1'
import sys,json,random,hashlib,time,argparse,concurrent.futures,multiprocessing
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
from runtime import configure_native_pyg,fold_seed
configure_native_pyg()
import numpy as np,torch
from torch.nn import functional as F
from torch_geometric.data import Data,Batch
from temporal_model import TemporalEncoder
OLD=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2/matched_gnn_ablation_v1')
P=json.loads((OLD/'config/protocol.json').read_text());SEEDS=[42,7,123,2026,31415]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix('.tmp');q.write_text(json.dumps(x,indent=2,allow_nan=False));q.replace(p)
def freeze():
 m=json.loads((ROOT/'manifest.json').read_text())
 for p,h in m.items():
  if sha(p)!=h:raise RuntimeError('Source drift '+p)
 return sha(ROOT/'manifest.json')
def graphs(v):
 out={}
 for t in P['cohort']:
  z=np.load(OLD/'inputs'/v/(t+'.npz'));out[t]=Batch.from_data_list([Data(x=torch.tensor(z[f'x{w}'],dtype=torch.float32),edge_index=torch.tensor(z[f'edge_index{w}'],dtype=torch.long),edge_attr=torch.tensor(z[f'edge_attr{w}'],dtype=torch.float32)) for w in range(4)])
 return out
@torch.no_grad()
def score(m,g):
 m.eval();return {t:float(m(x).sigmoid().item()) for t,x in g.items()}
def threshold(scores,training):
 vals=np.unique([scores[t] for t in training]);candidates=np.r_[0.,(vals[:-1]+vals[1:])/2,1.]
 def key(c):
  sens=np.mean([scores[t]>=c for t in training if P['labels'][t]==1]);spec=np.mean([scores[t]<c for t in training if P['labels'][t]==0]);return ((sens+spec)/2,-abs(c-.5),-c)
 return float(max(candidates,key=key))
def digest(m):
 h=hashlib.sha256()
 for k,v in m.state_dict().items():h.update(k.encode());h.update(v.detach().numpy().tobytes())
 return h.hexdigest()
def work(task):
 seed,v,f=task;torch.set_num_threads(1);torch.use_deterministic_algorithms(True);binding=freeze();out=ROOT/'folds'/str(seed)/v/f['pair'];out.parent.mkdir(parents=True,exist_ok=True);receipt=out.with_suffix('.json');weight=out.with_suffix('.pt')
 if receipt.exists():
  r=json.loads(receipt.read_text());assert r['manifest_sha256']==binding and r['weight_sha256']==sha(weight);return r
 s=fold_seed(seed,f['pair']);random.seed(s);np.random.seed(s);torch.manual_seed(s);g=graphs(v);m=TemporalEncoder();initial=digest(m);opt=torch.optim.Adam(m.parameters(),lr=.001,weight_decay=1e-4);sched=torch.optim.lr_scheduler.CosineAnnealingLR(opt,300);start=time.monotonic()
 for epoch in range(300):
  m.train()
  for t in f['training']:
   opt.zero_grad();loss=F.binary_cross_entropy_with_logits(m(g[t]),torch.tensor([float(P['labels'][t])]));assert torch.isfinite(loss);loss.backward();opt.step()
  sched.step()
 scores=score(m,g);cut=threshold(scores,f['training']);torch.save({'weights':m.state_dict(),'seed':s,'training':f['training'],'scores':scores},weight)
 r={'base_seed':seed,'variant':v,'pair':f['pair'],'positive':f['held_out_positive'],'negative':f['held_out_negative'],'training':f['training'],'scores':scores,'training_only_threshold':cut,'threshold_fit':'resubstitution_on_six_training_proteins_not_inner_CV','initial_digest':initial,'weight_sha256':sha(weight),'manifest_sha256':binding,'seconds':time.monotonic()-start,'epochs':300,'steps':1800};save(receipt,r);return r

def main():
 a=argparse.ArgumentParser();a.add_argument('--workers',type=int,default=4);args=a.parse_args();freeze();tasks=[(s,v,f) for s in SEEDS for f in P['folds'] for v in P['variants']];records=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn')) as ex:
  for future in concurrent.futures.as_completed([ex.submit(work,t) for t in tasks]):
   r=future.result();records.append(r);save(ROOT/'STATUS.json',{'completed':len(records),'expected':len(tasks),'phase':'training'});print(len(records),r['base_seed'],r['variant'],r['pair'],round(r['seconds'],1),flush=True)
 save(ROOT/'completion.json',records)
if __name__=='__main__':main()
