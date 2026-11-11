from pathlib import Path
import sys,os,json,random,argparse,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from common import *
from temporal_model import TemporalEncoder
from torch_geometric.data import Batch
MDEncoder=TemporalEncoder

def load_temporal(role,device):
 p=config();manifest=json.loads((ROOT/"config/temporal_input_manifest.json").read_text());graphs={}
 for row in manifest["graphs"]:
  if row["role"]!=role:continue
  path=ROOT/row["path"];require(sha(path)==row["sha256"],"Temporal graph drift")
  with np.load(path) as z:
   windows=[Data(x=torch.tensor(z[f"x{w}"],dtype=torch.float32),edge_index=torch.tensor(z[f"edge_index{w}"],dtype=torch.long),edge_attr=torch.tensor(z[f"edge_attr{w}"],dtype=torch.float32)) for w in range(4)]
  graphs[row["target"]]=Batch.from_data_list(windows).to(device)
 require(set(graphs)==set(p["cohort"]),"Incomplete temporal cohort");return p,graphs
import torch.nn.functional as F

def main(role,device,verify=False):
 torch.set_num_threads(1);torch.use_deterministic_algorithms(True)
 p,g=load_temporal(role,device);out=ROOT/'analysis'/'temporal'/role;out.mkdir(parents=True,exist_ok=True)
 frozen=json.loads((ROOT/'config/temporal_source_hashes.json').read_text())
 for path,digest in frozen.items():require(sha(ROOT/path)==digest,'Source drift '+path)
 rows=[];errors=[]
 for pos in [x for x in p['cohort'] if p['labels'][x]]:
  for neg in ['ABL1','AXL']:
   pair=pos+'__'+neg;training=[x for x in p['cohort'] if x not in [pos,neg]];seed=fold_seed(42,pair);weight=out/(pair+'.pt');receipt=out/(pair+'.json')
   if receipt.exists():
    rec=json.loads(receipt.read_text());require(sha(weight)==rec['checkpoint_sha256'],'Weight drift');state=torch.load(weight,map_location=device,weights_only=True);require(state['training']==training and state['pair']==[pos,neg],'Fold leakage');model=MDEncoder().to(device);model.load_state_dict(state['weights'])
   else:
    require(not verify,'Missing fold '+pair)
    random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
    if device=='cuda':torch.cuda.manual_seed_all(seed)
    model=MDEncoder().to(device);opt=torch.optim.Adam(model.parameters(),lr=.001,weight_decay=1e-4);sched=torch.optim.lr_scheduler.CosineAnnealingLR(opt,300);start=time.monotonic()
    for epoch in range(300):
     model.train()
     for target in training:
      opt.zero_grad();loss=F.binary_cross_entropy_with_logits(model(g[target]),torch.tensor([float(p['labels'][target])],device=device));require(bool(torch.isfinite(loss)),'Nonfinite loss');loss.backward();opt.step()
     sched.step()
    scores={x:predict(model,g[x]) for x in [pos,neg]}
    torch.save({'weights':{k:v.detach().cpu() for k,v in model.state_dict().items()},'training':training,'pair':[pos,neg],'seed':seed,'scores':scores,'source_hashes':frozen},weight)
    rec={'pair':pair,'positive':pos,'negative':neg,'positive_score':scores[pos],'negative_score':scores[neg],'prediction_only':{},'training':training,'training_counts':[1,5],'seed':seed,'checkpoint_sha256':sha(weight),'seconds':time.monotonic()-start};atomic_json(receipt,rec)
    model.load_state_dict(torch.load(weight,map_location=device,weights_only=True)['weights'])
   for target,key in [(pos,'positive_score'),(neg,'negative_score')]:
    err=abs(predict(model,g[target])-rec[key]);require(err<=(1e-5 if verify else 1e-7),'Replay failure');errors.append(err)
   rows.append(rec);print(json.dumps({'role':role,'pair':pair,'positive':rec['positive_score'],'negative':rec['negative_score'],'verification':verify}),flush=True)
 import pandas as pd
 ranking=sum(float(r['positive_score']>r['negative_score'])+.5*float(r['positive_score']==r['negative_score']) for r in rows)/12
 sens=sum(r['positive_score']>=.5 for r in rows)/12;spec=sum(r['negative_score']<.5 for r in rows)/12
 result={'role':role,'matched_pair_AUROC':ranking,'balanced_accuracy':(sens+spec)/2,'sensitivity':sens,'specificity':spec,'pairs':12,'independent_negative_targets':2,'max_replay_error':max(errors),'verified_cpu':verify,'source_hashes':frozen}
 if verify:atomic_json(out/'independent_cpu_verification.json',result)
 else:pd.DataFrame(rows).to_csv(out/'pairs.tsv',sep='\t',index=False);atomic_json(out/'completion.json',result)
 print(json.dumps(result),flush=True)
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--role',required=True,choices=['elob','vhl']);a.add_argument('--device',default='cuda');a.add_argument('--verify',action='store_true');args=a.parse_args();main(args.role,args.device,args.verify)
