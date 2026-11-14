"""Independent CPU checkpoint replay,initialization checks and whole-target split audit."""
from pathlib import Path
import sys,json,hashlib,random
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from contracts import *
from runtime import configure_native_pyg,fold_seed
configure_native_pyg()
import torch
from torch_geometric.data import Data,Batch
from temporal_model import TemporalEncoder

def digest(model):
 h=hashlib.sha256()
 for name,value in model.state_dict().items():h.update(name.encode());h.update(str(tuple(value.shape)).encode());h.update(value.detach().cpu().contiguous().numpy().tobytes())
 return h.hexdigest()
def main():
 torch.set_num_threads(1);torch.use_deterministic_algorithms(True);p=protocol();f=frozen();complete=json.loads((ROOT/'training_completion.json').read_text());require(complete['folds']==24,'All24folds required');by={(r['pair'],r['variant']):r for r in complete['records']};results=[];initial_by_pair={};scored=[]
 for fold in p['folds']:
  pair=fold['pair'];seed=fold_seed(p['seed'],pair)
  for variant in p['variants']:
   r=by[pair,variant];weight=ROOT/'results'/variant/(pair+'.pt');require(sha(weight)==r['checkpoint_sha256'],'Checkpointdrift');state=torch.load(weight,map_location='cpu',weights_only=True);held=[fold['held_out_positive'],fold['held_out_negative']];require(state['training']==fold['training'] and state['pair']==held and state['variant']==variant,'Storedsplitmismatch');require(not set(held)&set(state['training']) and len(state['training'])==6,'Heldoutleakage');require(state['source_manifest_sha256']==r['source_manifest_sha256']==f,'Sourcemismatch');require(state['seed']==r['seed']==seed,'Seedmismatch');require(r['epochs']==300 and r['optimizer_steps']==1800,'Trainingbudgetchanged');random.seed(seed);np.random.seed(seed);torch.manual_seed(seed);model=TemporalEncoder();initial=digest(model);require(initial==r['initial_state_digest']==state['initial_state_digest'],'Initialstate not reproducible');initial_by_pair.setdefault(pair,initial);require(initial_by_pair[pair]==initial,'Unmatchedinitialization');model.load_state_dict(state['weights']);require(digest(model)==r['final_state_digest'],'Finalstatesemanticdrift');model.eval();errors=[]
   for target,key in zip(held,['positive_score','negative_score']):
    z=np.load(ROOT/'inputs'/variant/(target+'.npz'));windows=[]
    for w in range(4):windows.append(Data(x=torch.from_numpy(z[f'x{w}'].astype(np.float32)),edge_index=torch.from_numpy(z[f'edge_index{w}'].astype(np.int64)),edge_attr=torch.from_numpy(z[f'edge_attr{w}'].astype(np.float32))))
    graph=Batch.from_data_list(windows)
    with torch.no_grad():score=float(torch.sigmoid(model(graph)).item())
    err=abs(score-r[key]);require(err<=1e-7,'CPUreplay mismatch');require(abs(score-state['scores'][target])<=1e-7,'Embeddedscore mismatch');errors.append(err);scored.append({'pair':pair,'variant':variant,'heldout_target':target,'user_label':p['labels'][target],'replayed_score':score,'score_error':err})
   results.append({'pair':pair,'variant':variant,'status':'CPU_VERIFIED','max_score_error':max(errors),'checkpoint_sha256':sha(weight),'same_seed_and_initial_state':True,'train_heldout_disjoint':True});print(variant,pair,'CPU_VERIFIED',flush=True)
 import pandas as pd
 pd.DataFrame(scored).to_csv(ROOT/'validation/independent_cpu_predictions.tsv',sep='\t',index=False);atomic_json(ROOT/'validation/independent_cpu_verification.json',{'status':'PASS','models':24,'heldout_scores':48,'records':results,'max_score_error':max(r['max_score_error'] for r in results),'same_initialization_all12_pairs':True,'source_manifest_sha256':f,'verification_script_sha256':sha(__file__),'labelled_targets':8,'independent_negative_targets':2,'biological_validation':False});print('PASS:24weights,48scores,allpairedinitializations andsplits independentlyreplayed.')
if __name__=='__main__':main()
