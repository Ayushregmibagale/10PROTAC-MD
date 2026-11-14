"""Prediction-only inference from all120 previously frozen GNNs; no fitting."""
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'lib'))
from runtime import configure_native_pyg
configure_native_pyg()
import numpy as np,pandas as pd,torch
from torch_geometric.data import Data,Batch
from temporal_model import TemporalEncoder
BASE=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2');ROBUST=BASE/'robustness_followup_v1';OLD=BASE/'matched_gnn_ablation_v1'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def graph(path):
 z=np.load(path);return Batch.from_data_list([Data(x=torch.tensor(z[f'x{w}'],dtype=torch.float32),edge_index=torch.tensor(z[f'edge_index{w}'],dtype=torch.long),edge_attr=torch.tensor(z[f'edge_attr{w}'],dtype=torch.float32)) for w in range(4)])
@torch.no_grad()
def predict(m,g):m.eval();return float(m(g).sigmoid().item())
def main():
 torch.set_num_threads(1);torch.use_deterministic_algorithms(True);check=json.loads((ROOT/'validation/feature_matching.json').read_text());assert check['status']=='PASS';graphs={}
 for row in check['unknown_graphs']:
  path=ROOT/row['path'];assert sha(path)==row['sha256'] and row['label'] is None;graphs[row['variant'],row['target']]=graph(path)
 records=json.loads((ROBUST/'completion.json').read_text());rows=[];replay=[]
 for r in records:
  seed,v,pair=r['base_seed'],r['variant'],r['pair'];path=ROBUST/'folds'/str(seed)/v/(pair+'.pt');assert sha(path)==r['weight_sha256'];state=torch.load(path,map_location='cpu',weights_only=True);assert not set(state['training'])&{'SIK3_full','RPS6KB2_full'};m=TemporalEncoder();m.load_state_dict(state['weights'])
  for t in ['SIK3_full','RPS6KB2_full']:
   g=graphs[v,t];score=predict(m,g);fresh=TemporalEncoder();fresh.load_state_dict(state['weights']);error=abs(predict(fresh,g)-score);assert error<=1e-7;rows.append({'seed':seed,'variant':v,'heldout_pair':pair,'target':t,'score':score,'fixed_threshold':.5,'fixed_degrader_call':int(score>=.5),'training_only_threshold':r['training_only_threshold'],'training_threshold_degrader_call':int(score>=r['training_only_threshold']),'experimental_label':None,'checkpoint_sha256':sha(path),'graph_sha256':next(x['sha256'] for x in check['unknown_graphs'] if x['target']==t and x['variant']==v),'replay_error':error})
  # Verify two original held-out predictions in the same runtime.
  err=max(abs(predict(m,graph(OLD/'inputs'/v/(t+'.npz')))-r['scores'][t]) for t in [r['positive'],r['negative']]);assert err<=1e-7;replay.append({'seed':seed,'variant':v,'pair':pair,'heldout_replay_max_error':err})
 df=pd.DataFrame(rows);df.to_csv(ROOT/'unknown_fold_predictions.tsv',sep='\t',index=False);summary=df.groupby(['target','variant']).agg(mean_score=('score','mean'),min_score=('score','min'),max_score=('score','max'),fixed_degrader_vote_fraction=('fixed_degrader_call','mean'),training_threshold_degrader_vote_fraction=('training_threshold_degrader_call','mean'),models=('score','size')).reset_index();summary.to_csv(ROOT/'unknown_prediction_summary.tsv',sep='\t',index=False);df.groupby(['target','variant','seed']).agg(mean_score=('score','mean'),fixed_degrader_vote_fraction=('fixed_degrader_call','mean'),training_threshold_degrader_vote_fraction=('training_threshold_degrader_call','mean')).reset_index().to_csv(ROOT/'unknown_predictions_by_seed.tsv',sep='\t',index=False)
 (ROOT/'validation/unknown_CPU_replay.json').write_text(json.dumps({'status':'PASS','frozen_models':120,'prediction_only_scores':len(rows),'unknown_labels_used':False,'training_run':False,'heldout_reference_replays':replay,'maximum_unknown_replay_error':float(df.replay_error.max())},indent=2));print(summary.to_string(index=False))
if __name__=='__main__':main()
