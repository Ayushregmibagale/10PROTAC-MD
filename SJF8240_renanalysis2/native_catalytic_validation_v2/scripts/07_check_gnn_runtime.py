"""Matched-input and finite-forward checks only; no optimization/training."""
from pathlib import Path
import sys,json,hashlib,importlib.util
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OLD=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis/mechanistic_v3');sys.path.insert(0,str(OLD/'lib'));from gnn_contract import configure_native_pyg
DISABLED=configure_native_pyg()
import torch
from torch_geometric.data import Data,Batch
from torch_geometric.nn import GATConv
spec=importlib.util.spec_from_file_location('original_temporal',OLD/'lib/temporal_model.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 torch.set_num_threads(1);torch.manual_seed(42);manifest=json.loads((ROOT/'validation/gnn_input_manifest.json').read_text());model=m.TemporalEncoder();model.c1=GATConv(34,16,heads=4,edge_dim=64,dropout=.3);model.eval();rows=[]
 for target in sorted({r['target'] for r in manifest['graphs']}):
  source=np.load(OLD/'inputs/temporal/vhl'/(target+'.npz'));variants={v:np.load(ROOT/'gnn_inputs'/v/(target+'.npz')) for v in ['native_features','zero_padded_control']}
  for variant,z in variants.items():
   windows=[]
   for w in range(4):
    x=z[f'x{w}'];ei=z[f'edge_index{w}'];ea=z[f'edge_attr{w}'];assert x.shape[1]==34 and np.isfinite(x).all() and np.isfinite(ea).all();assert np.array_equal(x[:,:29],source[f'x{w}']);assert np.array_equal(ei,source[f'edge_index{w}']) and np.array_equal(ea,source[f'edge_attr{w}']);assert ei.min()>=0 and ei.max()<len(x)
    if variant=='zero_padded_control':assert np.count_nonzero(x[:,29:])==0
    windows.append(Data(x=torch.tensor(x),edge_index=torch.tensor(ei),edge_attr=torch.tensor(ea)))
   with torch.no_grad():out=model(Batch.from_data_list(windows))
   assert out.shape==(1,) and torch.isfinite(out).all();row=next(r for r in manifest['graphs'] if r['target']==target and r['variant']==variant);assert sha(ROOT/row['path'])==row['sha256'];rows.append({'target':target,'variant':variant,'windows':4,'nodes_per_window':len(x),'input_width':34,'base29_and_edges_unchanged':True,'finite_untrained_forward':True})
 result={'status':'PASS','graphs':rows,'training_run':False,'optimizer_steps':0,'trained_predictions_saved':False,'parameters':sum(p.numel() for p in model.parameters()),'architecture':'Existing TemporalEncoder,onlyGATConv input width29to34; identical34-input architecture forbothvariants','model_source_sha256':sha(OLD/'lib/temporal_model.py'),'runtime_guard_sha256':sha(OLD/'lib/gnn_contract.py'),'disabled_optional_native_modules':DISABLED,'runtime_note':'An initial unguardedPyG import segfaulted due incompatible optional native extension. Existing process-local configure_native_pyg applied before import; no installed packages modified.','biological_validation':False};(ROOT/'validation/gnn_runtime_check.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS:16 matched34-feature graphs,original29features/edges preserved,finite untrained CPUforward;zero trainingsteps.')
if __name__=='__main__':main()
