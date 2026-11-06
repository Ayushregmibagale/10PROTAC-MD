#!/usr/bin/env python3
"""Synthetic architecture/runtime check; never fit or alter historical study data."""
import importlib.util,json,sys,argparse
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
if args.output.exists():raise SystemExit('Refusing to overwrite an existing control receipt')
r=Path(__file__).resolve().parents[1];w=r/'workflow';p=w/'scripts/04_gnn/01_sjf8240_e2_lopo.py'
s=importlib.util.spec_from_file_location('isolated_sjf_gnn',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
import torch,torch_geometric
from torch_geometric.data import Data
m.torch.set_num_threads(1);m.torch.manual_seed(392);m.torch.use_deterministic_algorithms(True)
e=torch.tensor([[0,1,1,2],[1,0,2,1]],dtype=torch.long)
a=Data(x=torch.randn(3,9),edge_index=e,edge_attr=torch.ones(4,3));b=Data(x=torch.randn(3,3),edge_index=e,edge_attr=torch.ones(4,3))
model=m.E2DualEncoderGNN(n_e2_scalars=4,hidden=4,heads=1,dropout=0.);opt=torch.optim.Adam(model.parameters(),lr=.001)
loss=m.train_step(model,opt,a,b,torch.zeros(1,4),1.,'cpu');prob=m.predict(model,a,b,torch.zeros(1,4),'cpu')
assert 0<=prob<=1 and all(torch.isfinite(x).all() for x in model.parameters())
result={'status':'SYNTHETIC_ARCHITECTURE_FORWARD_BACKWARD_PASSED','torch':torch.__version__,'torch_geometric':torch_geometric.__version__,'seed':392,'device':'cpu','optional_backends_disabled':m.DISABLED_OPTIONAL_GNN_BACKENDS,'synthetic_loss':loss,'synthetic_score':prob,'source_sha256':m.sha256(p),'contract_sha256':m.sha256(w/'scripts/lib/gnn_contract.py'),'scope':'Runtime/gradient check only; no study training or predictive performance claim'}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
