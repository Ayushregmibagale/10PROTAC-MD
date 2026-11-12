from pathlib import Path
import sys,json,argparse
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.root/'mechanistic_v3';sys.path.insert(0,str(root/'lib'));from common import *
from temporal_model import TemporalEncoder
from torch_geometric.data import Batch
rows=[];torch.set_num_threads(1)
for stage in ['baseline','temporal']:
 source=json.loads((root/'config'/('pair_source_hashes.json' if stage=='baseline' else 'temporal_source_hashes.json')).read_text())
 for path,h in source.items():require(sha(root/path)==h,'Source drift')
 for role in ['elob','vhl']:
  cfg,graphs=load_inputs(role,'cpu')
  if stage=='temporal':
   graphs={};manifest=json.loads((root/'config/temporal_input_manifest.json').read_text())
   for r in manifest['graphs']:
    if r['role']!=role:continue
    path=root/r['path'];require(sha(path)==r['sha256'],'Graph drift')
    with np.load(path) as z:graphs[r['target']]=Batch.from_data_list([Data(x=torch.tensor(z[f'x{w}'],dtype=torch.float32),edge_index=torch.tensor(z[f'edge_index{w}'],dtype=torch.long),edge_attr=torch.tensor(z[f'edge_attr{w}'],dtype=torch.float32)) for w in range(4)])
  for path in sorted((root/'analysis'/stage/role).glob('*.pt')):
   rec=json.loads(path.with_suffix('.json').read_text());require(sha(path)==rec['checkpoint_sha256'],'Weight drift');state=torch.load(path,map_location='cpu',weights_only=True);held=state['pair'];require(state['training']==[x for x in cfg['cohort'] if x not in held] and state['source_hashes']==source,'Fold/source drift');model=(MDEncoder() if stage=='baseline' else TemporalEncoder());model.load_state_dict(state['weights']);errors=[abs(predict(model,graphs[target])-score) for target,score in state['scores'].items()];require(max(errors)<1e-5,'Replay failure');rows.append({'stage':stage,'role':role,'fold':path.stem,'max_error':max(errors),'checkpoint_sha256':sha(path)})
require(len(rows)==48,'Missing pair models');atomic_json(a.output,{'status':'PASS','models':48,'records':rows});print('48 pair models passed')
