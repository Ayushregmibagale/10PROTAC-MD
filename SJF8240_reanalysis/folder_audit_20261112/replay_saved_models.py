"""Read-only CPU replay of every existing reviewed model; reports outside old runs."""
from pathlib import Path
import sys,json,argparse
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--version',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.root.resolve();version=a.version;run=root/version
sys.path.insert(0,str(run/'lib'));from common import *
torch.set_num_threads(1);rows=[]
for role in ['elob','vhl']:
 cfg,graphs=load_inputs(role,'cpu');verify_frozen()
 for held in cfg['cohort']:
  folder=run/'analysis'/role;weight=folder/f'fold_{held}.pt';rec=json.loads((folder/f'fold_{held}.json').read_text());require(sha(weight)==rec['checkpoint_sha256'],'Weight hash changed');state=torch.load(weight,map_location='cpu',weights_only=True);training=[x for x in cfg['cohort'] if x!=held];require(state['held_out']==held and state['training_targets']==training,'Fold exclusion mismatch')
  if cfg.get('E2_scalar_branch'):apply_e2_features(graphs,cfg,training,state['e2_normalization'])
  model=MDEncoder();model.load_state_dict(state['model_state_dict']);score=predict(model,graphs[held]);err=abs(score-rec['prediction']['score']);require(err<1e-5,'CPU replay mismatch');rows.append({'role':role,'held_out':held,'max_abs_error':err,'checkpoint_sha256':sha(weight)})
atomic_json(a.output,{'status':'PASS','version':version,'models':len(rows),'records':rows});print(version,len(rows),'passed')
