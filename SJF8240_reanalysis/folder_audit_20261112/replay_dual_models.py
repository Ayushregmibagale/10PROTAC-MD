from pathlib import Path
import sys,json,argparse,importlib.util
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.root.resolve();path=root/'scripts/04_train_dual_lopo.py';spec=importlib.util.spec_from_file_location('historical_dual',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);torch=m.torch;torch.set_num_threads(1);rows=[]
for role in ['elob','vhl']:
 for path in sorted((root/'analysis/gnn'/role/'seed_42').glob('fold_*.pt')):
  rec=json.loads(path.with_suffix('.json').read_text());state=torch.load(path,map_location='cpu',weights_only=True);held=state['held_out'];settings=state['settings'];cohort=settings['cohort'];labels=settings['labels']
  if state['training_targets']!=[x for x in cohort if x!=held]:raise m.ValidationError('Fold leakage')
  sg=m.load_graph(root/'analysis/static'/held/role/'graph.npz',held,labels[held],9);dpath=root/'analysis/dynamic/targets'/held/role/'graph.npz';dg=m.load_graph(dpath,held,labels[held],3) if dpath.exists() else None;model=m.DualEncoderGNN(hidden=64,heads=4,dropout=.3);model.load_state_dict(state['model_state_dict']);error=abs(m.predict(model,sg,dg)-rec['prediction']['score'])
  if error>1e-5:raise m.ValidationError('CPU replay mismatch')
  if m.sha(path)!=rec['checkpoint_sha256']:raise m.ValidationError('Weight drift')
  rows.append({'role':role,'held_out':held,'max_error':error,'checkpoint_sha256':m.sha(path)})
if len(rows)!=32:raise m.ValidationError('Missing historical folds')
m.atomic_json(a.output,{'status':'PASS','models':32,'historical_labels_preserved':True,'records':rows});print('32 historical dual models passed')
