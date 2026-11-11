"""Strict frozen inputs and model contracts for the corrected MD-only comparison."""
from pathlib import Path
import os,sys,json,hashlib
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'lib'))
from gnn_contract import configure_native_pyg,validate_feature_graph,fold_seed
DISABLED=configure_native_pyg()
from protac_qc import ValidationError,atomic_json
import numpy as np
import torch
from torch_geometric.data import Data
from model import MDEncoder,HIDDEN,HEADS,DROP
from sklearn.metrics import roc_auc_score,accuracy_score,balanced_accuracy_score,average_precision_score,confusion_matrix

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def require(ok,message):
 if not ok:raise ValidationError(message)
def verify_frozen():
 manifest=json.loads((ROOT/'config/execution_manifest.json').read_text())
 for path,digest in manifest.items():require(sha(ROOT/path)==digest,'Frozen source/configuration drift: '+path)
 return sha(ROOT/'config/execution_manifest.json')
def config():
 p=json.loads((ROOT/'config/protocol.json').read_text());labels=p['labels'];user=json.loads((ROOT/'config/user_labels.json').read_text())
 require(p['cohort']==sorted(labels) and len(labels)==8 and labels==user['labels'],'Eight-target user-confirmed cohort required')
 require({k for k,v in labels.items() if v==0}=={'ABL1','AXL'} and labels['EPHA2']==1 and all(type(v) is int and v in [0,1] for v in labels.values()),'User label/cohort mismatch')
 require((p['epochs'],p['hidden'],p['heads'],p['dropout'],p['seed'])==(300,HIDDEN,HEADS,DROP,42),'Frozen MD-only model settings differ')
 require(p['static_graph_branch'] is False and p['E2_scalar_branch'] is True and p['hyperparameter_tuning'] is False,'MD-only frozen branch required')
 require(sha(ROOT/'config/input_manifest.json')==p['input_manifest_sha256'] and sha(ROOT/'config/user_labels.json')==p['label_map_sha256'],'Label/input manifest drift')
 require(sha(ROOT/'inputs/e2_features.tsv')==p['E2_features_sha256'] and sha(ROOT/'provenance/e2_feature_provenance.json')==p['E2_provenance_sha256'],'E2 descriptor/provenance drift')
 return p

def load_graph(path,poi,label):
 with np.load(path,allow_pickle=False) as d:
  require(d['poi'].shape==(1,) and str(d['poi'][0])==poi,'Graph target mismatch')
  require(d['label'].shape==(1,) and np.isfinite(d['label'][0]) and d['label'][0]==label,'Graph label mismatch')
  require(d['edge_index'].dtype.kind in 'iu','Raw graph indices must be integer')
  graph=Data(x=torch.tensor(d['x'],dtype=torch.float32),edge_index=torch.tensor(d['edge_index'],dtype=torch.long),edge_attr=torch.tensor(d['edge_attr'],dtype=torch.float32),y=torch.tensor([label]),poi=poi,feature_schema='SJF8240-retained-paired-v1')
 validate_feature_graph(graph,3,schema='SJF8240-retained-paired-v1')
 require(bool(((graph.edge_attr[:,0]>=0)&(graph.edge_attr[:,0]<=1)).all()) and bool((graph.edge_attr[:,1:]>=0).all()),'Invalid MD edge features')
 return graph

def load_inputs(role,device='cpu'):
 p=config();manifest=json.loads((ROOT/'config/input_manifest.json').read_text());entries=[x for x in manifest['entries'] if x['role']==role]
 require(len(entries)==len(p['cohort'])+len(p['prediction_only']) and {x['poi'] for x in entries}==set(p['cohort']+p['prediction_only']),'Complete primary and prediction-only graph sets required')
 graphs={}
 for row in entries:
  path=ROOT/row['path'];require(path.resolve().is_relative_to((ROOT/'inputs').resolve()),'Unsafe input path')
  label=p['labels'].get(row['poi'],-1);require(row['label']==label and row['prediction_only']==(row['poi'] in p['prediction_only']),'Input manifest label/scope mismatch')
  require(sha(path)==row['sha256'],'Copied graph drift: '+row['poi'])
  graphs[row['poi']]=load_graph(path,row['poi'],label).to(device)
 return p,graphs

@torch.no_grad()
def predict(model,graph):
 model.eval();score=float(torch.sigmoid(model(graph)).item());require(np.isfinite(score) and 0<=score<=1,'Invalid inference score');return score

def metrics(rows):
 y=np.array([r['label'] for r in rows]);s=np.array([r['score'] for r in rows]);pred=(s>=.5).astype(int)
 require(y.shape==(8,) and np.isfinite(s).all() and set(y)=={0,1},'Complete finite eight-target predictions required')
 tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1]).ravel()
 return {'AUROC':float(roc_auc_score(y,s)),'accuracy':float(accuracy_score(y,pred)),'balanced_accuracy':float(balanced_accuracy_score(y,pred)),'average_precision':float(average_precision_score(y,s)),'sensitivity':float(tp/(tp+fn)),'specificity':float(tn/(tn+fp)),'confusion_matrix_labels_0_1':[[int(tn),int(fp)],[int(fn),int(tp)]],'positive_count':int((y==1).sum()),'negative_count':int((y==0).sum()),'majority_class_accuracy_reference':.75}

def balanced_class_weights(labels,training_targets):
 require(len(set(training_targets))==len(training_targets) and training_targets,'Unique nonempty training targets required')
 counts={y:sum(labels[poi]==y for poi in training_targets) for y in [0,1]}
 require(all(counts.values()),'Both classes required within training fold')
 return {y:len(training_targets)/(2*counts[y]) for y in [0,1]}

def balanced_bce(logit,label,weights):
 target=torch.tensor([float(label)],device=logit.device)
 return torch.nn.functional.binary_cross_entropy_with_logits(logit,target,weight=torch.tensor([weights[label]],device=logit.device))

def fit_e2_scaler(raw,training_targets):
 require(len(set(training_targets))==len(training_targets) and len(training_targets)>1,'Unique training targets required for E2 scaling')
 values=np.array([raw[poi] for poi in training_targets],dtype=float)
 require(values.shape==(len(training_targets),2) and np.isfinite(values).all(),'Finite two-feature training E2 matrix required')
 mean=values.mean(axis=0);scale=values.std(axis=0);scale=np.where(scale>1e-12,scale,1.)
 return {'training_targets':list(training_targets),'mean':mean.tolist(),'scale':scale.tolist()}

def apply_e2_features(graphs,p,training_targets,saved_scaler=None):
 import pandas as pd
 frame=pd.read_csv(ROOT/'inputs/e2_features.tsv',sep='\t');require(not frame.poi.duplicated().any() and set(frame.poi)==set(p['cohort']),'Complete unique E2 target coverage required')
 raw={row['poi']:[float(row[c]) for c in p['E2_columns']] for row in frame.to_dict('records')}
 require(all(np.isfinite(values).all() for values in raw.values()),'Missing E2 features cannot be silently imputed')
 scaler=fit_e2_scaler(raw,training_targets)
 if saved_scaler is not None:
  require(saved_scaler['training_targets']==scaler['training_targets'] and saved_scaler['feature_columns']==p['E2_columns'] and saved_scaler['feature_source_sha256']==p['E2_features_sha256'],'Saved E2 scaler lineage mismatch')
  require(np.allclose(saved_scaler['mean'],scaler['mean'],rtol=0,atol=1e-12) and np.allclose(saved_scaler['scale'],scaler['scale'],rtol=0,atol=1e-12),'E2 scaler not fitted on training targets only')
  scaler=saved_scaler
 mean=np.array(scaler['mean']);scale=np.array(scaler['scale'])
 for poi,graph in graphs.items():graph.e2=torch.tensor(((np.array(raw[poi])-mean)/scale).reshape(1,2),dtype=torch.float32,device=graph.x.device)
 return {**scaler,'feature_columns':p['E2_columns'],'feature_source_sha256':p['E2_features_sha256']}

def original_gat_backbone_matches():
 import ast
 def body(path):
  node=next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.ClassDef) and n.name=='MDEncoder')
  init=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
  forward=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name=='forward')
  statements=[n for n in init.body if not (isinstance(n,ast.Assign) and any(isinstance(t,ast.Attribute) and t.attr=='head' for t in n.targets))]+forward.body[:-1]
  return ast.dump(ast.Module(body=statements,type_ignores=[]),include_attributes=False)
 return body(ROOT/'lib/model.py')==body(ROOT/'provenance/model_before_E2.py')
