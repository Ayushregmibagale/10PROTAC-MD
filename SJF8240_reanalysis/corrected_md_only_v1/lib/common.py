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
 require(p['static_graph_branch'] is False and p['E2_scalar_branch'] is False and p['hyperparameter_tuning'] is False,'MD-only frozen branch required')
 require(sha(ROOT/'config/input_manifest.json')==p['input_manifest_sha256'] and sha(ROOT/'config/user_labels.json')==p['label_map_sha256'],'Label/input manifest drift')
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
 require(len(entries)==10 and {x['poi'] for x in entries}==set(p['cohort']+p['prediction_only']),'Complete primary and prediction-only graph sets required')
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
