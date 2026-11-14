"""Freeze and verify matched graph,fold and label contracts before training."""
from pathlib import Path
import sys,json,ast
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'));from contracts import *
def main():
 p=protocol();m=json.loads((ROOT/'config/input_manifest.json').read_text());user=json.loads((ROOT/'config/user_labels.json').read_text())['labels'];require(p['labels']==user and len(user)==8 and set(k for k,v in user.items() if v==0)=={'ABL1','AXL'},'Label scope mismatch');require(len(m['entries'])==16,'Incomplete16graphs');by={}
 for r in m['entries']:
  path=ROOT/r['path'];require(sha(path)==r['sha256'],'Inputchanged');require(r['user_label']==user[r['target']],'Graphlabelmismatch');by[r['target'],r['variant']]=np.load(path)
 for target in p['cohort']:
  control=by[target,'zero_padded_control'];aug=by[target,'native_features']
  for w in range(4):
   a=aug[f'x{w}'];b=control[f'x{w}'];require(a.shape==b.shape and a.shape[1]==34 and np.isfinite(a).all() and np.isfinite(b).all(),'Graphfeature shape/values');require(np.array_equal(a[:,:29],b[:,:29]) and not np.any(b[:,29:]),'Unmatchedbase/control');require(np.all(a[:,29:][a[:,23]==0]==0),'NewfeaturesoutsidePOIlysinenodes')
   for key in ['edge_index','edge_attr']:require(np.array_equal(aug[f'{key}{w}'],control[f'{key}{w}']),'Unmatchededges')
   ei=aug[f'edge_index{w}'];ea=aug[f'edge_attr{w}'];require(ei.dtype.kind in 'iu' and ei.shape[0]==2 and ea.shape==(ei.shape[1],3) and np.isfinite(ea).all() and (not ei.size or (ei.min()>=0 and ei.max()<len(a))),'Invalidedges')
 for fold in p['folds']:
  held={fold['held_out_positive'],fold['held_out_negative']};train=set(fold['training']);require(not train&held and train|held==set(p['cohort']) and len(train)==6,'Foldleakage');require(sum(user[x]==0 for x in train)==1 and sum(user[x]==1 for x in train)==5,'Unexpectedtrainingclasscounts')
 old=Path(json.loads((ROOT/'config/source_provenance.json').read_text())['original_model_path']);require((ROOT/'lib/temporal_model.py').read_text()==old.read_text().replace('GATConv(29,16','GATConv(34,16'),'Architecturechangedbeyondinputwidth');f=frozen();atomic_json(ROOT/'validation/preflight.json',{'status':'PASS','graphs':16,'targets':8,'pairs_per_variant':12,'training_counts_each_fold':{'positive':5,'negative':1},'initial_seed_pairing_required':True,'source_manifest_sha256':f,'metadata_labels_used_for_training_only':True,'source_experiment_modified':False});print('PASS:16 matchedgraphs,12 whole-target heldoutpairs,unchangedmodel apartfrom34inputs.')
if __name__=='__main__':main()
