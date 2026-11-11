from pathlib import Path
import sys,ast,json
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from common import load_graph,load_inputs,MDEncoder

def test_exact_original_md_only_architecture():
 def node(path):return next(x for x in ast.parse(path.read_text()).body if isinstance(x,ast.ClassDef) and x.name=='MDEncoder')
 assert ast.dump(node(ROOT/'lib/model.py'),include_attributes=False)==ast.dump(node(ROOT/'provenance/04h_md_only_ablation_original.py'),include_attributes=False)

def test_all_primary_targets_have_md_and_correct_user_labels():
 p,graphs=load_inputs('vhl');assert len(p['cohort'])==10 and all(x in graphs for x in p['cohort'])
 assert p['labels']['EPHA2']==1 and {k for k,v in p['labels'].items() if v==0}=={'ABL1','AXL','SIK3_full','RPS6KB2_full'}
 assert all(graphs[k].y.item()==p['labels'][k] for k in p['cohort'])

def test_historical_epha2_label_cannot_enter_corrected_fit(tmp_path):
 path=tmp_path/'bad.npz';np.savez(path,poi=np.array(['EPHA2']),label=np.array([0]),x=np.zeros((2,3)),edge_index=np.array([[0,1],[1,0]]),edge_attr=np.ones((2,3)))
 with pytest.raises(ValueError,match='label mismatch'):load_graph(path,'EPHA2',1)

def test_fractional_raw_indices_cannot_be_coerced(tmp_path):
 path=tmp_path/'bad.npz';np.savez(path,poi=np.array(['EPHA2']),label=np.array([1]),x=np.zeros((2,3)),edge_index=np.array([[0.,1.2],[1.,0.]]),edge_attr=np.ones((2,3)))
 with pytest.raises(ValueError,match='indices must be integer'):load_graph(path,'EPHA2',1)
