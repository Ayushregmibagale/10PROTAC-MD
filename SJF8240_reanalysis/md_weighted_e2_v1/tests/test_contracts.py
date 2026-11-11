from pathlib import Path
import sys,ast,json
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from common import load_graph,load_inputs,MDEncoder

def test_original_gat_backbone_preserved_with_two_scalar_head_inputs():
 from common import original_gat_backbone_matches
 assert original_gat_backbone_matches() and MDEncoder().head[0].in_features==130

def test_all_primary_targets_have_md_and_correct_user_labels():
 p,graphs=load_inputs('vhl');assert len(p['cohort'])==8 and all(x in graphs for x in p['cohort'])
 assert p['labels']['EPHA2']==1 and {k for k,v in p['labels'].items() if v==0}=={'ABL1','AXL'}
 assert all(graphs[k].y.item()==p['labels'][k] for k in p['cohort'])

def test_historical_epha2_label_cannot_enter_corrected_fit(tmp_path):
 path=tmp_path/'bad.npz';np.savez(path,poi=np.array(['EPHA2']),label=np.array([0]),x=np.zeros((2,3)),edge_index=np.array([[0,1],[1,0]]),edge_attr=np.ones((2,3)))
 with pytest.raises(ValueError,match='label mismatch'):load_graph(path,'EPHA2',1)

def test_fractional_raw_indices_cannot_be_coerced(tmp_path):
 path=tmp_path/'bad.npz';np.savez(path,poi=np.array(['EPHA2']),label=np.array([1]),x=np.zeros((2,3)),edge_index=np.array([[0.,1.2],[1.,0.]]),edge_attr=np.ones((2,3)))
 with pytest.raises(ValueError,match='indices must be integer'):load_graph(path,'EPHA2',1)
