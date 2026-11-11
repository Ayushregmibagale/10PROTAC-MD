from pathlib import Path
import importlib.util,sys
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/lib'))
spec=importlib.util.spec_from_file_location('training',ROOT/'scripts/04_train_dual_lopo.py')
training=importlib.util.module_from_spec(spec);spec.loader.exec_module(training)

def write_graph(path,label):
 np.savez(path,poi=np.array(['MET']),label=np.array([label]),x=np.zeros((2,9)),edge_index=np.array([[0,1],[1,0]],dtype=np.int64),edge_attr=np.ones((2,3)))

def test_fractional_label_cannot_be_coerced_into_training_class(tmp_path):
 path=tmp_path/'graph.npz';write_graph(path,1.2)
 with pytest.raises(ValueError,match='label mismatch'):training.load_graph(path,'MET',1,9)

def test_wrong_target_cannot_enter_fold(tmp_path):
 path=tmp_path/'graph.npz';write_graph(path,1)
 with pytest.raises(ValueError,match='target mismatch'):training.load_graph(path,'AXL',1,9)

def test_valid_raw_graph_loads_with_exact_label(tmp_path):
 path=tmp_path/'graph.npz';write_graph(path,1)
 graph=training.load_graph(path,'MET',1,9)
 assert graph.y.item()==1 and graph.x.shape==(2,9)
