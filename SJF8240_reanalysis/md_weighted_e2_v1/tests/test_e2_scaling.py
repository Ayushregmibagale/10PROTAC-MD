from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from common import fit_e2_scaler,load_inputs,apply_e2_features
import numpy as np
def test_held_out_feature_cannot_affect_fit_statistics():
 raw={'a':[1.,2.],'b':[3.,4.],'held':[1e9,-1e9]}
 a=fit_e2_scaler(raw,['a','b']);raw['held']=[-5.,999.]
 assert a==fit_e2_scaler(raw,['a','b']) and a['mean']==[2.,3.]
def test_zero_variance_uses_explicit_unit_scale():
 assert fit_e2_scaler({'a':[1.,2.],'b':[1.,4.]},['a','b'])['scale'][0]==1.
def test_actual_training_scalars_normalize_with_no_held_out_fit():
 p,graphs=load_inputs('vhl');train=p['cohort'][1:];scaler=apply_e2_features(graphs,p,train)
 values=np.concatenate([graphs[x].e2.numpy() for x in train],axis=0)
 assert np.allclose(values.mean(axis=0),0,atol=1e-6) and np.allclose(values.std(axis=0),1,atol=1e-6)
 assert scaler['training_targets']==train and all(g.e2.shape==(1,2) for g in graphs.values())
