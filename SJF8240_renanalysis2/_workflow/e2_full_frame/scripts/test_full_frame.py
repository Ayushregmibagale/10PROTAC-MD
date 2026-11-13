from pathlib import Path
import sys,numpy as np,pytest
sys.path.insert(0,str(Path(__file__).parent))
from extract import paired_delta

def test_nearest_point_then_fraction_of_frames():
 distances=np.array([[20.,70.],[80.,10.],[30.,90.],[90.,20.]])
 result=paired_delta(distances,np.array([True,True,False,False]));assert np.array_equal(result,[0,0])

def test_known_direction_and_empty_state():
 d=np.array([[20.],[30.],[70.],[80.]]);assert paired_delta(d,[1,1,0,0])[0]==1;assert paired_delta(d,[0,0,1,1])[0]==-1;assert paired_delta(d,[1,1,1,1]) is None

def test_strict_cutoff():
 assert paired_delta(np.array([[50.],[49.999]]),[1,0])[0]==-1

def test_nonfinite_and_time_length_rejected():
 with pytest.raises(ValueError):paired_delta(np.array([[np.nan],[1.]]),[1,0])
 with pytest.raises(ValueError):paired_delta(np.array([[1.],[2.]]),[1])

def test_exact_protein_rank_and_ties():
 from analyze import exact_rank
 assert abs(exact_rank(np.arange(8),[0,0,1,1,1,1,1,1])['rank_two_sided_exact_p']-2/28)<1e-12
 assert exact_rank(np.zeros(8),[0,0,1,1,1,1,1,1])['rank_two_sided_exact_p']==1

def test_zero_is_not_positive_in_target_summary():
 import pandas as pd
 from analyze import summarize
 frame=pd.DataFrame([{'variant':'v','poi':'p','label':1,'resid':1,'delta_reach_50A':0.,'paired_poses':3},{'variant':'v','poi':'p','label':1,'resid':2,'delta_reach_50A':.1,'paired_poses':3}]);row=summarize(frame).iloc[0];assert row.positive_lysines==1 and row.zero_lysines==1 and row.fraction_positive==.5 and row.mean_delta==.05
