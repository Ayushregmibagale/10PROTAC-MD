"""Fresh8-case standard-MDAnalysis-unwrap check against full-frame tree packets."""
from pathlib import Path
import sys,json,concurrent.futures
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));from extract import work,require,atomic_json,sha256,CAN

def replay(case):return work(case,independent=True,output=ROOT/'validation/standard_unwrap/poses')
def main():
 labels=json.loads((ROOT/'config/protocol.json').read_text())['labels'];records=json.loads((CAN/'analysis/dynamic/completion.json').read_text())['records'];cases=[(target,sorted(r['pose_id'] for r in records if r['poi']==target)[0]) for target in sorted(labels)]
 with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:rows=list(pool.map(replay,cases))
 result=[]
 for r in rows:
  a=ROOT/'analysis/poses'/r['poi']/r['pose']/'packet.npz';b=ROOT/'validation/standard_unwrap/poses'/r['poi']/r['pose']/'packet.npz';require(a.exists(),'Main packet not yet ready');require(r['independent_MDA_unwrap'],'Wrong replay method')
  with np.load(a) as x,np.load(b) as y:
   errors={};flips={}
   for key in ['times_ns','lys_resids','legacy_mask','true_vhl_mask']:require(np.array_equal(x[key],y[key]),'Replay correspondence/mask mismatch')
   for key in ['legacy_min_e2_A','corrected_min_e2_A']:
    errors[key]=float(abs(x[key]-y[key]).max());flips[key]=int(((x[key]<50)!=(y[key]<50)).sum());require(errors[key]<.002 and flips[key]==0,'Standard unwrap changed geometry/threshold results')
  result.append({'poi':r['poi'],'pose':r['pose'],'frames':201,'distance_errors_A':errors,'threshold_classification_changes':flips,'main_packet_sha256':sha256(a),'replay_packet_sha256':sha256(b)})
 atomic_json(ROOT/'validation/independent_unwrap_replay.json',{'status':'PASS','cases':8,'frames':1608,'selection':'First lexicographic retained pose per target, selected without score access','scope':'Fresh read with standard MDAnalysis unwrap; shared reviewed Kabsch/imaging helpers. Numerical cross-check, not biological validation.','records':result});print('8 full-frame standard-unwrap replays passed')
if __name__=='__main__':main()
