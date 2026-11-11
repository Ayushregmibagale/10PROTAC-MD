"""Fresh trajectory replay for one preselected first pose per target."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import extract_lysines as e
labels=json.loads((ROOT/'config/protocol.json').read_text())['labels'];cases=[]
for target in labels:
 pose=sorted((ROOT/'analysis/lysines/poses'/target).iterdir())[0].name;cases.append((target,pose))
e.ROOT=ROOT/'validation/trajectory_replays';rows=[]
for target,pose in cases:
 r=e.work((target,pose));old=np.load(ROOT/'analysis/lysines/poses'/target/pose/'packet.npz');new=np.load(e.ROOT/'analysis/lysines/poses'/target/pose/'packet.npz');assert set(old.files)==set(new.files)
 errors={}
 for key in old.files:
  if old[key].dtype.kind in 'fc':
   err=float(np.max(np.abs(old[key]-new[key]))) if old[key].size else 0.;assert err<1e-6;errors[key]=err
  else:assert np.array_equal(old[key],new[key])
 rows.append({'target':target,'pose':pose,'max_error':max(errors.values()),'arrays':len(old.files),'packet_identical':r['packet_sha256']==e.sha256(ROOT/'analysis/lysines/poses'/target/pose/'packet.npz')})
e.atomic_json(ROOT/'validation/trajectory_replay.json',{'status':'PASS','cases':8,'selection':'first lexicographic pose per confirmed target, no score-based selection','records':rows,'method':'fresh read and recomputation using frozen extraction implementation; not independent physical method'})
