import json,hashlib
import pytest
from test_static_integrity import builder,fixture

def test_build_and_resume_preserve_paired_validity_accounting(tmp_path,monkeypatch):
 data,_,_=fixture();valid=tmp_path/'valid.pdb';valid.write_bytes(data);empty=tmp_path/'empty.pdb';empty.write_bytes(b'')
 listing=tmp_path/'poses.txt';listing.write_text(str(valid)+'\n'+str(empty)+'\n')
 cfg={'poi':'fixture','kin_chain':'D','kin_min':0,'historical_label':0,'pose_list':'poses.txt','pose_list_sha256':builder.sha(listing),'first_pose_sha256':builder.sha(valid)}
 monkeypatch.setattr(builder,'ROOT',tmp_path)
 monkeypatch.setattr(builder,'validate_roles',lambda groups,poi:{'test_fixture_only':True})
 builder.build(cfg)
 folder=tmp_path/'analysis/static/fixture';receipt=json.loads((folder/'receipt.json').read_text())
 assert receipt['input_pose_count']==2 and receipt['included_pose_count']==1 and receipt['excluded_pose_count']==1
 assert receipt['outputs']['elob']['poses']==receipt['outputs']['vhl']['poses']==1
 excluded=json.loads((folder/'excluded_poses.json').read_text())
 assert len(excluded)==1 and excluded[0]['path']==str(empty) and excluded[0]['exclusion_reason']=='EMPTY_FILE'
 builder.build(cfg)  # Successful, hash-verified reuse of both completed arms.
 (folder/'vhl/graph.npz').write_bytes(b'corrupt')
 with pytest.raises(ValueError,match='graph drift'):builder.build(cfg)
