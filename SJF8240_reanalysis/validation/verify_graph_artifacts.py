"""Audit complete retained graph packets and provenance without regenerating them."""
from pathlib import Path
import json,hashlib,datetime,collections
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as stream:
  for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def require(condition,message):
 if not condition:raise ValueError(message)
def validate_npz(path,width):
 with np.load(path,allow_pickle=False) as d:
  x=d['x'];ei=d['edge_index'];ea=d['edge_attr']
  require(x.ndim==2 and x.shape[1]==width and len(x)>0 and np.isfinite(x).all(),'Node schema invalid: '+str(path))
  require(ei.shape[0]==2 and ei.ndim==2 and ei.dtype.kind in 'iu' and (not ei.size or (ei.min()>=0 and ei.max()<len(x))),'Index schema invalid')
  require(ea.shape==(ei.shape[1],3) and np.isfinite(ea).all(),'Edge schema invalid')
  require(((ea[:,0]>=0)&(ea[:,0]<=1)).all(),'Contact frequencies outside [0,1]')
def validate_contact_packet(path):
 with np.load(path,allow_pickle=False) as d:
  required={'joint','sum_left','sum_right','freq_left','freq_right','left_resids','right_resids','times_ns'}
  require(set(d.files)==required,'Unexpected contact packet schema')
  data={k:d[k].copy() for k in d.files}
  require(all(np.isfinite(value).all() for value in data.values()),'Nonfinite contact packet')
  left=data['left_resids'];right=data['right_resids'];joint=data['joint'];fl=data['freq_left'];fr=data['freq_right']
  require(left.ndim==right.ndim==1 and len(np.unique(left))==len(left) and len(np.unique(right))==len(right),'Duplicate/invalid residue keys')
  require(left.dtype.kind in 'iu' and right.dtype.kind in 'iu','Residue identifiers must be integers')
  require(joint.shape==(len(left),len(right)) and fl.shape==data['sum_left'].shape==left.shape and fr.shape==data['sum_right'].shape==right.shape,'Contact packet shapes differ')
  require(data['times_ns'].shape==(201,) and np.allclose(data['times_ns'],np.linspace(0,20,201),rtol=0,atol=5.1e-6),'Packet time window mismatch')
  for value in [joint,fl,fr]:require(((value>=0)&(value<=1)).all() and np.allclose(value*201,np.rint(value*201),atol=1e-10,rtol=0),'Contact frequency not a 201-frame event count')
  require((joint<=np.minimum(fl[:,None],fr[None,:])+1e-12).all() and (joint>=np.maximum(0,fl[:,None]+fr[None,:]-1)-1e-12).all(),'Joint events violate marginal probability bounds')
  require((data['sum_left']>=0).all() and (data['sum_right']>=0).all(),'Negative contact distances')
  return data

def main():
 inventory=json.loads((ROOT/'config/source_inventory.json').read_text())['records'];static=[];files=0;accepted=0;excluded=0
 for cfg in inventory:
  poi=cfg['poi'];folder=ROOT/'analysis/static'/poi;receipt=json.loads((folder/'receipt.json').read_text());provenance_path=folder/'pose_provenance.json';provenance=json.loads(provenance_path.read_text())
  require(receipt['status']=='STATIC_PAIRED_COMPLETE' and receipt['poi']==poi,'Static target incomplete')
  require(sha(provenance_path)==receipt['pose_provenance_sha256'],'Provenance drift')
  listing=ROOT/cfg['pose_list'];require(sha(listing)==cfg['pose_list_sha256'],'Pose list drift')
  require([p['path'] for p in provenance]==listing.read_text().splitlines() and len(provenance)==cfg['poses'],'Configured source coverage mismatch')
  require(sha(cfg['first_pose'])==cfg['first_pose_sha256'],'Original reference structure changed')
  omitted=[p for p in provenance if not p.get('included',True)];included=len(provenance)-len(omitted)
  require(receipt.get('excluded_pose_count',0)==len(omitted),'Exclusion count mismatch')
  if omitted:
   ledger=folder/'excluded_poses.json';require(sha(ledger)==receipt['excluded_poses_sha256'],'Exclusion ledger drift');require(json.loads(ledger.read_text())==omitted,'Ledger/provenance mismatch')
   for record in omitted:
    require(record['exclusion_reason'] in {'EMPTY_FILE','TRUNCATED_ATOM_RECORD','INCOMPLETE_REFERENCE_NODE_SET'},'Unreviewed exclusion')
    require(sha(record['path'])==record['sha256'],'Excluded original file changed')
  for role in ['elob','vhl']:
   path=folder/role/'graph.npz';require(receipt['outputs'][role]['poses']==included and sha(path)==receipt['outputs'][role]['graph_sha256'],'Paired graph count/hash mismatch');validate_npz(path,9)
  seen=collections.Counter(p['sha256'] for p in provenance if p.get('included',True))
  static.append({'poi':poi,'configured_files':len(provenance),'accepted_complete_files':included,'excluded_files':len(omitted),'coordinate_overflow_records':sum(len(p.get('coordinate_overflow_records',[])) for p in provenance),'exact_duplicate_content_files_beyond_first':sum(n-1 for n in seen.values()),'receipt_sha256':sha(folder/'receipt.json')})
  files+=len(provenance);accepted+=included;excluded+=len(omitted)
 dynamic_path=ROOT/'analysis/dynamic/completion.json';dynamic=json.loads(dynamic_path.read_text());require(dynamic['status']=='DYNAMIC_COMPLETE' and dynamic['cases']==174 and len(dynamic['records'])==174,'Dynamic case coverage mismatch')
 identities=set();counts=collections.Counter()
 for row in dynamic['records']:
  key=(row['poi'],row['pose_id']);require(key not in identities,'Duplicate dynamic case');identities.add(key);counts[row['poi']]+=1
  folder=ROOT/'analysis/dynamic/poses'/row['poi']/row['pose_id'];require(json.loads((folder/'receipt.json').read_text())==row,'Dynamic receipt mismatch')
  require(row['frames']==201 and row['new_MD'] is False,'Dynamic frame/no-new-MD contract mismatch')
  for name,digest in row['outputs_sha256'].items():require(Path(name).name==name and sha(folder/name)==digest,'Dynamic pose graph drift');validate_contact_packet(folder/name)
 for poi,count in counts.items():
  for role in ['elob','vhl']:
   name=poi+'_full' if poi in ['SIK3','RPS6KB2'] else poi
   folder=ROOT/'analysis/dynamic/targets'/name/role;receipt=json.loads((folder/'receipt.json').read_text());require(receipt['poses']==count and receipt['role']==role and receipt['poi']==name and receipt['frames_per_pose']==201 and receipt['equal_pose_weight'] and receipt['joint_frequency_not_product_of_marginals'],'Dynamic aggregation mismatch');require(sha(folder/'graph.npz')==receipt['graph_sha256'],'Dynamic target graph drift');validate_npz(folder/'graph.npz',3)
   packets=[validate_contact_packet(ROOT/'analysis/dynamic/poses'/poi/row['pose_id']/(role+'.npz')) for row in dynamic['records'] if row['poi']==poi]
   first=packets[0]
   require(all(np.array_equal(p['left_resids'],first['left_resids']) and np.array_equal(p['right_resids'],first['right_resids']) for p in packets),'Per-pose residue correspondence mismatch')
   joint=np.mean([p['joint'] for p in packets],axis=0);left=np.mean([p['sum_left'] for p in packets],axis=0);right=np.mean([p['sum_right'] for p in packets],axis=0)
   i,j=np.nonzero(joint>=.01);n=len(left)
   expected_ei=np.array([np.concatenate([i,n+j]),np.concatenate([n+j,i])],dtype=np.int64)
   features=np.column_stack([joint[i,j],left[i]/10,right[j]/10]).astype(np.float32);expected_ea=np.concatenate([features,features])
   expected_x=np.zeros((len(left)+len(right),3),dtype=np.float32)
   for role_index,ids,start in [(0,first['left_resids'],0),(1,first['right_resids'],len(left))]:
    expected_x[start:start+len(ids),role_index]=1;expected_x[start:start+len(ids),2]=(ids-ids.min())/max(int(ids.max()-ids.min()),1)
   with np.load(folder/'graph.npz',allow_pickle=False) as graph:
    require(np.array_equal(graph['edge_index'],expected_ei) and np.array_equal(graph['edge_attr'],expected_ea) and np.array_equal(graph['x'],expected_x),'Target graph does not reproduce equal-pose packet aggregation')
 report={'status':'RETAINED_GRAPH_PACKET_AUDIT_PASS','checked_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'validator_sha256':sha(__file__),'static_targets':len(static),'configured_source_files':files,'accepted_complete_source_files':accepted,'excluded_unusable_source_files':excluded,'per_target':static,'dynamic_cases':174,'dynamic_targets':len(counts),'dynamic_case_counts':dict(counts),'dynamic_frames_per_case':201,'dynamic_case_frames_total':201*174,'dynamic_completion_sha256':sha(dynamic_path),'full_original_trajectory_rehash_performed':False,'all_derived_graphs_hash_checked':True,'dynamic_target_aggregation_independently_reproduced':True,'contact_packet_joint_marginal_bounds_checked':True,'new_simulation_or_training':False}
 destination=ROOT/'validation/graph_packet_verification.json';require(not destination.exists(),'Preserve existing validation receipt');destination.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='per_target'},indent=2))
if __name__=='__main__':main()
