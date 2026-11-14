"""Prepare matched34-feature temporal inputs; no fitting or model replacement."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OLD=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis/mechanistic_v3')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 completion=json.loads((ROOT/'completion.json').read_text());labels=json.loads((ROOT/'config/labels.json').read_text())['labels'];manifest=[];features=['mean_NZ_to_SY8C3_A/100','NZ_reachable20_fraction','NZ_exposed_fraction','exposed_reachable_scaffold_compatible_fraction','joint_trueVHLbridge_productive_fraction']
 for target in labels:
  source=OLD/'inputs/temporal/vhl'/(target+'.npz');base=np.load(source);cases=[r for r in completion['records'] if r['poi']==target];oldmeta=np.load(OLD/'analysis/lysines/poses'/target/cases[0]['pose']/'packet.npz');lys_indices=oldmeta['lys_node_indices'];lys_resids=oldmeta['lys_resids'];packets=[]
  for r in cases:
   path=ROOT/'results/poses'/target/r['pose']/'packet.npz';assert sha(path)==r['packet_sha256'];z=np.load(path);p={k:z[k] for k in z.files};z.close();assert np.array_equal(p['lys_resids'],lys_resids),'Residue identity drift';packets.append(p)
  variants={'native_features':{},'zero_padded_control':{}}
  for window in range(4):
   data=[]
   for p in packets:
    tt=p['times_ns'];mask=(tt>=5*window)&((tt<5*(window+1)) if window<3 else tt<=20);assert mask.sum()==(50 if window<3 else 51);d=p['donor_distance_A'][0,mask,0,:];ex=p['nz_exposed_fraction'][mask];clear=p['target_scaffold_max_overlap_A'][0,mask]<=.8;bridge=p['bridge'][mask];joint=(d<=20)&(ex>=.1)&clear[:,None];data.append(np.stack([d.mean(0)/100,(d<=20).mean(0),ex.mean(0),joint.mean(0),(joint&bridge[:,None]).mean(0)],axis=1))
   original=base[f'x{window}'];assert original.shape[1]==29;nl=int(original[:,0].sum());assert np.all(original[nl:,1]==1);assert np.all(original[nl+lys_indices,23]==1);added=np.zeros((len(original),5),np.float32);added[nl+lys_indices]=np.mean(data,axis=0);assert np.isfinite(added).all()
   for variant,addition in [('native_features',added),('zero_padded_control',np.zeros_like(added))]:
    variants[variant][f'x{window}']=np.c_[original,addition].astype(np.float32)
    for name in ['edge_index','edge_attr']:variants[variant][f'{name}{window}']=base[f'{name}{window}'].copy()
  for variant,arrays in variants.items():
   dest=ROOT/'gnn_inputs'/variant;dest.mkdir(exist_ok=True,parents=True);path=dest/(target+'.npz');np.savez_compressed(path,**arrays);manifest.append({'target':target,'user_label':labels[target],'variant':variant,'path':str(path.relative_to(ROOT)),'sha256':sha(path),'source_graph_sha256':sha(source),'poses':len(cases),'node_features':34})
 result={'graphs':manifest,'added_features':features,'base29_features_and_edges_unchanged':True,'matched_architecture_input_width':34,'windows_ns':[[0,5],[5,10],[10,15],[15,20]],'full201frame_data':True,'perlysine_equal_pose_weight':True,'training_run':False,'biological_release':False,'label_based_feature_selection':False,'note':'Provisional engineered lysine-node features; no full scaffold nodes added. Zero-padded control matches input dimensionality. Historical molecule provenance retained.'};(ROOT/'validation/gnn_input_manifest.json').write_text(json.dumps(result,indent=2)+'\n');print('Prepared16 matched graphs across8 targets;5 new lysine features;no fitting.')
if __name__=='__main__':main()
