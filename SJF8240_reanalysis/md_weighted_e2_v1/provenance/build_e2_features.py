"""Reproduce audited ligand-imaged E2 scalars from existing frames; no simulations."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
BASE=Path('/home/ayushregmibagale/10PROTAC-MD');REANALYSIS=BASE/'SJF8240_reanalysis';E2=REANALYSIS/'analysis/e2/ligand_connected_e2_cohort_v1';CONTACT=BASE/'prosettac_robust/existing_result_audit/contact_core_v1/poses';ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 labels=json.loads((REANALYSIS/'config/current_md_labels.json').read_text())['labels'];dynamic=json.loads((REANALYSIS/'analysis/dynamic/completion.json').read_text())['records'];cases=[];rows=[]
 for poi in sorted(labels):
  target=[]
  for record in dynamic:
   if record['poi']!=poi:continue
   folder=E2/'poses'/poi/record['pose_id'];summary_path=folder/'summary.json';frame_path=folder/'frames.tsv';contact_path=CONTACT/poi/record['pose_id']/'contacts_A.tsv';summary=json.loads(summary_path.read_text())
   assert sha(frame_path)==summary['frames_sha256']
   source_hashes=summary['source_input_sha256']
   for path,digest in record['source_sha256'].items():assert source_hashes[path]==digest,'MD/E2 source lineage mismatch'
   expected_contact=next(d for path,d in source_hashes.items() if path.endswith('/contacts_A.tsv'));assert sha(contact_path)==expected_contact
   f=pd.read_csv(frame_path,sep='\t');c=pd.read_csv(contact_path,sep='\t');assert len(f)==len(c)==201
   assert np.allclose(f.time_ns,c.time_ns,atol=5.1e-6,rtol=0) and np.allclose(f.time_ns,np.linspace(0,20,201),atol=5.1e-6,rtol=0)
   bridge=(c.true_VHL_heavy_A.to_numpy()<5.)&(c.POI_heavy_A.to_numpy()<5.);minimum=f.ligand_image_min_E2_A.to_numpy();reach=f.ligand_image_reach_80A.to_numpy();assert np.isfinite(minimum).all() and np.isfinite(reach).all() and ((reach>=0)&(reach<=1)).all()
   mean=float(minimum.mean());delta=float(reach[bridge].mean()-reach[~bridge].mean()) if bridge.any() and (~bridge).any() else None
   audited=summary['comparisons']['ligand_image'];state=audited['conditional']['true_VHL_POI_heavy']
   assert abs(mean-audited['mean_min_E2_A'])<1e-10 and state['bridge_frames']==int(bridge.sum())
   assert (delta is None and state['delta_80A_paired'] is None) or (delta is not None and abs(delta-state['delta_80A_paired'])<1e-12)
   item={'poi':poi,'pose_id':record['pose_id'],'mean_min_E2_A':mean,'delta_80A_paired':delta,'bridge_frames':int(bridge.sum()),'nonbridge_frames':int((~bridge).sum()),'frames_sha256':sha(frame_path),'contacts_sha256':sha(contact_path),'summary_sha256':sha(summary_path),'summary_path':str(summary_path),'frames_path':str(frame_path),'contacts_path':str(contact_path)};target.append(item);cases.append(item)
  paired=[x['delta_80A_paired'] for x in target if x['delta_80A_paired'] is not None];assert target and paired
  rows.append({'poi':poi,'e2_mean_min_dist_A':float(np.mean([x['mean_min_E2_A'] for x in target])),'delta_reach_80A':float(np.mean(paired)),'MD_poses':len(target),'paired_state_poses':len(paired)})
 table=pd.DataFrame(rows)
 reported_means=pd.read_csv(REANALYSIS/'analysis/e2/coherent_e2_publication_v1/per_target_imaging.tsv',sep='\t').set_index('poi');reported_deltas=pd.read_csv(REANALYSIS/'analysis/e2/coherent_e2_publication_v1/paired_state_sensitivity.tsv',sep='\t');reported_deltas=reported_deltas[(reported_deltas.image_convention=='ligand_image')&(reported_deltas.state_mask=='true_VHL_POI_heavy')&(reported_deltas.threshold_A==80)].set_index('poi')
 for row in rows:
  assert abs(row['e2_mean_min_dist_A']-reported_means.loc[row['poi'],'ligand_image_mean_min_E2_A'])<1e-10
  assert abs(row['delta_reach_80A']-reported_deltas.loc[row['poi'],'mean_within_pose_bridge_minus_nonbridge'])<1e-12
 table.to_csv(ROOT/'e2_features.tsv',sep='\t',index=False)
 manifest={'status':'E2_FEATURES_INDEPENDENTLY_REPRODUCED_FOR_EXPLORATORY_GNN','feature_columns':['e2_mean_min_dist_A','delta_reach_80A'],'definitions':{'e2_mean_min_dist_A':'Equal-pose mean of time-averaged minimum distance over modeled E2 points and POI lysine NZ atoms, ligand-connected coherent image, Angstrom; unconditional across201frames.','delta_reach_80A':'Equal mean of within-pose (bridge minus nonbridge) fractions of modeled E2 points with a lysine NZ within strict80A; use only poses with both states. Bridge is simultaneous true-VHL/ligand and POI/ligand heavy distances <5A.'},'image_convention':'ligand_image','bridge_mask':'true_VHL_POI_heavy','targets':8,'poses':len(cases),'frames':201*len(cases),'all_targets_have_finite_two_feature_values':bool(np.isfinite(table[['e2_mean_min_dist_A','delta_reach_80A']]).all().all()),'both_graph_roles_use_identical_E2_features':True,'feature_table_sha256':sha(ROOT/'e2_features.tsv'),'builder_sha256':sha(__file__),'per_case':cases,'old_ABL1_delta_minus_00466_not_used':True,'missing_conditional_states':'Undefined per-pose contrasts remain missing, never0; target contrasts average paired-state poses only, with coverage counts retained.','interpretation':'Modeled geometry proxy; not lysine burial, E2 catalysis, degradation probability or a scientifically certified input release. User explicitly requested exploratory E2-GNN fitting.','prediction_only_extensions':'No corrected audited E2 features for SIK3/RPS6KB2 here; E2 model must not fabricate or zero-fill them.','new_MD_or_docking':False}
 (ROOT/'feature_provenance.json').write_text(json.dumps(manifest,indent=2)+'\n');print(table.to_string(index=False))
if __name__=='__main__':main()
