from pathlib import Path
import json
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def main():
 rows=[];windows=[];perlys=[]
 labels=json.loads((ROOT/'config/protocol.json').read_text())['labels']
 for target,label in labels.items():
  packets=[dict(np.load(p)) for p in sorted((ROOT/'analysis/lysines/poses'/target).glob('*/packet.npz'))];joint=[];sasa=[];bs=[];nb=[];ls=[];paired=0
  for z in packets:
   exposed=z['complex_sasa_A2']>5;b=z['bridge'];joint.append(((z['min_e2_A']<80)&exposed&b[:,None]).mean());sasa.append(z['complex_sasa_A2'].mean())
   if b.any():bs.append(exposed[b].mean());ls.append(z['complex_sasa_A2'][b].mean(0))
   if (~b).any():nb.append(exposed[~b].mean())
   if b.any() and (~b).any():paired+=1
  for j,resid in enumerate(packets[0]['lys_resids']):perlys.append({'target':target,'label':label,'resid':int(resid),'bridge_mean_NZ_sasa_A2':float(np.mean(ls,axis=0)[j]) if ls else np.nan,'poses_with_bridge':len(ls)})
  rows.append({'target':target,'label':label,'poses':len(packets),'NZ_sasa_pose_mean_A2':float(np.mean(sasa)),'NZ_sasa_pose_SD_A2':float(np.std(sasa)),'joint_fraction_pose_mean':float(np.mean(joint)),'joint_fraction_pose_SD':float(np.std(joint)),'bridge_exposed_fraction':float(np.mean(bs)) if bs else np.nan,'nonbridge_exposed_fraction':float(np.mean(nb)) if nb else np.nan,'bridge_supported_poses':len(bs),'nonbridge_supported_poses':len(nb),'paired_state_supported_poses':paired,'lysines_mean_bridge_sasa_below5':int((np.mean(ls,axis=0)<5).sum()) if ls else None})
  for w in range(4):
   vals=[]
   for z in packets:
    m=(z['times']>=w*5)&((z['times']<(w+1)*5) if w<3 else z['times']<=20);vals.append([z['complex_sasa_A2'][m].mean(),((z['complex_sasa_A2'][m]>5)&(z['min_e2_A'][m]<80)&z['bridge'][m,None]).mean()])
   a=np.array(vals);windows.append({'target':target,'label':label,'window':w,'start_ns':w*5,'end_ns':(w+1)*5,'NZ_sasa_A2':a[:,0].mean(),'joint_fraction':a[:,1].mean()})
 dest=ROOT/'results';pd.DataFrame(rows).to_csv(dest/'pose_stability_and_bridge_exposure.tsv',sep='\t',index=False);pd.DataFrame(windows).to_csv(dest/'time_window_summary.tsv',sep='\t',index=False);pd.DataFrame(perlys).to_csv(dest/'per_lysine_bridge_SASA.tsv',sep='\t',index=False)
 out='''# Per-lysine geometric findings

141 existing poses; 41 sampled frames each = 5,781 frames. Heavy-atom solvent occlusion, 1.4 A water probe, 960 sphere points. All original inputs unchanged. Time windows are 0–5, 5–10, 10–15 and 15–20 ns, last including 20. Each pose has equal weight. E2 quantities reproduce previously audited ligand-connected images (maximum absolute distance replay error recorded per pose).

These measurements do not support the earlier claim that bridging buries all 18 ABL1 lysines. That claim incorrectly interpreted an E2 reach-fraction change as burial. NZ solvent-accessible area now measures exposure directly. The geometries remain highly overlapping across the two labelled classes. Exposed means NZ SASA >5 square angstroms; reachable means minimum distance to a modeled E2 point <80 angstroms. Thresholds were fixed before inspecting per-target exposure results or fitting the enriched model and are broad geometric proxies, not a ubiquitination reaction criterion.

'''
 for r in rows:out+=f'- {r["target"]}: mean NZ SASA {r["NZ_sasa_pose_mean_A2"]:.2f} A² (pose SD {r["NZ_sasa_pose_SD_A2"]:.2f}); exposed fraction in bridge frames {r["bridge_exposed_fraction"]:.3f}; joint exposed/reachable/bridge occupancy {r["joint_fraction_pose_mean"]:.3f}.\n'
 out+='''
The accompanying tables retain every lysine, pose variability, bridge-state coverage and four time windows. Bridge and nonbridge descriptive averages can have different contributing poses; their difference is not a paired causal estimate. No feature threshold or direction was selected to maximize class separation; no new conventional classifier was fitted.

SASA represents solvent access and steric occlusion by retained complex heavy atoms. It does not model approach of the full E2–ubiquitin complex, catalytic orientation, CRL2 motion or ubiquitin transfer. Hydrogen atoms are omitted from occluders consistently. C-alpha residue IDs are simulation numbering, not verified ubiquitination-site annotations. Only two independently labelled negatives remain, regardless of pose/frame count.
'''
 (dest/'MECHANISTIC_FINDINGS.md').write_text(out);print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()
