"""Fixed protein-level summaries and exact rank tests; no fitted predictor."""
from pathlib import Path
import sys,json,itertools,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata,mannwhitneyu
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));from extract import paired_delta,require,verify_source,sha256,atomic_json,CAN,OLD
VARIANTS={'legacy_geometry_legacy_mask':('legacy_min_e2_A','legacy_mask'),'legacy_geometry_true_VHL':('legacy_min_e2_A','true_vhl_mask'),'corrected_geometry_legacy_mask':('corrected_min_e2_A','legacy_mask'),'corrected_geometry_true_VHL':('corrected_min_e2_A','true_vhl_mask')}
def exact_rank(values,labels):
 x=np.asarray(values,float);y=np.asarray(labels,int);require(len(x)==8 and sum(y)==6 and np.isfinite(x).all(),'Eight finite protein summaries required');r=rankdata(x);u=float(r[y==1].sum()-21);null=[];diffnull=[]
 for n in itertools.combinations(range(8),2):
  positive=np.ones(8,bool);positive[list(n)]=False;null.append(float(r[positive].sum()-21));diffnull.append(float(x[positive].mean()-x[~positive].mean()))
 null=np.asarray(null);diffnull=np.asarray(diffnull);diff=float(x[y==1].mean()-x[y==0].mean());p=min(1.,2*min(float((null>=u-1e-12).mean()),float((null<=u+1e-12).mean())))
 return {'rank_U':u,'rank_two_sided_exact_p':p,'positive_mean':float(x[y==1].mean()),'negative_mean':float(x[y==0].mean()),'difference':diff,'abs_difference_permutation_p':float((abs(diffnull)>=abs(diff)-1e-12).mean()),'label_allocations':28}
def summarize(frame):
 rows=[]
 for (variant,target),g in frame.groupby(['variant','poi']):
  values=g.delta_reach_50A.to_numpy();require(np.isfinite(values).all(),'Nonfinite target contrasts');rows.append({'variant':variant,'poi':target,'label':int(g.label.iloc[0]),'lysines':len(g),'positive_lysines':int((values>0).sum()),'negative_lysines':int((values<0).sum()),'zero_lysines':int((values==0).sum()),'positive_after_1e12_tolerance':int((values>1e-12).sum()),'fraction_positive':float((values>0).mean()),'mean_delta':float(values.mean()),'paired_poses':int(g.paired_poses.iloc[0])})
 return pd.DataFrame(rows)
def statistical_table(summary):
 rows=[]
 for variant,g in summary.groupby('variant'):
  if len(g)!=8:continue
  for col in ['mean_delta','fraction_positive']:rows.append({'variant':variant,'summary':col,**exact_rank(g[col],g.label)})
 return pd.DataFrame(rows)
def main():
 verify_source()
 for rel,digest in json.loads((ROOT/'config/reporting_manifest.json').read_text()).items():require(sha256(ROOT/rel)==digest,'Reporting source changed')
 cfg=json.loads((ROOT/'config/protocol.json').read_text());labels=cfg['labels'];completion=json.loads((ROOT/'analysis/completion.json').read_text());require(completion['cases']==141 and completion['frames']==28341,'Incomplete full-frame cohort')
 oldposepath=OLD/'analysis/per_lysine_e2_all_pois.tsv';oldsummarypath=OLD/'analysis/per_lysine_e2_summary.tsv';oldpose=pd.read_csv(oldposepath,sep='\t');oldsummary=pd.read_csv(oldsummarypath,sep='\t');perpose=[];coverage=[];packets={};legacy_errors=[];coarse_errors=[]
 for record in completion['records']:
  target,pose=record['poi'],record['pose'];p=ROOT/'analysis/poses'/target/pose/'packet.npz';require(sha256(p)==record['packet_sha256'],'Packet drift')
  with np.load(p) as z:packet={k:z[k].copy() for k in z.files}
  require(packet['times_ns'].shape==(201,) and len(set(packet['lys_resids']))==len(packet['lys_resids']),'Frame/residue mismatch');packets[(target,pose)]=packet
  coarse=CAN/'mechanistic_v3/analysis/lysines/poses'/target/pose/'packet.npz'
  with np.load(coarse) as z:
   require(np.array_equal(z['lys_resids'],packet['lys_resids']) and np.array_equal(z['bridge'],packet['true_vhl_mask'][::5]),'Coarse residue/mask mismatch');error=float(abs(z['min_e2_A']-packet['corrected_min_e2_A'][::5]).max());require(error<.002,'Coarse distance replay mismatch');coarse_errors.append(error)
  matched=all(packet[k].any() and (~packet[k]).any() for k in ['legacy_mask','true_vhl_mask'])
  for variant,(distance,mask) in VARIANTS.items():
   b=packet[mask];delta=paired_delta(packet[distance],b);coverage.append({'poi':target,'pose':pose,'variant':variant,'bridge_frames':int(b.sum()),'nonbridge_frames':int((~b).sum()),'paired':delta is not None,'matched_support':matched})
   if delta is not None:
    for resid,val in zip(packet['lys_resids'],delta):perpose.append({'poi':target,'pose':pose,'label':labels[target],'variant':variant,'resid':int(resid),'delta_reach_50A':float(val),'matched_support':matched})
  old=oldpose[(oldpose.poi==target)&(oldpose.pose_id==pose)].set_index('resid');require(set(old.index)==set(packet['lys_resids']),'Legacy residue/case scope mismatch');delta=paired_delta(packet['legacy_min_e2_A'],packet['legacy_mask']);expected=old.loc[packet['lys_resids']].delta_reach_50A.to_numpy()
  if delta is None:require(np.isnan(expected).all(),'Legacy state eligibility mismatch')
  else:require(np.isfinite(expected).all(),'Legacy missing contrast');error=float(abs(delta-expected).max());require(error<=.0002,'Legacy per-pose reproduction failed '+target+'/'+pose);legacy_errors.append(error)
 perpose=pd.DataFrame(perpose);cover=pd.DataFrame(coverage);perlys=perpose.groupby(['variant','poi','label','resid'],as_index=False).agg(delta_reach_50A=('delta_reach_50A','mean'),paired_poses=('pose','nunique'));summary=summarize(perlys);stats=statistical_table(summary);matched=perpose[perpose.matched_support].groupby(['variant','poi','label','resid'],as_index=False).agg(delta_reach_50A=('delta_reach_50A','mean'),paired_poses=('pose','nunique'));matchedsummary=summarize(matched);matchedstats=statistical_table(matchedsummary)
 require(len(summary)==32 and all(summary.groupby('variant').size()==8),'Incomplete target coverage');results=ROOT/'results';results.mkdir(exist_ok=True)
 for name,table in [('per_pose_per_lysine.tsv',perpose),('state_coverage.tsv',cover),('per_lysine.tsv',perlys),('target_summary.tsv',summary),('protein_statistics.tsv',stats),('matched_support_target_summary.tsv',matchedsummary),('matched_support_statistics.tsv',matchedstats)]:table.to_csv(results/name,sep='\t',index=False)
 # Direct comparison with published rounded legacy summaries.
 a=perlys[perlys.variant=='legacy_geometry_legacy_mask'].merge(oldsummary[['poi','resid','delta_reach_50A']],on=['poi','resid'],suffixes=('_replayed','_published'),validate='one_to_one');a['absolute_error']=abs(a.delta_reach_50A_replayed-a.delta_reach_50A_published);require(len(a)==132 and a.absolute_error.max()<=.0002,'Published summary failed reproduction');a.to_csv(results/'legacy_reproduction.tsv',sep='\t',index=False)
 current=summary[summary.variant=='corrected_geometry_true_VHL'];legacy=summary[summary.variant=='legacy_geometry_legacy_mask'];comparison=legacy.merge(current,on=['poi','label'],suffixes=('_legacy','_corrected'));comparison.to_csv(results/'legacy_vs_corrected.tsv',sep='\t',index=False)
 # Exactly two primary protein summaries: report Holm correction as well.
 primary=stats[stats.variant=='corrected_geometry_true_VHL'].copy();p=primary.rank_two_sided_exact_p.to_numpy();order=np.argsort(p);adjust=np.empty(2);adjust[order[0]]=min(1,2*p[order[0]]);adjust[order[1]]=max(adjust[order[0]],p[order[1]]);primary['Holm_adjusted_two_summary_p']=adjust;primary.to_csv(results/'primary_statistics.tsv',sep='\t',index=False)
 report='# Full-frame true-VHL per-lysine reach50 comparison\n\nAll141 retained cases and201 frames per case (28,341 total) were recomputed. Same50 A threshold and eight current-labelled proteins; ABL1/AXL are negatives, EPHA2 is positive. SIK3/RPS6KB2 remain unknown and excluded. No GNN fitting or new MD.\n\n'
 report+='| Target | Legacy positive lysines | Corrected true-VHL positive lysines | Legacy mean delta | Corrected mean delta | Corrected paired poses |\n|---|---:|---:|---:|---:|---:|\n'
 for row in comparison.to_dict('records'):report+=f'| {row["poi"]} | {row["positive_lysines_legacy"]}/{row["lysines_legacy"]} | {row["positive_lysines_corrected"]}/{row["lysines_corrected"]} | {row["mean_delta_legacy"]:+.5f} | {row["mean_delta_corrected"]:+.5f} | {row["paired_poses_corrected"]} |\n'
 report+='\n## Protein-level evidence\n\n'
 for row in primary.to_dict('records'):report+=f'- {row["summary"]}: exact two-sided rank p={row["rank_two_sided_exact_p"]:.6f}; Holm-adjusted p={row["Holm_adjusted_two_summary_p"]:.6f}; mean degrader-minus-nondegrader difference={row["difference"]:+.6f}.\n'
 report+='''
The independent units are8 proteins, not132 lysines or28,341 frames. Exact rank p-values enumerate all28 allocations of two negatives, with averaged ranks for ties and twice the smaller inclusive tail. Difference-statistic permutation p-values are separately labelled in the tables. These are retrospective, observational comparisons, not external validation; labels may be confounded by protein identity. A small p-value alone would not establish a degradation predictor. Positive means strictly delta>0; a1e-12 numerical-zero sensitivity count is retained.

## Meaning and controls

For each lysine/frame, take the nearest modeled E2 point distance. Reach50 is the fraction of bridge/nonbridge FRAMES with that distance<50 A. Delta is bridge minus nonbridge, calculated within each pose and then equally averaged over poses containing both states. It is not the fraction of E2 points and not evidence of solvent burial or ubiquitin-transfer chemistry. The reference file's663 concatenated rows contain51 unique coordinates; duplicate rows do not affect a nearest-point distance.

All four geometry/mask combinations use identical201-frame grids. The corrected geometry uses ElonginC fitting and ligand-connected whole-component imaging. True-VHL bridging uses periodic ligand/VHL and ligand/target heavy-atom contacts<5 A. The legacy mask uses the original metrics.csv convention, whose VHL label actually selected ELOB. The geometry factor bundles reference/image convention; the mask factor bundles identity/atom-selection policy. We do not attribute a change to one subcomponent alone.

The primary comparisons use every eligible paired-state pose for each mask. Matched-support sensitivity uses only poses eligible under both masks, preventing different pose subsets from driving that particular comparison. Coverage tables retain all poses, including those with no valid contrast. No unsupported contrast is imputed as zero. All8 targets must be present before protein-level tests are computed. Distances at the old41 sampled times are checked against the saved mechanistic_v3 packets.

The legacy per-pose and published per-lysine tables are independently checked to their rounding precision. Original files, labels, old GNN scores and older analyses remain unchanged. This result tests the specific reach50 hypothesis; it does not retrain or retrospectively validate the0.750 temporal GNN, which did not explicitly receive per-lysine delta50.

## Files

results/target_summary.tsv and legacy_vs_corrected.tsv show all targets and conventions. results/per_lysine.tsv and per_pose_per_lysine.tsv retain each lysine. results/protein_statistics.tsv and primary_statistics.tsv contain complete exact tests; matched_support_*.tsv contains pose-matched sensitivity. results/state_coverage.tsv records state support. validation/ records frozen-source, legacy/coarse reproduction and independent unwrap replay. Local analysis/poses/ contains full distance packets with source hashes.
'''
 (ROOT/'RESULTS.md').write_text(report)
 atomic_json(ROOT/'validation/reproduction.json',{'status':'PASS','cases':141,'frames':28341,'per_lysine_rows':132,'max_legacy_per_pose_delta_error':max(legacy_errors),'max_published_legacy_delta_error':float(a.absolute_error.max()),'max_coarse_distance_error_A':max(coarse_errors),'oldpose_sha256':sha256(oldposepath),'oldsummary_sha256':sha256(oldsummarypath),'source_manifest_sha256':verify_source()})
 print(report)
if __name__=='__main__':main()
