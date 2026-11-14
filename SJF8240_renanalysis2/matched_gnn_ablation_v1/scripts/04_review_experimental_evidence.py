"""Read-only eight-label evidence review;no label changes or differential refit."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'));from contracts import *
BASE=Path('/home/ayushregmibagale/SJF8240_LF_and_PS');AUDIT=Path('/home/ayushregmibagale/10PROTAC-MD/prosettac_robust/existing_result_audit/user_experimental_summary_v2')
UID={'ABL1':'P00519','AXL':'P30530','DDR2':'Q16832','EPHA2':'P29317','EPHB2':'P29323','MAPK14_full':'Q16539','MET':'P08581','RIPK2_full':'O43353'}
def main():
 labels=protocol()['labels'];files=[('LF_saved_new',BASE/'E250916_LF/LIMMA_results_DMSO_vs_SJF8240_Foretinib_new.csv.gz','adj.P.Val'),('PS_light',BASE/'E260106_PS/LIMMA_LIGHT_results_all_contrasts.csv.gz','adj_P_Val'),('PS_heavy_newly_synthesized',BASE/'E260106_PS/LIMMA_HEAVY_results_all_contrasts.csv.gz','adj_P_Val')];rows=[];sources={};rawtables={}
 for dataset,path,pcol in files:
  sources[str(path)]=sha(path);df=pd.read_csv(path);rawtables[dataset]=df
  for target,uid in UID.items():
   # Exact canonical group identity;ambiguous/multiaccession groups are not silently merged.
   selected=df[(df.Protein_group.astype(str)==uid)&(df.Contrast=='SJF8240_vs_DMSO')];require(len(selected)<=1,'Ambiguouscanonical result:'+target)
   if selected.empty:rows.append({'dataset':dataset,'target':target,'uid':uid,'user_label':labels[target],'result_available':False,'log2FC':np.nan,'BH_adjusted_P':np.nan,'interpretation':'Missing/inconclusive;not a negative label'});continue
   r=selected.iloc[0];fc=float(r.logFC);q=float(r[pcol]);interpretation='Observed decrease' if fc<0 and q<.05 else ('Observed increase' if fc>0 and q<.05 else 'No statistically established change in this table;not equivalence')
   rows.append({'dataset':dataset,'target':target,'uid':uid,'user_label':labels[target],'result_available':True,'log2FC':fc,'BH_adjusted_P':q,'interpretation':interpretation,'source_path':str(path),'source_sha256':sources[str(path)]})
 current=pd.DataFrame(rows);prior=pd.read_csv(AUDIT/'named_saved_tables_SJF_vs_DMSO.tsv',sep='\t');lf=current[current.dataset=='LF_saved_new']
 for r in lf.itertuples():
  saved=prior[(prior.poi==r.target)&(prior.display_dataset=='LF_saved_new')];require(len(saved)==1,'PriorLFview missing');require(abs(float(saved.logFC.iloc[0])-r.log2FC)<1e-12 and abs(float(saved.BH_adjusted_P.iloc[0])-r.BH_adjusted_P)<1e-12,'LFtable changed since priorreview')
 allversions=pd.read_csv(AUDIT/'canonical_results_all_versions.tsv',sep='\t');versions=allversions[allversions.poi.isin(labels)].copy();versions=versions.rename(columns={'historical_label':'historical_label_in_prior_review'});versions['current_user_label']=versions.poi.map(labels);sources[str(AUDIT/'canonical_results_all_versions.tsv')]=sha(AUDIT/'canonical_results_all_versions.tsv');versions.to_csv(ROOT/'results/experimental_all_saved_versions.tsv',sep='\t',index=False);current.to_csv(ROOT/'results/experimental_evidence_eight_targets.tsv',sep='\t',index=False)
 text='''# Experimental-label and chemistry review for the matched GNN ablation

The eight training/evaluation labels remain exactly the user's confirmed six-degrader/two-nondegrader map. This read-only review is not used to change labels,features,thresholds or training. It is separate from the new matched ablation.

The user-confirmed LF/PS conditions are RPE1,500nM,12-hour treatment. They differ from the published cell-line endpoint reviewed previously. The PS heavy channel is newly synthesized protein;neither PS channel alone establishes a degradation-rate constant. Authoritative LF processing/version and independent replication provenance remain unresolved. Named_new is a saved version,not proof that it is the final as-run pipeline.

## Version-explicit LF SJF8240-versus-DMSO view

| Target | User label | log2FC | BH-adjustedP | Interpretation |
|---|---|---:|---:|---|
'''
 for r in lf.itertuples():text+=f'| {r.target} | {r.user_label} | {r.log2FC:.3f} | {r.BH_adjusted_P:.4g} | {r.interpretation} |\n'
 text+='''
Five of the six user-positive labels show significant decreases in this named LF table. EPHA2 is inconclusive in this RPE1 table;this does not override the user's positive label or establish a non-degrader. ABL1 is inconclusive near zero;AXL increases versus DMSO. These observations are supportive/contextual evidence,not independent confirmation of each exact binary label or proof of direct ubiquitination. A nonsignificant result is not evidence of equivalence or no degradation.

The complete saved-version and contrast table is retained separately. In particular,AXL can decrease versus foretinib while increasing versus DMSO;the comparator matters. LF versions differ and are not averaged or cherry-picked to match the model. PS light/heavy channels are reported separately. No R pipeline was rerun,no multi-gigabyte raw spectra were reprocessed,andno new wet-lab result was generated. The earlier review's source summaries were checked against the named source files;LFvalues reproduce to1e-12.

The supplied folders also contain later PELSA datasets,but their generic processing READMEs do not establish binding for these eight targets. This review does not claim direct binding validation or analyze those newer assays into model labels. No messages or external submissions were sent.

## Chemistry and interpretation

The retained computational ligand has a historical connectivity/coordinate-stereochemistry provenance that differs from supplier-matched SJF8240. The new ablation preserves that actual historical model and all old results;it cannot validate a supplier-SJF8240 classifier merely by improving AUROC. Geometry repair does not change molecular connectivity. Historical and supplier-matched branches remain separate. Chemical identity audit:../native_catalytic_validation_v2/FINDINGS.md and../../SJF8240_kinome_screen_analysis/repair_v3/chemical_identity_audit.json.

The independent units are eight targets with only two negatives. Assay context,protein identity and chemical provenance limit biological interpretation even if the controlled input ablation improves ranking. Experimental sources were not used to tune the ablation or relabel EPHA2.
'''
 (ROOT/'EXPERIMENTAL_EVIDENCE_REVIEW.md').write_text(text);atomic_json(ROOT/'validation/experimental_review_receipt.json',{'status':'READ_ONLY_REVIEW_COMPLETE','source_sha256':sources,'labels_changed':False,'new_differential_fit':False,'replication_assumptions_not_inferred':True,'LFprior_reproduction_max_error':0.,'review_script_sha256':sha(__file__)});print('Reviewed8labels againstnamedLFandseparatePSchannels;labels unchanged.')
if __name__=='__main__':main()
