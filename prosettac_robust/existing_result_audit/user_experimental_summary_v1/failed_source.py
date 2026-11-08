"""Version-explicit summaries of the frozen targeted LF/PS extraction."""
from pathlib import Path
import json,sys
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,ValidationError
src=ROOT/'existing_result_audit/user_experimental_evidence_v1';out=ROOT/'existing_result_audit/user_experimental_summary_v1'
if out.exists():raise ValidationError('Output already exists')
out.mkdir();d=pd.read_csv(src/'target_differential_results.tsv',sep='\t');q=pd.read_csv(src/'target_quantification.tsv',sep='\t')
# Keep all versions and contrasts; do not collapse independent gene-only groups into the canonical accession.
canonical=d[d.accession_match & ~d.multi_accession_group].copy()
canonical['dataset']=canonical.experiment+' / '+canonical.file
key=['dataset','poi','Contrast']
if canonical.duplicated(key).any():raise ValidationError('Ambiguous canonical result')
canonical['down_BH005_and_log2FC_le_minus03']=(canonical.BH_adjusted_P<.05)&(canonical.logFC<=-.3)
canonical['not_significant_is_negative_label']=False
canonical.to_csv(out/'canonical_results_all_versions.tsv',sep='\t',index=False)
selected=[('LF_saved_new','E250916_LF','LIMMA_results_DMSO_vs_SJF8240_Foretinib_new.csv.gz'),('PS_LIGHT_saved','E260106_PS','LIMMA_LIGHT_results_all_contrasts.csv.gz'),('PS_HEAVY_saved','E260106_PS','LIMMA_HEAVY_results_all_contrasts.csv.gz')]
rows=[]
for label,exp,file in selected:
 t=canonical[(canonical.experiment==exp)&(canonical.file==file)&(canonical.Contrast=='SJF8240_vs_DMSO')]
 for poi in json.loads((src/'summary.json').read_text())['targets']:
  hits=t[t.poi==poi];row={'display_dataset':label,'poi':poi,'result_available':len(hits)==1}
  if len(hits):row.update(hits.iloc[0][['logFC','BH_adjusted_P','down_BH005_and_log2FC_le_minus03','Protein_group']].to_dict())
  rows.append(row)
pd.DataFrame(rows).to_csv(out/'named_saved_tables_SJF_vs_DMSO.tsv',sep='\t',index=False)
qual=[]
for (exp,file,poi,condition),t in q[q.accession_match & ~q.multi_accession_group].groupby(['experiment','file','poi','Condition']):
 for channel,precursor in [('LFQ','N_precursors'),('LFQ_L','N_precursors_L'),('LFQ_H','N_precursors_H')]:
  if channel not in t or t[channel].isna().all():continue
  vals=pd.to_numeric(t[channel],errors='coerce');valid=vals.notna()&np.isfinite(vals)&(vals>0)
  accepted=valid&(t[precursor]>=2) if precursor in t else pd.Series(False,index=t.index)
  qual.append({'experiment':exp,'file':file,'poi':poi,'Condition':condition,'channel':channel,'rows':len(t),'samples':t.Sample.nunique(),'valid_positive_samples':t.loc[valid,'Sample'].nunique(),'positive_with_at_least_two_precursors_samples':t.loc[accepted,'Sample'].nunique()})
pd.DataFrame(qual).to_csv(out/'target_sample_support.tsv',sep='\t',index=False)
summary=json.loads((src/'summary.json').read_text());paths=summary['source_sha256'];lf=str(Path('/home/ayushregmibagale/SJF8240_LF_and_PS/E250916_LF/LIMMA_results_DMSO_vs_SJF8240_Foretinib.csv.gz'));ps=lf.replace('E250916_LF','E260106_PS')
report={'status':'VERSION_EXPLICIT_DESCRIPTIVE_SUMMARY_NOT_LABEL_RELEASE','LF_old_result_copy_in_PS_folder_identical_sha256':paths[lf]==paths[ps],'selected_display_tables_not_declared_authoritative':selected,'canonical_available_counts':{label:sum(x['result_available'] for x in rows if x['display_dataset']==label and x['poi']!='CDK14_structure_reference') for label,_,_ in selected},'display_threshold':{'BH_adjusted_P_lt':.05,'log2FC_le':-.3,'purpose':'Explicit descriptive annotation, not chosen training labels'},'source_sha256':{str(p):sha256(p) for p in [src/'summary.json',src/'target_differential_results.tsv',src/'target_quantification.tsv']},'audit_sha256':sha256(Path(__file__)),'as_run_LF_version_resolved':False,'original_files_modified':False,'labels_released':False,'limits':['All LF versions remain visible; selecting a favorable version by performance is prohibited.','A matching accession is not proof of unique peptide specificity; raw precursor validation is not completed.','Saved result FDR thresholds vary across scripts (raw P<0.05 versus BH 0.20); descriptive table uses its stated threshold only.','LF script comment says two samples per condition while its filter is total rows >=4; this does not enforce per-condition support.','A copied LF result inside PS folder must not be counted as independent experimental replication.']}
report['output_sha256']={p.name:sha256(p) for p in out.iterdir() if p.is_file()}
(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['canonical_available_counts','LF_old_result_copy_in_PS_folder_identical_sha256']},indent=2))
