"""Read-only accession/gene audit of user-supplied LF and PS experiments."""
from pathlib import Path
import re,json,sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,ValidationError
source_root=Path('/home/ayushregmibagale/SJF8240_LF_and_PS')
cohort_path=ROOT/'workflow/config/sjf8240_publication_cohort.json';cohort=json.loads(cohort_path.read_text())
targets=dict(cohort['protein_groups']);targets.update(SIK3='Q9Y2K2',RPS6KB2='Q9UBS0',CDK14_structure_reference='O94921')
# Prediction-target accessions must agree with the retained canonical sequence record.
sequence_path=ROOT/'existing_result_audit/scientific_signoff_v1/all19_sequences.json'
# Retain the explicit accession mapping as an audit input, not an experimental identity claim.
out=ROOT/'existing_result_audit/user_experimental_evidence_v1'
if out.exists():raise ValidationError('Immutable audit output exists')
out.mkdir();(out/'sources').mkdir()
sources={str(cohort_path):sha256(cohort_path)};inventory=[];statistics=[];matches=[];quants=[];samples=[]
def tokens(value):return set(re.split(r'[;,\s]+',str(value).strip())) if pd.notna(value) else set()
def canonical(value):return {x.split('-')[0] for x in tokens(value)}
def finite(x):return bool(pd.notna(x) and np.isfinite(float(x)))
def matching(df,poi,accession):
 genes=df['Genes'].map(tokens);accessions=df['Protein_group'].map(canonical)
 gene=poi.removesuffix('_full') if poi!='CDK14_structure_reference' else 'CDK14'
 return accessions.map(lambda x:accession in x)|genes.map(lambda x:gene in x)
for exp in ['E250916_LF','E260106_PS']:
 root=source_root/exp
 for p in sorted(root.iterdir()):
  if not p.is_file() or p.name.startswith('.'):continue
  inventory.append({'experiment':exp,'path':str(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns,'hashed':p.name.endswith('.csv.gz') or p.suffix=='.R'})
  if p.suffix=='.R':
   sources[str(p)]=sha256(p)
   # Freeze only primary quantification/differential-analysis scripts, never execute them.
   if p.name in ['DAAPROTACvsDMSOandForetinbvsDMSO.R','Diann_protein_quant_LF_E250916_new.R','DIANN_protein_quant.R','DAA_LIGHT_for_SJF8240_vs_DMSO_Foretinib_vs_DMSO_FDR.R','DAA_HEAVY_for_SJF8240_vs_DMSO_FORETINIB_vs_DMSO_p_value.R']:
    (out/'sources'/f'{exp}__{p.name}').write_bytes(p.read_bytes())
  if not p.name.endswith('.csv.gz'):continue
  sources[str(p)]=sha256(p);df=pd.read_csv(p)
  stats={'experiment':exp,'file':p.name,'rows':len(df),'columns':list(df.columns)}
  if 'Contrast' in df:stats['contrast_counts']={str(k):int(v) for k,v in df.Contrast.value_counts(dropna=False).items()}
  statistics.append(stats)
  if not {'Genes','Protein_group'}<=set(df):continue
  if p.name.startswith('LIMMA'):
   if not {'Contrast','logFC'}<=set(df):continue
   fdr='adj_P_Val' if 'adj_P_Val' in df else 'adj.P.Val'
   for poi,accession in targets.items():
    for _,row in df.loc[matching(df,poi,accession)].iterrows():
     raw=canonical(row.Protein_group);gene=poi.removesuffix('_full') if poi!='CDK14_structure_reference' else 'CDK14'
     matches.append({'experiment':exp,'file':p.name,'poi':poi,'expected_accession':accession,'historical_label':cohort['historical_labels'].get(poi),'Protein_group':row.Protein_group,'Genes':row.Genes,'accession_match':accession in raw,'gene_match':gene in tokens(row.Genes),'protein_group_canonical_accessions':';'.join(sorted(raw)),'multi_accession_group':len(raw)>1,'Contrast':row.Contrast,'logFC':float(row.logFC) if finite(row.logFC) else None,'P_value':float(row['P.Value']) if finite(row.get('P.Value')) else None,'BH_adjusted_P':float(row[fdr]) if finite(row.get(fdr)) else None})
  if p.name.startswith('protein_quant'):
   mapping=df[['Run','Sample','Condition','Replicate']].drop_duplicates()
   for rec in mapping.to_dict('records'):samples.append(dict(experiment=exp,file=p.name,**rec))
   for poi,accession in targets.items():
    hit=df.loc[matching(df,poi,accession)]
    cols=[c for c in ['LFQ','LFQ_H','LFQ_L','N_precursors','N_precursors_H','N_precursors_L','N_precursors_proteotypic','N_precursors_proteotypic_H','N_precursors_proteotypic_L'] if c in df]
    for _,row in hit.iterrows():
     rec={k:row[k] for k in ['Run','Sample','Condition','Replicate','Protein_group','Genes']};rec.update(experiment=exp,file=p.name,poi=poi,expected_accession=accession,accession_match=accession in canonical(row.Protein_group),multi_accession_group=len(canonical(row.Protein_group))>1)
     rec.update({c:float(row[c]) if finite(row[c]) else None for c in cols});quants.append(rec)
if any(sha256(Path(p))!=v for p,v in sources.items()):raise ValidationError('Source changed during read-only audit')
pd.DataFrame(matches).to_csv(out/'target_differential_results.tsv',sep='\t',index=False)
pd.DataFrame(quants).to_csv(out/'target_quantification.tsv',sep='\t',index=False)
pd.DataFrame(samples).to_csv(out/'sample_mapping.tsv',sep='\t',index=False)
(out/'table_inventory.json').write_text(json.dumps(statistics,indent=2)+'\n')
(out/'file_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
report={'status':'USER_EXPERIMENTAL_EVIDENCE_EXTRACTED_NOT_RELABELLED','conditions_user_confirmed':{'cell_line':'RPE1','dose_nM':500,'treatment_hours':12,'applies_to':['E250916_LF','E260106_PS'],'PS_heavy':'newly synthesized proteins','source':'User message in this session','pulse_duration_hours':None,'foretinib_dose_separately_confirmed':False},'targets':targets,'historical_labels_unchanged':True,'source_sha256':sources,'audit_sha256':sha256(Path(__file__)),'differential_target_rows':len(matches),'quantification_target_rows':len(quants),'raw_reports_reprocessed':False,'limma_refitted':False,'new_MD_or_GNN_training':False,'limits':['Conditions are user-confirmed metadata, not independently recovered laboratory notebook entries.','PS pulse timing and isotope incorporation require protocol confirmation before a kinetic interpretation.','Saved result tables have multiple versions and output-overwriting scripts; no exact as-run receipt is inferred.','PS heavy measures newly synthesized pool under the confirmed labeling design; reduced heavy intensity alone does not establish a change in degradation rate.','Non-significant, unquantified and ambiguous protein-group observations are not experimental negative labels.','Original LF/PS folders and historical PROTAC-MD outputs remain unchanged.']}
report['output_sha256']={p.name:sha256(p) for p in out.iterdir() if p.is_file()}
(out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'target_result_rows':len(matches),'target_quant_rows':len(quants),'sources_hashed':len(sources)},indent=2))
