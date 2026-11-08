"""Join retained assay workbook by accession; no experimental relabeling."""
from pathlib import Path
import json,sys,numbers
import openpyxl
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,ValidationError
source=ROOT/'sjf8240_review/primary_sources/Bondeson2018_publisher_mmc4.xlsx'
cohort=ROOT/'workflow/config/sjf8240_publication_cohort.json'
identity=ROOT/'existing_result_audit/integration_19_v2/CDK17_source_identity.json'
paths=[source,cohort,identity];hashes={str(p):sha256(p) for p in paths}
c=json.loads(cohort.read_text());w=openpyxl.load_workbook(source,read_only=True,data_only=True)
s=w['Sheet1'];rows=list(s.iter_rows(values_only=True));out=ROOT/'existing_result_audit/assay_accession_v1'
if out.exists():raise ValidationError('Immutable output exists')
out.mkdir();records=[]
for poi,accession in list(c['protein_groups'].items())+[('CDK14_structure_reference','O94921')]:
 hits=[(i,row) for i,row in enumerate(rows,1) if str(row[0]).strip()==accession]
 if len(hits)>1:raise ValidationError('Ambiguous accession join')
 entry={'poi':poi,'accession':accession,'historical_label':c['historical_labels'].get(poi),'row_count':len(hits)}
 if hits:
  number,row=hits[0];low,high=row[6:8]
  numeric=lambda x:isinstance(x,numbers.Real) and not isinstance(x,bool)
  entry.update(sheet='Sheet1',excel_row=number,source_gene_name=row[1],compound1_binding_column_raw=row[3],proteomics_100nM_remaining=low if numeric(low) else None,proteomics_1uM_remaining=high if numeric(high) else None,observed_peptides=row[5],quantified_both=bool(numeric(low) and numeric(high)),source_missing_marker=None if numeric(low) else str(low))
 else:entry['quantified_both']=False
 records.append(entry)
report={'status':'ACCESSION_JOIN_COMPLETE_NO_TRAINING_RELEASE','workbook_sha256':hashes[str(source)],'records':records,'historical_targets':17,'historical_targets_with_both_abundance_measurements':sum(x['quantified_both'] for x in records if x['poi']!='CDK14_structure_reference'),
 'CDK17_structure_correspondence':'The measured CDK17/Q00537 row is distinct from the simulated CDK14/O94921 (workbook PFTK1) row. O94921 has binding-column data but no proteomics quantification; the CDK17 abundance label cannot be transferred to it.',
 'binding_column_interpretation':'Raw workbook column retained without converting to Kd/IC50 or assigning new high/low categories. Binding-column measurements are not protein-abundance outcomes.',
 'source_sha256':hashes,'audit_sha256':sha256(Path(__file__)),'label_or_structure_changed':False,'new_MD_or_GNN_training':False}
if any(sha256(Path(p))!=v for p,v in hashes.items()):raise ValidationError('Source drift')
(out/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'quantified_historical_targets':report['historical_targets_with_both_abundance_measurements'],'CDK17_and_CDK14':[r for r in records if r['poi'].startswith('CDK17') or r['poi'].startswith('CDK14')]},indent=2))
