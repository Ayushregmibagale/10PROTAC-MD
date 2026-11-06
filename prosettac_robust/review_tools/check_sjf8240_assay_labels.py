#!/usr/bin/env python3
"""Audit the historical labels against a named published assay, without relabeling."""
import csv,json,hashlib
from pathlib import Path
import pandas as pd
r=Path(__file__).resolve().parents[1];review=r/'sjf8240_review';source=review/'primary_sources/Bondeson2018_publisher_mmc4.xlsx';d=pd.read_excel(source,header=None)
cohort=json.loads((r/'workflow/config/sjf8240_publication_cohort.json').read_text())['historical_labels'];rows=[]
for poi,label in cohort.items():
 gene=poi.removesuffix('_full');hit=d[d[1].astype(str).str.strip().eq(gene)];row={'poi':poi,'gene':gene,'historical_label':label,'endpoint':'VHL compound 1, MDA-MB-231, 24 h; 100 nM or 1 uM; primary proteomics endpoint'}
 if len(hit)!=1:row.update(status='NOT_FOUND_IN_THIS_ASSAY_TABLE',derived_proteomics_label=None)
 else:
  x=hit.iloc[0];vals=[pd.to_numeric(x[i],errors='coerce') for i in [6,7,10,11]]
  if not all(pd.notna(a) for a in vals):row.update(status='INSUFFICIENT_QUANTIFICATION',derived_proteomics_label=None)
  else:
   low,high,foretinib,epimer=map(float,vals);loss=min(low,high)<=.8;control_ok=foretinib>=.9 and epimer>=.9;derived=int(loss) if (not loss or control_ok) else None
   row.update(uniprot=str(x[0]),observed_peptides=int(x[5]) if pd.notna(x[5]) else None,VHL_100nM_remaining=low,VHL_1uM_remaining=high,foretinib_1uM_remaining=foretinib,epimer_1uM_remaining=epimer,derived_proteomics_label=derived,status='AGREES_FOR_THIS_ASSAY' if derived==label else 'ASSAY_DISCORDANCE_REQUIRES_REVIEW')
   if gene=='ABL1':row.update(status='PROTEOMICS_WESTERN_DISCORDANCE_IN_PAPER',note='The paper reports no compound-1 degradation by western blot despite the proteomics decrease. Endpoint priority must be frozen independently of predictions.')
 rows.append(row)
summary={'status':'ASSAY_SPECIFIC_LABEL_REVIEW_NOT_A_TRAINING_RELEASE','source_url':'https://ars.els-cdn.com/content/image/1-s2.0-S2451945617303549-mmc4.xlsx','article_url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5777153/','DOI':'10.1016/j.chembiol.2017.09.010','workbook_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'historical_cohort_count':len(rows),'source_rows_found':sum('uniprot'in x for x in rows),'historical_labels_unchanged':True,'criteria':'At least 20% lower abundance at either concentration, excluding reductions accompanied by >10% loss in foretinib or VHL-epimer controls; endpoint-specific and not a universal biological label','rows':rows,'prohibited':'Never derive experimental labels from docking success, MD bridging, graph availability or classifier predictions.'}
(review/'primary_assay_label_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print([(x['poi'],x['status'],x['derived_proteomics_label']) for x in rows])
