"""Descriptive fold-prior audit, not an additional predictive model."""
from pathlib import Path
import json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=json.loads((ROOT/'config/protocol.json').read_text());rows=[]
 for role in ['elob','vhl']:
  for held in p['cohort']:
   r=json.loads((ROOT/'analysis'/role/f'fold_{held}.json').read_text());counts=r['training_class_counts'];assert held not in r['training_targets'] and len(r['training_targets'])==9
   rows.append({'role':role,'held_out':held,'evaluation_label':p['labels'][held],'label_is_hypothetical':held in p['hypothetical_targets'],'training_positives':counts['positive'],'training_negatives':counts['negative'],'training_positive_fraction':counts['positive']/9,'held_out_score':r['prediction']['score'],'predicted_class':r['prediction']['pred']})
 out=ROOT/'results/fold_class_priors.tsv';pd.DataFrame(rows).to_csv(out,sep='\t',index=False)
 text='''# Interpretation of the v2 cohort test

The unweighted-loss and full target-LOPO procedure are intentionally unchanged from v1. Holding out a positive leaves five positives and four negatives (55.6% positive) in training. Holding out a negative leaves six positives and three negatives (66.7% positive). A predictor dominated by training prevalence could therefore rank negative-held-out targets above positive-held-out targets. This is a conditional explanation, not proof that the fitted GNN uses only prevalence. The saved fold_class_priors.tsv lets this be inspected directly without another fitted model or any score adjustment.

Moving from 6:2 to 6:4 does not remove every imbalance-related issue: positives remain the majority in every training fold. Cohort composition, graph geometry and the number of per-epoch graph updates also change. The hypothetical labels of SIK3/RPS6KB2 are assumptions used in this version, not new experimental truth. Their held-out scores test consistency with those assumptions; they do not validate the assumptions biologically.

Use the full ten-target metrics together with the original-eight subset, hypothetical-target table and rank-pair decomposition. A higher aggregate metric may be driven by the two assumed negatives. Success or failure of this specific frozen GNN cannot prove general geometric separability or impossibility for all MD contact-graph methods. No additional fitting, threshold change, score inversion or model selection occurs in this descriptive audit.
'''
 (ROOT/'INTERPRETATION.md').write_text(text);sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest();(ROOT/'validation/fold_prior_audit.json').write_text(json.dumps({'status':'DESCRIPTIVE_FOLD_COUNTS_VERIFIED','new_predictive_model_fitted':False,'rows':len(rows),'fold_prior_table_sha256':sha(out),'interpretation_sha256':sha(ROOT/'INTERPRETATION.md'),'script_sha256':sha(Path(__file__))},indent=2)+'\n')
 print('Descriptive audit saved for20folds; no new predictive model fitted.')
if __name__=='__main__':main()
