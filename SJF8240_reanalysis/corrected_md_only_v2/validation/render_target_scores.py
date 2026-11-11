"""Render the ten held-out target score pairs without model fitting."""
from pathlib import Path
import json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=json.loads((ROOT/'config/protocol.json').read_text());frame=pd.read_csv(ROOT/'results/paired_predictions.tsv',sep='\t')
 assert len(frame)==10 and set(frame.poi)==set(p['cohort'])
 text='# All ten held-out target scores\n\nScores are uncalibrated model outputs for class1 (degraded), with fixed threshold0.5. Every target is excluded from its own nine-target fitting cohort. SIK3/RPS6KB2 labels below are hypotheses, not known experimental outcomes.\n\n| Target | Evaluation label | Label status | ELOB score | ELOB class | VHL score | VHL class |\n|---|---:|---|---:|---:|---:|---:|\n'
 for _,row in frame.sort_values('poi').iterrows():
  status='hypothetical negative; experimental label unknown' if row.poi in p['hypothetical_targets'] else 'user-confirmed label'
  for role in ['elob','vhl']:
   receipt=json.loads((ROOT/'analysis'/role/f'fold_{row.poi}.json').read_text());assert row.poi not in receipt['training_targets'] and len(receipt['training_targets'])==9
  text+=f"| {row.poi} | {int(row.label)} | {status} | {row.score_elob:.8f} | {int(row.pred_elob)} | {row.score_vhl:.8f} | {int(row.pred_vhl)} |\n"
 path=ROOT/'results/TARGET_SCORES.md';path.write_text(text);sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest();(ROOT/'validation/target_score_table.json').write_text(json.dumps({'status':'ALL_TEN_TARGET_SCORES_AND_HOLDOUTS_CHECKED','paired_prediction_sha256':sha(ROOT/'results/paired_predictions.tsv'),'rendered_table_sha256':sha(path),'renderer_sha256':sha(Path(__file__)),'new_fitting':False},indent=2)+'\n');print(text)
if __name__=='__main__':main()
