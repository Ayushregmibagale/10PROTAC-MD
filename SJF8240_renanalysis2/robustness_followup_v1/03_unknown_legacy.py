"""Replay prediction-only targets using historical compatible MD-only checkpoints."""
from pathlib import Path
import sys,json,hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent;OLD=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis/corrected_md_only_v1');sys.path.insert(0,str(OLD/'lib'))
from common import *
ROOT=Path(__file__).resolve().parent
import pandas as pd

def main():
 torch.set_num_threads(1);binding=verify_frozen();rows=[]
 for role in ['vhl','elob']:
  p,g=load_inputs(role,'cpu');require(p['prediction_only']==['RPS6KB2_full','SIK3_full'] or set(p['prediction_only'])=={'RPS6KB2_full','SIK3_full'},'Unknown scope changed')
  for held in p['cohort']:
   path=OLD/'analysis'/role/('fold_'+held+'.pt');receipt=json.loads(path.with_suffix('.json').read_text());require(sha(path)==receipt['checkpoint_sha256'],'Checkpoint drift');state=torch.load(path,map_location='cpu',weights_only=True);require(not set(state['training_targets'])&set(p['prediction_only']),'Unknown used in training');m=MDEncoder();m.load_state_dict(state['model_state_dict'])
   for target in p['prediction_only']:
    score=predict(m,g[target]);saved=next(x for x in receipt['prediction_only'] if x['poi']==target)['score'];require(abs(score-saved)<=1e-7,'Unknown replay mismatch');rows.append({'model':'legacy_corrected_md_only_v1','role':role,'held_out_labelled_target':held,'unknown_target':target,'score':score,'predicted_degrader_at05':int(score>=.5),'experimental_label':None,'checkpoint_sha256':sha(path),'replay_absolute_error':abs(score-saved)})
 df=pd.DataFrame(rows);df.to_csv(ROOT/'unknown_legacy_fold_predictions.tsv',sep='\t',index=False);summary=df.groupby(['model','role','unknown_target']).agg(mean_score=('score','mean'),min_score=('score','min'),max_score=('score','max'),degrader_vote_fraction=('predicted_degrader_at05','mean')).reset_index();summary.to_csv(ROOT/'unknown_legacy_summary.tsv',sep='\t',index=False)
 (ROOT/'UNKNOWN_TARGETS.md').write_text('# SIK3 and RPS6KB2 prediction-only assessment\n\n'+summary.to_markdown(index=False)+'\n\nThese are predictions from the older corrected MD-only v1 model, independently replayed on CPU. Neither unknown was labelled or used in fitting. This is not experimental degradation evidence. ELOB is an alternative historical arm, not true VHL.\n\nThe AUROC 0.8333 temporal catalytic model cannot yet score these targets: their matching four-window 29-feature graphs and five catalytic-geometry annotations have not been built. Existing aggregated 3-feature graphs are not compatible. No missing catalytic features were silently zero-filled. Newest-model unknown classification remains unavailable until matching inputs are prepared and verified from existing trajectories.\n');print(summary.to_string(index=False))
if __name__=='__main__':main()
