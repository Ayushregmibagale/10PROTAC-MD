"""Read-only final contracts and artifact hash checks."""
from pathlib import Path
import sys,json,hashlib,ast
import pandas as pd
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 for path in ROOT.rglob('*.py'):ast.parse(path.read_text())
 if (ROOT/'software_manifest.json').exists():
  for rel,h in json.loads((ROOT/'software_manifest.json').read_text()).items():assert sha(ROOT/rel)==h,'Software drift '+rel
 completion=json.loads((ROOT/'COMPLETION.json').read_text());assert completion['status']=='COMPLETE' and completion['known_targets']==8 and completion['unknown_targets']==2 and completion['unknown_scores']==240 and not completion['new_training'] and not completion['unknown_labels_used'] and completion['EPHA2_label']==1
 controls=json.loads((ROOT/'validation/feature_matching.json').read_text());assert controls['status']=='PASS' and len(controls['eight_original_targets_reextracted'])==8
 for row in controls['unknown_graphs']:assert sha(ROOT/row['path'])==row['sha256'] and row['label'] is None
 for row in json.loads((ROOT/'extraction_completion.json').read_text()):
  for rel,h in row['outputs'].items():assert sha(ROOT/'extracted'/row['target']/row['pose']/rel)==h
 for row in json.loads((ROOT/'reach50_completion.json').read_text()):assert sha(ROOT/'reach50'/row['target']/row['pose']/'packet.npz')==row['packet_sha256']
 for name in ['unknown_CPU_replay.json','independent_geometry.json','residue_numbering.json']:assert json.loads((ROOT/'validation'/name).read_text())['status']=='PASS'
 scores=pd.read_csv(ROOT/'unknown_fold_predictions.tsv',sep='\t');assert len(scores)==240 and scores.experimental_label.isna().all();assert scores.groupby(['target','variant']).size().eq(60).all();assert scores.replay_error.max()<=1e-7
 targets=pd.read_csv(ROOT/'mechanistic/ten_target_summary.tsv',sep='\t');assert len(targets)==10 and targets.poses.sum()==174 and (targets.group=='unknown').sum()==2
 lys=pd.read_csv(ROOT/'mechanistic/per_lysine_equal_pose_summary.tsv',sep='\t');assert 'canonical_mapping_verified' in lys and 'delta_ensemble_point_reach50' in lys and 'delta_legacy_reach50' in lys
 assert len(pd.read_csv(ROOT/'mechanistic/per_pose_window.tsv',sep='\t'))==174*4
 if (ROOT/'SAVED_FILES.json').exists():
  for rel,h in json.loads((ROOT/'SAVED_FILES.json').read_text()).items():assert sha(ROOT/rel)==h,'Saved artifact drift '+rel
 print('PASS:174poses,10targets,4matchinggraphs,240unknownscores,8exactmatchingcontrols,12independentframechecks;no training or unknown labels')
if __name__=='__main__':main()
