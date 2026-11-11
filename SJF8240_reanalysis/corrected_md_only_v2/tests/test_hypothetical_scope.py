from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from common import config,ROOT
import numpy as np
def test_hypothesis_does_not_overwrite_unknown_source_labels_or_features():
 p=config();rows=json.loads((ROOT/'config/input_manifest.json').read_text())['entries']
 for row in rows:
  if row['poi'] not in p['hypothetical_targets']:continue
  assert row['label_hypothetical'] and row['experimental_label'] is None
  with np.load(row['source_path'],allow_pickle=False) as old,np.load(ROOT/row['path'],allow_pickle=False) as new:
   assert old['label'][0]==-1 and new['label'][0]==0
   for key in old.files:
    if key!='label':assert old[key].dtype==new[key].dtype and np.array_equal(old[key],new[key])
def test_all_original_v1_hyperparameters_and_unweighted_loss_are_preserved():
 p=config();old=json.loads((ROOT/'provenance/v1_protocol.json').read_text())
 for key in ['epochs','seed','hidden','heads','dropout','learning_rate','weight_decay','optimizer','scheduler','loss','decision_threshold']:
  assert p[key]==old[key]
 assert p['loss']=='unweighted BCEWithLogits' and not p['E2_scalar_branch']
