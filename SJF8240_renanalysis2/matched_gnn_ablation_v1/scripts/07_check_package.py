"""Integrity checks for allmatchedfolds,checkpoints,predictions andsource."""
from pathlib import Path
import sys,json,ast
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'));from contracts import *
def main():
 for rel,h in json.loads((ROOT/'config/software_manifest.json').read_text()).items():
  require(sha(ROOT/rel)==h,'Software/configurationchanged:'+rel)
  if rel.endswith('.py'):ast.parse((ROOT/rel).read_text())
 f=frozen();p=protocol();complete=json.loads((ROOT/'training_completion.json').read_text());require(complete['folds']==24,'Incomplete24folds');require(len({(r['variant'],r['pair']) for r in complete['records']})==24,'Duplicatefold')
 for r in complete['records']:
  path=ROOT/'results'/r['variant']/(r['pair']+'.pt');require(sha(path)==r['checkpoint_sha256'],'Weightdrift');require(r['source_manifest_sha256']==f,'Foldsource drift');require(json.loads(path.with_suffix('.replay.json').read_text())['status']=='PASS','Immediatefoldreplaymissing')
 v=json.loads((ROOT/'validation/independent_cpu_verification.json').read_text());require(v['status']=='PASS' and v['models']==24 and v['heldout_scores']==48 and v['same_initialization_all12_pairs'],'CPUverificationincomplete')
 for name,h in json.loads((ROOT/'validation/comparison_receipt.json').read_text())['result_sha256'].items():require(sha(ROOT/'results'/name)==h,'Prediction/evidenceoutputdrift:'+name)
 require(json.loads((ROOT/'validation/experimental_review_receipt.json').read_text())['labels_changed'] is False,'Labelsmodifiedduringreview');print('PASS:24foldweights,48CPUreplayedheldoutscores,12matchedinitializations,frozenprotocol andoutputs.')
if __name__=='__main__':main()
