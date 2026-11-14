"""Check saved reproducibility and software results,not scientific deployment."""
from pathlib import Path
import json,hashlib,ast
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 for rel,h in json.loads((ROOT/'config/software_manifest.json').read_text()).items():
  if sha(ROOT/rel)!=h:raise SystemExit('Source/configuration changed:'+rel)
  if rel.endswith('.py'):ast.parse((ROOT/rel).read_text())
 for folder in [ROOT,ROOT/'chemical_tail_sensitivity']:
  for rel,h in json.loads((folder/'config/execution_manifest.json').read_text()).items():
   if sha(folder/rel)!=h:raise SystemExit('Frozen extraction drift:'+str(folder/rel))
  c=json.loads((folder/'completion.json').read_text())
  if c['cases']!=141 or c['frames']!=28341:raise SystemExit('Incomplete MDcohort')
  for r in c['records']:
   if sha(folder/'results/poses'/r['poi']/r['pose']/'packet.npz')!=r['packet_sha256']:raise SystemExit('Packet drift')
  if json.loads((folder/'validation/independent_verification.json').read_text())['status']!='PASS':raise SystemExit('Independent MDreplay missing')
  for name,h in json.loads((folder/'validation/comparison_receipt.json').read_text())['outputs_sha256'].items():
   if sha(folder/'results'/name)!=h:raise SystemExit('Result drift:'+name)
 g=json.loads((ROOT/'validation/gnn_input_manifest.json').read_text())
 for r in g['graphs']:
  if sha(ROOT/r['path'])!=r['sha256']:raise SystemExit('GNNinput drift')
 if json.loads((ROOT/'validation/gnn_runtime_check.json').read_text())['status']!='PASS':raise SystemExit('GNNruntime check missing')
 if json.loads((ROOT/'validation/kinome_receipt.json').read_text())['poses']!=723:raise SystemExit('Incomplete kinomeannotation')
 if json.loads((ROOT/'validation/final_tail_replay.json').read_text())['status']!='PASS':raise SystemExit('Tail geometry replay absent')
 print('PASS:two141-pose/28,341-frame analyses,16 GNNinputs,723 kinomeannotations,source/resultchecks. No biologicalclassifier release.')
if __name__=='__main__':main()
