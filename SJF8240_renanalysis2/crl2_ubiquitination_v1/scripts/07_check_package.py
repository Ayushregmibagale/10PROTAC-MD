"""Check immutable code/configuration,completed cohort and numerical receipts."""
from pathlib import Path
import json,hashlib,ast
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 for rel,h in json.loads((ROOT/'config/software_manifest.json').read_text()).items():
  if sha(ROOT/rel)!=h:raise SystemExit('Source/configuration changed: '+rel)
  if rel.endswith('.py'):ast.parse((ROOT/rel).read_text())
 for rel,h in json.loads((ROOT/'config/execution_manifest.json').read_text()).items():
  if sha(ROOT/rel)!=h:raise SystemExit('Extraction manifest drift: '+rel)
 c=json.loads((ROOT/'completion.json').read_text())
 if c['cases']!=141 or c['frames']!=28341:raise SystemExit('Incomplete cohort')
 for r in c['records']:
  p=ROOT/'results/poses'/r['poi']/r['pose']/'packet.npz'
  if sha(p)!=r['packet_sha256']:raise SystemExit('Packet changed: '+str(p))
 a=json.loads((ROOT/'validation/analysis_receipt.json').read_text())
 for name,h in a['outputs_sha256'].items():
  if sha(ROOT/'results'/name)!=h:raise SystemExit('Analysis output drift: '+name)
 if json.loads((ROOT/'validation/independent_verification.json').read_text())['status']!='PASS':raise SystemExit('Independent verification incomplete')
 print('PASS: code,configuration,141 packets,28,341 frames and numerical result receipts. Scientific validation remains separate.')
if __name__=='__main__':main()
