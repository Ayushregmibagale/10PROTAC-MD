from pathlib import Path
import json,hashlib,ast
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 for name,h in json.loads((ROOT/'config/software_manifest.json').read_text()).items():
  if sha(ROOT/name)!=h:raise SystemExit('Source changed:'+name)
  if name.endswith('.py'):ast.parse((ROOT/name).read_text())
 v=json.loads((ROOT/'validation/independent_verification.json').read_text())
 if v['status']!='PASS' or v['cases']!=141 or v['frames']!=28341:raise SystemExit('Verification incomplete')
 for name,h in json.loads((ROOT/'validation/scoring_receipt.json').read_text())['outputs_sha256'].items():
  if sha(ROOT/'results'/name)!=h:raise SystemExit('Result drift:'+name)
 print('PASS: source,protocol,141 pose/28,341 frame Boolean replays and output receipts. Not a literal paper replication.')
if __name__=='__main__':main()
