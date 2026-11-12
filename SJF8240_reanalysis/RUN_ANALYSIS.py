"""Public entry point: validated, compound-specific inputs; never runs archived experiments."""
from pathlib import Path
import argparse,subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parent
if (ROOT/'robust_analysis_v1').is_dir():ROOT=ROOT.parent
PACKAGE=ROOT/'folder_audit_20261112/robust_analysis_v1'
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['validate','batch','extract','train','replay']);a,rest=p.parse_known_args()
 registry=json.loads((PACKAGE.parent/'reviewed_entrypoints.json').read_text())
 for name,digest in registry['package_source_sha256'].items():
  if hashlib.sha256((PACKAGE/name).read_bytes()).hexdigest()!=digest:raise SystemExit('Reviewed source changed: '+name)
 script={'validate':'validate_study.py','batch':'validate_batch.py','extract':'build_dynamic.py','train':'run_study.py','replay':'run_study.py'}[a.action]
 if a.action=='replay':
  if '--device' in rest and rest[rest.index('--device')+1]!='cpu':raise SystemExit('Replay requires an explicit CPU environment')
  rest+=['--verify','--device','cpu']
 raise SystemExit(subprocess.call([sys.executable,str(PACKAGE/script),*rest]))
