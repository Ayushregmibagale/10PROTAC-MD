"""Read-only source/registry/frozen-manifest check; never executes discovered scripts."""
from pathlib import Path
import argparse,json,hashlib,ast

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(root,output):
 root=Path(root).resolve();registry=json.loads((root/'folder_audit_20261112/reviewed_entrypoints.json').read_text());rows=[];errors=[]
 for p in sorted(root.rglob('*.py')):
  rel=str(p.relative_to(root))
  if any(x in p.parts for x in ['.runtime','runtime','__pycache__','.pytest_cache']):continue
  try:ast.parse(p.read_text());syntax='PASS'
  except Exception as exc:syntax=str(exc);errors.append({'path':rel,'error':'syntax'})
  expected=registry['project_source_sha256'].get(rel);status='MATCH' if expected==sha(p) else 'UNREVIEWED_OR_CHANGED';rows.append({'path':rel,'sha256':sha(p),'syntax':syntax,'review':status})
  if status!='MATCH':errors.append({'path':rel,'error':status})
 for version in ['corrected_md_only_v1','corrected_md_only_v2','md_weighted_v1','md_weighted_e2_v1','mechanistic_v3']:
  for name in (['pair_source_hashes.json','temporal_source_hashes.json'] if version=='mechanistic_v3' else ['execution_manifest.json']):
   manifest=root/version/'config'/name
   for rel,h in json.loads(manifest.read_text()).items():
    if not (root/version/rel).is_file() or sha(root/version/rel)!=h:errors.append({'path':version+'/'+rel,'error':'frozen source changed'})
 report={'status':'PASS' if not errors else 'FAIL','scope':'Project-owned Python, excluding installed runtime dependencies; no discovered script is executed','files':len(rows),'errors':errors,'records':rows};Path(output).write_text(json.dumps(report,indent=2)+'\n');print(report['status'],len(rows),'project files;',len(errors),'errors');return bool(errors)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();raise SystemExit(main(a.root,a.output))
