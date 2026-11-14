from pathlib import Path
import json,hashlib,os
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def atomic_json(p,value):
 p=Path(p);p.parent.mkdir(exist_ok=True,parents=True);tmp=p.with_name(p.name+'.pending');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def frozen():
 manifest=json.loads((ROOT/'config/execution_manifest.json').read_text())
 for rel,h in manifest.items():require(sha(ROOT/rel)==h,'Frozen code/config/input drift:'+rel)
 return sha(ROOT/'config/execution_manifest.json')
def protocol():return json.loads((ROOT/'config/protocol.json').read_text())
