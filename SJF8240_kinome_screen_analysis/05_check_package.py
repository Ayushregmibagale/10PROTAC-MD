from pathlib import Path
import hashlib,json,ast
ROOT=Path(__file__).resolve().parent
for name,expected in json.loads((ROOT/'config/package_manifest.json').read_text()).items():
 p=ROOT/name
 assert hashlib.sha256(p.read_bytes()).hexdigest()==expected, 'Changed package file: '+name
 if p.suffix=='.py':ast.parse(p.read_text(),filename=str(p))
print('PASS: packaged source hashes and syntax. Numerical analysis completion is separate.')
