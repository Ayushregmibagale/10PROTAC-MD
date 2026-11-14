"""Verify the four pinned public reference inputs against the native-audit ledger."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];PRIOR=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2/crl2_ubiquitination_v1')
def main():
 audit=json.loads((ROOT/'validation/native_reference_audit.json').read_text())['sources_sha256']
 for original,digest in audit.items():
  name=Path(original).name;path=(ROOT/'references'/name) if name in ['8R5H.cif','SY8.cif'] else PRIOR/'references'/name
  if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:raise SystemExit('Reference missing/changed:'+str(path))
 print('PASS:four public reference inputs bound to the native audit.')
if __name__=='__main__':main()
