"""Verify four pinned public PDB inputs; --download fetches missing files only."""
from pathlib import Path
import argparse,json,hashlib,urllib.request,os
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--download',action='store_true');a=p.parse_args();ledger=json.loads((ROOT/'config/reference_sources.json').read_text())
 for code,row in ledger.items():
  dest=ROOT/'references'/(code+'.cif')
  if not dest.exists():
   if not a.download:raise SystemExit('Missing reference '+code+'; use --download')
   with urllib.request.urlopen(row['url'],timeout=45) as response:data=response.read()
   if hashlib.sha256(data).hexdigest()!=row['sha256']:raise SystemExit('Remote reference changed; explicit review/new version required: '+code)
   tmp=dest.with_suffix('.pending');tmp.write_bytes(data);os.replace(tmp,dest)
  if hashlib.sha256(dest.read_bytes()).hexdigest()!=row['sha256']:raise SystemExit('Reference hash mismatch: '+code)
 print('PASS: all four pinned PDB references.')
if __name__=='__main__':main()
