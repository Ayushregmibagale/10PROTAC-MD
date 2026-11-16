"""Verify staged pilot hashes, bound-head sterics and declared PTMs; no submission."""
import argparse,json,hashlib,importlib.util
from pathlib import Path

def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',type=Path,required=True);a=p.parse_args();a.bundle=a.bundle.resolve()
 here=Path(__file__).parent;steric=load('prepared_sterics',here/'04_audit_prepared_coverage.py');mods=load('declared_modifications',here/'12_check_declared_protein_modifications.py')
 manifest=json.loads((a.bundle/'INPUT_MANIFEST.json').read_text());rows=[]
 for row in manifest['controls']:
  folder=a.bundle/row['directory'];folder.resolve().relative_to(a.bundle)
  for name,digest in row['input_file_sha256'].items():
   if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=digest:raise ValueError('Staged input changed: '+row['compound']+'/'+name)
  result={'compound':row['compound'],'target':row['target'],'hash_integrity':'PASS','poi_head_sterics':steric.sterics(folder/'poi_head.sdf',folder/'protein.pdb',chains={'D'}),'vhl_head_sterics':steric.sterics(folder/'vhl_head.sdf',folder/'vbc.pdb',chains={'A','B','C'}),'qualification_status':'REVIEW_REQUIRED_NOT_RELEASED'}
  if (folder/'protein_modification_spec.json').exists():
   spec=json.loads((folder/'protein_modification_spec.json').read_text());result['declared_protein_modifications']=mods.inspect(folder/'protein.pdb',spec['protein_chains'],spec['protein_modifications'])
  rows.append(result)
 report={'status':'SIX_INPUT_HASH_AND_STERIC_AUDIT_COMPLETE','production_ready':False,'controls':rows,'limitations':['Steric values are diagnostics, not binding or complete structure validation','No relaxation, docking, MD or remote submission','Full MS143 PTM preservation/selection remains an end-to-end qualification requirement'],'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (a.bundle/'INPUT_AUDIT.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({'status':report['status'],'controls':[{k:r[k] for k in ['compound','target','hash_integrity']} for r in rows]},indent=2))
if __name__=='__main__':main()
