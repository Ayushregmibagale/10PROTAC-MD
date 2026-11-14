"""Independent NumPy distance replay for previous top10 and receipt checks."""
from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent));import check_engagement as c

def raw_segments(path):
 segments={};index=0;lig={}
 for line in Path(path).read_text().splitlines():
  if line.startswith('TER'):index+=1
  if not line.startswith(('ATOM  ','HETATM')):continue
  name=line[12:16].strip();element=line[76:78].strip().upper() or name[0]
  if element in ['H','D']:continue
  xyz=[float(line[30:38]),float(line[38:46]),float(line[46:54])]
  if line[21]=='X':lig[name]=(element,xyz)
  elif line.startswith('ATOM  '):segments.setdefault((index,line[21]),[]).append((name,element,xyz))
 return segments,lig

def main():
 receipt=json.loads((HERE/'completion.json').read_text())
 for name,digest in receipt['outputs_sha256'].items():assert c.sha(HERE/name)==digest,name
 d=pd.read_csv(HERE/'previous_top10_comparison.tsv',sep='\t');checks=[]
 for row in d.itertuples():
  assert c.sha(row.path)==row.pdb_sha256
  segments,lig=raw_segments(row.path)
  # Kinase is the larger CA-bearing TER segment within chain A, not ELOB.
  target=max((g for (segment,chain),g in segments.items() if chain=='A'),key=lambda g:sum(a[0]=='CA' for a in g))
  xyz=np.array([v[1] for v in lig.values()]);p=np.array([a[2] for a in target]);ds=np.sqrt(np.sum((xyz[:,None,:]-p[None,:,:])**2,axis=2));vr=np.array([c.RADII[v[0]] for v in lig.values()])[:,None]+np.array([c.RADII[a[1]] for a in target])[None,:]-ds
  assert abs(vr.max()-row.ligand_POI_max_overlap_A)<1e-10
  internal=np.linalg.norm(np.array(lig['O11'][1])-np.array(lig['C16'][1]));assert 1.52+1.7-internal>.8
  checks.append({'gene':row.gene,'ligand_POI_overlap_A':float(vr.max()),'O11_C16_distance_A':float(internal),'source_sha256':row.pdb_sha256})
 # Numeric threshold and atom-specific radius sanity checks, distinct from ranking.
 assert abs(c.overlap(np.array([[0.,0,0]]),['C'],np.array([[3.,0,0]]),['C'])['max_overlap_A']-.4)<1e-12
 assert c.overlap(np.array([[0.,0,0]]),['O'],np.array([[3.,0,0]]),['O'])['max_overlap_A']<.1
 mapping=json.loads((HERE/'head_mapping_receipts.json').read_text());assert all(not r['heads']['VHL']['reference_engagement_verifiable'] and not r['heads']['POI']['reference_engagement_verifiable'] for r in mapping)
 result={'status':'PASS','method':'Independent raw PDB segment parser and direct NumPy broadcasting for previous top10 kinase contacts; explicit linker O11-C16 distance replay; output hashes; overlap formula checks. Not an independent full rescreen or biological validation.','poses_replayed':checks,'all_reference_pools_untrusted':len(mapping),'verification_script_sha256':c.sha(__file__)}
 (HERE/'independent_verification.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS: 10 raw-source distance replays, linker overlap, all output hashes and threshold arithmetic.')
if __name__=='__main__':main()
