"""Paired static graph regeneration from all retained poses; no docking."""
from pathlib import Path
import argparse,json,re,hashlib,time,sys
import numpy as np
from scipy.spatial.distance import cdist
ROOT=Path(__file__).resolve().parents[1]
AUDIT=Path('/home/ayushregmibagale/10PROTAC-MD/prosettac_robust')
sys.path.insert(0,str(ROOT/'scripts/lib'))
from protac_qc import ValidationError
PATTERN=re.compile(rb'^(?:ATOM  |HETATM)(?:.{6} CA |.{15}X).*$',re.M)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
COORDINATE = re.compile(rb'[+-]?\d+\.\d{3}')
OVERFLOW_RECORD = re.compile(rb'^\s*([+-]?\d+\.\d{3})\s*([+-]?\d+\.\d{3})\s*([+-]?\d+\.\d{3})\s+\d+\.\d{2}\s+[+-]?\d+\.\d{2}\s+[A-Za-z]{1,2}\s*(?:[0-9][+-])?\s*$')
def coordinates(line):
 fields=[line[30:38].strip(),line[38:46].strip(),line[46:54].strip()]
 if all(COORDINATE.fullmatch(x) for x in fields):return [float(x) for x in fields],False
 # Historical writer overflowed an 8-column coordinate field. Exact three-decimal
 # tokens plus the complete occupancy/B-factor/element suffix disambiguate it.
 match=OVERFLOW_RECORD.fullmatch(line[30:])
 if match is None:raise ValidationError('Ambiguous/malformed PDB coordinate record: '+repr(line))
 return [float(x) for x in match.groups()],True

def parse(data,cfg,recoveries=None):
 groups={k:[] for k in ['elob','vhl','poi','ligand']}
 for line in PATTERN.findall(data):
  name=line[12:16].strip().decode();chain=chr(line[21]);res=int(line[22:26]);resname=line[17:20].strip().decode();icode=chr(line[26]);element=line[76:78].strip().decode() if len(line)>=78 else name[0]
  if chain=='X' and element not in ['H','D'] and not name.startswith('H'):role='ligand'
  elif name=='CA' and chain=='C':role='vhl'
  elif name=='CA' and chain=='A' and 1<=res<=103:role='elob'
  elif name=='CA' and chain==cfg['kin_chain'] and (chain=='D' or res>=cfg['kin_min']):role='poi'
  else:continue
  xyz,recovered=coordinates(line)
  if recovered and recoveries is not None:recoveries.append({'record':line.decode(),'xyz_A':xyz})
  groups[role].append(((chain,res,icode,name,element,resname),xyz))
 for key in groups:
  groups[key].sort(key=lambda x:x[0])
  if not groups[key]:raise ValidationError('Empty role '+key)
  keys=[x[0] for x in groups[key]]
  if len(set(keys))!=len(keys):raise ValidationError('Duplicate node identifiers '+key)
  if not np.isfinite(np.array([x[1] for x in groups[key]])).all():raise ValidationError('Nonfinite coordinates')
 if len(groups['ligand'])!=79:raise ValidationError('Historical ligand heavy-atom count differs from 79')
 return groups

def validate_roles(groups,poi):
 from Bio.SeqUtils import seq1
 from Bio.Align import PairwiseAligner
 reference=json.loads((ROOT/'config/role_reference.json').read_text())
 expected={'elob':'Q15370','vhl':'P40337','poi':reference['targets'][poi]}
 aligner=PairwiseAligner();aligner.mode='local';aligner.match_score=2;aligner.mismatch_score=-1;aligner.open_gap_score=-10;aligner.extend_gap_score=-.5
 proof={}
 for role,uid in expected.items():
  sequence=''.join(seq1(x[0][5],undef_code='X') for x in groups[role]);target=reference['sequences'][uid]
  aligned=aligner.align(target,sequence)[0]
  pairs=[(i,j) for (a,b),(c,d) in zip(*aligned.aligned) for i,j in zip(range(a,b),range(c,d))]
  identity=sum(target[i]==sequence[j] for i,j in pairs)/max(len(pairs),1);coverage=len(pairs)/len(sequence)
  if identity<.9 or coverage<.9:raise ValidationError('Sequence identity does not support assigned role '+role)
  proof[role]={'accession':uid,'identity':identity,'observed_coverage':coverage,'CA_nodes':len(sequence)}
 return proof

class Accumulator:
 def __init__(self,groups,role):
  self.role=role;self.keys=[(r,x[0]) for r in [role,'poi','ligand'] for x in groups[r]]
  self.sizes=[len(groups[r]) for r in [role,'poi','ligand']];n=len(self.keys)
  self.count=np.zeros((n,n),dtype=np.uint32);self.total=np.zeros((n,n));self.minimum=np.full((n,n),np.inf);self.n=0
  self.x=np.zeros((n,9),dtype=np.float32);start=0;self.roles=np.concatenate([np.full(k,i) for i,k in enumerate(self.sizes)])
  elements={'C':0,'N':1,'O':2,'F':3,'S':4}
  for ri,r in enumerate([role,'poi','ligand']):
   g=groups[r];nums=np.array([x[0][1] for x in g]) if r!='ligand' else np.arange(len(g));pos=(nums-nums.min())/max(int(nums.max()-nums.min()),1)
   for j,(key,_) in enumerate(g):self.x[start+j,ri]=1;self.x[start+j,3]=pos[j];self.x[start+j,4+elements[key[4]]]=1
   start+=len(g)
  self.cross=self.roles[:,None]!=self.roles[None,:]
 def add(self,groups):
  keys=[(r,x[0]) for r in [self.role,'poi','ligand'] for x in groups[r]]
  if keys!=self.keys:raise ValidationError('Node/residue correspondence changed across retained poses')
  xyz=np.array([x[1] for r in [self.role,'poi','ligand'] for x in groups[r]])
  d=cdist(xyz,xyz);contact=(d<8.)&self.cross
  self.count+=contact;self.total+=np.where(contact,d,0);self.minimum=np.minimum(self.minimum,np.where(contact,d,np.inf));self.n+=1
 def save(self,dest,poi,label):
  dest.mkdir(parents=True,exist_ok=False);freq=self.count/self.n;i,j=np.nonzero(freq>=.005)
  edge=np.column_stack([freq[i,j],self.total[i,j]/self.count[i,j]/10,self.minimum[i,j]/10]).astype(np.float32)
  np.savez_compressed(dest/'graph.npz',x=self.x,edge_index=np.array([i,j],dtype=np.int64),edge_attr=edge,label=np.array([label]),poi=np.array([poi]),node_identifiers=np.array([json.dumps(x) for x in self.keys]))
  return {'poses':self.n,'nodes_by_role':self.sizes,'nodes':len(self.keys),'directed_edges':len(i),'graph_sha256':sha(dest/'graph.npz')}
def build(cfg,limit=None):
 poi=cfg['poi'];out=ROOT/('analysis/static' if limit is None else 'analysis/static_diagnostic')/poi
 if limit is None and (out/'receipt.json').exists():
  receipt=json.loads((out/'receipt.json').read_text())
  if receipt.get('status')!='STATIC_PAIRED_COMPLETE' or receipt.get('poi')!=poi:raise ValidationError('Invalid completion receipt: '+poi)
  if sha(out/'pose_provenance.json')!=receipt['pose_provenance_sha256']:raise ValidationError('Pose provenance drift: '+poi)
  provenance=json.loads((out/'pose_provenance.json').read_text())
  if sha(ROOT/cfg['pose_list'])!=cfg['pose_list_sha256'] or [p['path'] for p in provenance]!=(ROOT/cfg['pose_list']).read_text().splitlines():raise ValidationError('Completed pose list drift: '+poi)
  for role in ['elob','vhl']:
   if receipt['outputs'][role]['poses']!=len(provenance) or sha(out/role/'graph.npz')!=receipt['outputs'][role]['graph_sha256']:raise ValidationError('Completed graph drift: '+poi)
  print(poi+' verified complete',flush=True);return
 if out.exists():raise ValidationError('Partial output requires explicit review: '+str(out))
 paths=[Path(x) for x in (ROOT/cfg['pose_list']).read_text().splitlines()]
 if sha(ROOT/cfg['pose_list'])!=cfg['pose_list_sha256']:raise ValidationError('Pose list drift')
 if limit is not None:paths=paths[:limit]
 started=time.monotonic();acc=None;provenance=[]
 from concurrent.futures import ThreadPoolExecutor
 def read(p):
  before=p.stat();data=p.read_bytes();after=p.stat()
  if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise ValidationError('Source changed during read')
  return p,data,after
 # Bounded prefetch avoids loading an entire retained docking ensemble into RAM.
 pool=ThreadPoolExecutor(max_workers=4)
 from collections import deque
 def buffered():
  iterator=iter(paths);pending=deque()
  for _ in range(8):
   p=next(iterator,None)
   if p is not None:pending.append(pool.submit(read,p))
  while pending:
   result=pending.popleft().result();p=next(iterator,None)
   if p is not None:pending.append(pool.submit(read,p))
   yield result
 for n,(p,data,after) in enumerate(buffered(),1):
  recoveries=[]
  try:g=parse(data,cfg,recoveries)
  except (ValueError,IndexError) as exc:raise ValidationError(str(p)+': '+str(exc)) from exc
  if acc is None:
   proof=validate_roles(g,cfg['poi'])
   acc={role:Accumulator(g,role) for role in ['elob','vhl']}
  for x in acc.values():x.add(g)
  provenance.append({'path':str(p),'bytes':len(data),'mtime_ns':after.st_mtime_ns,'sha256':hashlib.sha256(data).hexdigest(),'coordinate_overflow_records':recoveries})
  if n%1000==0:print(json.dumps({'poi':poi,'processed':n,'total':len(paths),'elapsed_s':round(time.monotonic()-started,1)}),flush=True)
 pool.shutdown()
 out.mkdir(parents=True,exist_ok=False)
 label=cfg['historical_label'] if cfg['historical_label'] is not None else -1
 outputs={role:x.save(out/role,poi,label) for role,x in acc.items()}
 (out/'pose_provenance.json').write_text(json.dumps(provenance,separators=(',',':'))+'\n')
 receipt={'status':'STATIC_PAIRED_COMPLETE' if limit is None else 'DIAGNOSTIC_ONLY','poi':poi,'outputs':outputs,'sequence_role_proof':proof,'role_control':'elob and vhl use same pipeline, different E3 role; not exact historical replay','pose_provenance_sha256':sha(out/'pose_provenance.json'),'source_sha256':sha(__file__),'elapsed_s':time.monotonic()-started,'new_docking':False,'coordinate_overflow_records':sum(len(p['coordinate_overflow_records']) for p in provenance)}
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--poi',nargs='+');p.add_argument('--limit',type=int);p.add_argument('--workers',type=int,default=1);a=p.parse_args();records=json.loads((ROOT/'config/source_inventory.json').read_text())['records']
 records=[cfg for cfg in records if not a.poi or cfg['poi'] in a.poi]
 if a.workers>1:
  if a.limit is not None:raise ValidationError('Diagnostic parallel selection requires explicit separate outputs')
  from concurrent.futures import ProcessPoolExecutor
  with ProcessPoolExecutor(max_workers=a.workers) as pool:list(pool.map(build,records))
 else:
  for cfg in records:build(cfg,a.limit)
