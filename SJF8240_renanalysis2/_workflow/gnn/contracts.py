"""Fail-closed contracts for one compound; no SJF8240 cohort defaults."""
from pathlib import Path
import hashlib,json,os,tempfile,math
import numpy as np
class ContractError(ValueError):pass
def require(ok,message):
 if not ok:raise ContractError(message)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def atomic_json(path,value):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);name=None
 try:
  with tempfile.NamedTemporaryFile(mode='w',dir=path.parent,delete=False) as f:
   name=f.name;json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
  os.replace(name,path);name=None
 finally:
  if name is not None and Path(name).exists():Path(name).unlink()
def bound_file(base,row):
 p=Path(row['path']);p=p if p.is_absolute() else base/p;require(p.is_file(),'Missing input '+str(p));require(sha(p)==row['sha256'],'Input hash mismatch '+str(p));return p.resolve()
def graph_arrays(path,width,windows=1):
 with np.load(path,allow_pickle=False) as z:
  result=[]
  for w in range(windows):
   suffix='' if windows==1 else str(w);x=z['x'+suffix];ei=z['edge_index'+suffix];ea=z['edge_attr'+suffix]
   require(x.ndim==2 and len(x)>0 and x.shape[1]==width and x.dtype.kind in 'fi' and np.isfinite(x).all(),'Invalid node tensor')
   require(ei.ndim==2 and ei.shape[0]==2 and ei.dtype.kind in 'iu','Invalid edge indices')
   require(ea.shape==(ei.shape[1],3) and np.isfinite(ea).all(),'Invalid edge features')
   require(not ei.size or (ei.min()>=0 and ei.max()<len(x)),'Edge index out of bounds')
   require(np.all((ea[:,0]>=0)&(ea[:,0]<=1)) and np.all(ea[:,1:]>=0),'Invalid contact frequency/distance')
   require(np.all(x[:,:2].sum(1)==1) and np.all((x[:,:2]==0)|(x[:,:2]==1)),'Invalid node protein roles')
   require(np.all((x[:,2]>=0)&(x[:,2]<=1)),'Invalid residue position')
   if width==29:
    require(np.all(x[:,3:23].sum(1)==1) and np.all((x[:,3:23]==0)|(x[:,3:23]==1)),'Invalid amino-acid one-hot')
    require(np.all((x[:,23]==0)|(x[:,23]==1)) and np.all(x[:,24:]>=0),'Invalid lysine features')
    require(np.all(x[x[:,23]==0,24:]==0),'Lysine descriptors assigned to non-POI-lysine nodes')
    require(np.all(x[:,27:29]<=1),'Invalid lysine occupancies')
    require(np.all(x[:,23]<=x[:,1]) and np.all(x[:,23]<=x[:,11]),'Lysine mask must identify POI lysine residues')
   result.append((x.copy(),ei.copy(),ea.copy()))
  if windows>1:
   require(all(np.array_equal(result[0][0][:,:24],a[0][:,:24]) for a in result),'Window node identities/order changed')
  return result

def validate_study(path):
 path=Path(path).resolve();base=path.parent;c=json.loads(path.read_text());require(c.get('schema_version')==1,'Study schema version required')
 require(isinstance(c.get('compound_id'),str) and c['compound_id'].strip(),'Compound identity required');require(c.get('units')=={'coordinates':'angstrom','time':'ns'},'Explicit angstrom/ns units required')
 require(c['model'] in ['original_md','lysine_temporal'],'Unsupported model');require(c['evaluation'] in ['leave_pair_out','lopo'],'Unsupported evaluation')
 t=c['training'];require(type(t['epochs']) is int and t['epochs']>0 and type(t['seed']) is int,'Invalid budget/seed');require(t['loss'] in ['unweighted','balanced'],'Explicit supported loss required');require(t['threshold']==.5,'Threshold must be frozen at 0.5')
 for k in ['lr','weight_decay']:require(type(t[k]) in [int,float] and math.isfinite(t[k]) and t[k]>=0,'Invalid optimizer setting')
 require(t['lr']>0,'Positive learning rate required');roles=c['roles'];require(roles and set(roles)<= {'vhl','elob'},'Only reviewed VHL/ELOB graph arms supported')
 for role,entry in roles.items():
  require(entry['uniprot']=={'vhl':'P40337','elob':'Q15370'}[role],'Role identity mismatch: '+role)
  bound_file(base,entry['identity_evidence'])
 rows=c['targets'];ids=[r['target_id'] for r in rows];require(len(ids)==len(set(ids)) and all(ids),'Unique target IDs required');groups=[r['protein_group'] for r in rows];require(len(groups)==len(set(groups)) and all(groups),'Duplicate biological target/alias group')
 labels={};paths={};width,windows=(3,1) if c['model']=='original_md' else (29,4)
 for r in rows:
  target=r['target_id'];label=r['label'];status=r['label_status'];require(set(r['graphs'])==set(roles),'Missing or extra role graph '+target)
  require(status in ['confirmed_experimental','confirmed_user','unknown'],'Hypothetical labels cannot be training evidence')
  if status=='unknown':require(label is None,'Unknown target must have null label')
  else:
   require(type(label) is int and label in [0,1],'Invalid binary label');evidence=bound_file(base,r['label_evidence']);values=json.loads(evidence.read_text());require(values['labels'].get(target)==label,'Label evidence disagrees with study');labels[target]=label
  for role,g in r['graphs'].items():
   require(g.get('compound_id')==c['compound_id'] and g.get('target_id')==target and g.get('role')==role,'Graph scope mismatch')
   p=bound_file(base,g);graph_arrays(p,width,windows)
   with np.load(p,allow_pickle=False) as data:
    for field,expected in [('poi',target),('compound_id',c['compound_id']),('role',role)]:
     if field in data:require(data[field].shape==(1,) and str(data[field][0])==expected,'Embedded graph identity mismatch')
    if 'label' in data:require(data['label'].shape==(1,) and int(data['label'][0])==(label if label is not None else -1),'Graph label differs from study')
   paths[(target,role)]=p
 require(set(labels.values())=={0,1},'Both labelled classes required')
 require(min(sum(v==y for v in labels.values()) for y in [0,1])>=2,'At least two targets per class required for grouped CV')
 return c,labels,paths

def folds(labels,evaluation):
 targets=sorted(labels)
 if evaluation=='leave_pair_out':held=[(p,n) for p in targets if labels[p]==1 for n in targets if labels[n]==0]
 else:held=[(x,) for x in targets]
 rows=[]
 for group in held:
  training=[x for x in targets if x not in group];require(set(labels[x] for x in training)=={0,1},'Single-class training fold');key='__'.join(group);seed_key=json.dumps(list(group),separators=(',',':'));rows.append({'id':hashlib.sha256(seed_key.encode()).hexdigest()[:16],'name':key,'held_out':list(group),'training':training})
 return rows

def balanced_weights(labels,training):
 n=len(training);counts={v:sum(labels[x]==v for x in training) for v in [0,1]};require(min(counts.values())>0,'Both classes required');return {x:n/(2*counts[labels[x]]) for x in training}

def rank_metrics(records):
 require(records,'No predictions');sensitivity=float(np.mean([r['positive_score']>=.5 for r in records]));specificity=float(np.mean([r['negative_score']<.5 for r in records]));rank=float(np.mean([float(r['positive_score']>r['negative_score'])+.5*float(r['positive_score']==r['negative_score']) for r in records]));return {'matched_pair_AUROC':rank,'sensitivity':sensitivity,'specificity':specificity,'balanced_accuracy':(sensitivity+specificity)/2,'pairs':len(records)}

def verify_sequence_identity(observed,expected_uniprot,references,min_coverage=.95):
 from Bio.Align import PairwiseAligner
 require(observed and expected_uniprot in references,'Missing reference sequence')
 aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-2,open_gap_score=-8,extend_gap_score=-.5);scores={}
 for uid,reference in references.items():
  alignment=aligner.align(reference,observed)[0];scores[uid]=alignment.counts().identities/len(observed)
 best=max(scores.values());require(scores[expected_uniprot]>=min_coverage,'Observed sequence does not match declared identity');require([uid for uid,score in scores.items() if abs(score-best)<1e-12]==[expected_uniprot],'Wrong or ambiguous protein identity');return scores[expected_uniprot]

def joint_contact_observation(left,right,cutoff):
 left=np.asarray(left,float);right=np.asarray(right,float);require(left.ndim==right.ndim==1 and np.isfinite(left).all() and np.isfinite(right).all(),'Invalid per-frame contact distances');return np.outer(left<cutoff,right<cutoff)
