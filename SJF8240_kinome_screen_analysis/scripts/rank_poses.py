"""Read-only ranking of retained clustered pRosettaC candidates; no degradation model."""
from pathlib import Path
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
import json,hashlib,argparse,concurrent.futures,time
from collections import defaultdict
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner
ROOT=Path(__file__).resolve().parents[1]
def require(ok,msg):
 if not ok:raise ValueError(msg)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic_json(path,value):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix('.pending.json');q.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(q,p)
def parse(path):
 raw=Path(path).read_bytes();groups=defaultdict(lambda:{'xyz':[],'resids':[],'names':[],'ca':[]});lig=[];elements=[];seg=0;seen=set();score=None;header=None
 for line in raw.decode('ascii').splitlines():
  if line.startswith('TER'):seg+=1;continue
  if line.startswith('label '):header=line.split();continue
  if line.startswith('pose ') and header:
   fields=line.split()
   if len(fields)==len(header):score=float(fields[header.index('total')])
  if not line.startswith(('ATOM  ','HETATM')):continue
  if line[16:17] not in [' ','A']:continue
  name=line[12:16].strip();element=line[76:78].strip().upper() or name.lstrip('0123456789')[0]
  if element in ['H','D']:continue
  key=(seg,line[21],line[22:27],name);require(key not in seen,'Duplicate heavy atom');seen.add(key)
  xyz=(float(line[30:38]),float(line[38:46]),float(line[46:54]));require(np.isfinite(xyz).all(),'Nonfinite atom')
  resid=int(line[22:26]);chain=line[21]
  if chain=='X':lig.append(xyz);elements.append(element);continue
  if line.startswith('HETATM'):continue
  g=groups[(seg,chain)];g['xyz'].append(xyz);g['resids'].append(resid);g['names'].append(name)
  if name=='CA':g['ca'].append((resid,line[17:20].strip(),xyz))
 out=[]
 for key,g in groups.items():
  if not g['ca']:continue
  g.update(xyz=np.asarray(g['xyz']),resids=np.asarray(g['resids']),names=np.asarray(g['names']),key=key);g['sequence']=''.join(seq1(x[1],custom_map={'HID':'H','HIE':'H','HIP':'H','CYX':'C','MSE':'M'},undef_code='X') for x in g['ca']);out.append(g)
 return {'groups':out,'ligand':np.asarray(lig),'ligand_elements':sorted(elements),'score':score,'sha256':hashlib.sha256(raw).hexdigest()}
def role_identity(group,uid,refs):
 sequence=refs[uid];good=[]
 for resid,name,_ in group['ca']:
  aa=seq1(name,custom_map={'HID':'H','HIE':'H','HIP':'H','CYX':'C','MSE':'M'},undef_code='X');good.append(0<resid<=len(sequence) and sequence[resid-1]==aa)
 return float(np.mean(good))
def fit(ref,mobile,points):
 ref=np.asarray(ref);mobile=np.asarray(mobile);a=ref-ref.mean(0);b=mobile-mobile.mean(0);u,s,vt=np.linalg.svd(a.T@b);r=vt.T@np.diag([1,1,np.linalg.det(vt.T@u.T)])@u.T
 return (np.asarray(points)-ref.mean(0))@r.T+mobile.mean(0),float(np.sqrt(((a@r.T-b)**2).sum(1).mean()))
def ca_xyz(group):return np.array([x[2] for x in group['ca']])
def select_by_sequence(groups,sequence):
 matches=[g for g in groups if g['sequence']==sequence];require(len(matches)==1,'Missing or ambiguous sequence group');return matches[0]
def distances(a,b):return cKDTree(b).query(a)[0]
def model_context(pool,cfg):
 directory=Path(pool['directory']);files={name:directory/name for name in ['Side0_0001_ABC.pdb','PT0_0001.pdb','PT1_0001.pdb','protein.pdb']}
 proteinframe=directory/'Side1_0001_A.pdb'
 if not proteinframe.exists():proteinframe=files['protein.pdb']
 files['target_frame']=proteinframe;parameters=directory/'Protac_params.txt';require(parameters.is_file() and 'SJF8240' in parameters.read_text(),'Compound input identity is not SJF8240');files['protac_parameters']=parameters;refs=json.loads(Path(cfg['role_reference']).read_text())['sequences'];side=parse(files['Side0_0001_ABC.pdb']);roles={}
 for uid,role in [('Q15370','elob'),('Q15369','eloc'),('P40337','vhl')]:
  matches=[g for g in side['groups'] if role_identity(g,uid,refs)>=.95];require(len(matches)==1,'E3 identity unverified '+role);roles[role]=matches[0]
 protein=parse(proteinframe);require(len(protein['groups'])==1,'Multiple target groups in reference');target=protein['groups'][0]
 # Verify the prepared model against its accession-specific local AF structure.
 af=Path(cfg['kinase_structure_root'])/(pool['uid']+'.pdb');require(af.exists(),'No accession-specific sequence reference');files['accession_sequence_reference']=af;afgroups=parse(af)['groups'];require(afgroups,'No AF reference sequence');fullseq=max(afgroups,key=lambda g:len(g['ca']))['sequence'];aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-2,open_gap_score=-8,extend_gap_score=-.5);identity=aligner.align(fullseq,target['sequence'])[0].counts().identities/len(target['sequence']);require(identity>=.95,'Prepared target does not match accession reference')
 full=parse(cfg['full_e3_reference']);froles={}
 for uid,role in [('Q15370','elob'),('Q15369','eloc'),('P40337','vhl')]:froles[role]=next(g for g in full['groups'] if role_identity(g,uid,refs)>=.95)
 anchor=ca_xyz(froles['eloc']);require(len(anchor)==len(roles['eloc']['ca']),'Anchor correspondence');vhlinframe,anchorerror=fit(anchor,ca_xyz(roles['eloc']),froles['vhl']['xyz']);require(anchorerror<=cfg['anchor_rmsd_A'],'Reference E3 fit failed')
 pt0=parse(files['PT0_0001.pdb']);pt1=parse(files['PT1_0001.pdb'])
 # Heads may use a non-X chain; they contain no CA-bearing protein groups.
 def head_xyz(path):
  xyz=[]
  for l in path.read_text().splitlines():
   if l.startswith(('ATOM  ','HETATM')) and not l[12:16].strip().startswith('H'):xyz.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
  require(xyz,'Empty reference head');return np.array(xyz)
 vnear=distances(vhlinframe,head_xyz(files['PT0_0001.pdb']));pnear=distances(target['xyz'],head_xyz(files['PT1_0001.pdb']));vsite=np.unique(froles['vhl']['resids'][vnear<cfg['pocket_definition_A']]);psite=np.unique(target['resids'][pnear<cfg['pocket_definition_A']]);require(len(vsite)>=2 and len(psite)>=2,'Insufficient head-reference pocket region')
 return {'roles':roles,'target':target,'full_vhl':froles['vhl'],'anchor_reference':anchor,'vhl_pocket_resids':vsite,'poi_pocket_resids':psite,'identity':float(identity),'source_sha256':{k:sha(p) for k,p in files.items()},'reference_head_contacts_A':{'VHL':float(vnear.min()),'POI':float(pnear.min())}}
def pose_row(path,context,cfg,pool):
 d=parse(path);lig=d['ligand'];require(lig.shape==(79,3),'Expected complete79-heavy-atom PROTAC');g={role:select_by_sequence(d['groups'],r['sequence']) for role,r in context['roles'].items() if role!='vhl'};refs=json.loads(Path(cfg['role_reference']).read_text())['sequences'];observed_vhl=[g for g in d['groups'] if len(g['ca'])>=20 and role_identity(g,'P40337',refs)>=.95];require(len(observed_vhl)==1,'Missing or ambiguous true-VHL fragment');g['vhl']=observed_vhl[0];poi=select_by_sequence(d['groups'],context['target']['sequence']);vhl=context['full_vhl'];vfull,rmsd=fit(context['anchor_reference'],ca_xyz(g['eloc']),vhl['xyz']);vca,_=fit(context['anchor_reference'],ca_xyz(g['eloc']),ca_xyz(vhl));vca_map={r:xyz for (r,_,_),xyz in zip(vhl['ca'],vca)};cerror=float(np.sqrt(np.mean([np.sum((vca_map[r]-xyz)**2) for r,_,xyz in g['vhl']['ca']])))
 dv=distances(lig,vfull);dk=distances(lig,poi['xyz']);observed=distances(lig,g['vhl']['xyz']);vp=vfull[np.isin(vhl['resids'],context['vhl_pocket_resids'])];pp=poi['xyz'][np.isin(poi['resids'],context['poi_pocket_resids'])];require(len(vp)>0 and len(pp)>0,'Pocket residue mapping');dvp=distances(lig,vp);dkp=distances(lig,pp);e3=np.vstack([vfull,g['elob']['xyz'],g['eloc']['xyz']]);all_lig=distances(lig,e3);inter=distances(e3,poi['xyz']);clash=cfg['severe_clash_A'];nv=int((dvp<cfg['contact_A']).sum());nk=int((dkp<cfg['contact_A']).sum());bridge=bool(dv.min()<cfg['contact_A'] and dk.min()<cfg['contact_A']);pockets=bool(dvp.min()<cfg['contact_A'] and dkp.min()<cfg['contact_A'] and min(nv,nk)>=cfg['min_ligand_contact_atoms']);fitok=bool(rmsd<=cfg['anchor_rmsd_A'] and cerror<=cfg['vhl_backbone_rmsd_A']);clashfree=bool(min(all_lig.min(),dk.min(),inter.min())>=clash)
 return {'uid':pool['uid'],'gene':pool['gene'],'pool':pool['pool'],'path':str(path),'pdb_sha256':d['sha256'],'source_bytes':path.stat().st_size,'source_mtime_ns':path.stat().st_mtime_ns,'cluster':path.parent.name,'pose':path.name,'source_label':'unknown' if pool['uid'] not in cfg['known_labels'] else cfg['known_labels'][pool['uid']],'target_CA':len(poi['ca']),'observed_VHL_CA':len(g['vhl']['ca']),'reference_VHL_CA':len(vhl['ca']),'VHL_reference_fit_RMSD_A':rmsd,'VHL_fragment_backbone_RMSD_A':cerror,'ligand_VHL_full_min_A':float(dv.min()),'ligand_VHL_observed_min_A':float(observed.min()),'ligand_POI_min_A':float(dk.min()),'ligand_VHL_pocket_min_A':float(dvp.min()),'ligand_POI_pocket_min_A':float(dkp.min()),'VHL_pocket_ligand_contact_atoms':nv,'POI_pocket_ligand_contact_atoms':nk,'balanced_contact_atoms':min(nv,nk),'geometric_mean_contact_atoms':float(np.sqrt(nv*nk)),'ligand_scaffold_min_A':float(all_lig.min()),'full_E3_POI_min_A':float(inter.min()),'bridge':bridge,'pocket_proximity':pockets,'full_reference_fit_ok':fitok,'no_severe_clashes':clashfree,'passes_full_screen':bool(bridge and pockets and fitok and clashfree),'observed_bridge':bool(observed.min()<cfg['contact_A'] and dk.min()<cfg['contact_A']),'ligand_element_signature':''.join(d['ligand_elements']),'rosetta_total_score':d['score']}
def work(pool):
 cfg=json.loads((ROOT/'config/protocol.json').read_text());out=ROOT/'results/pools'/pool['key'];out.mkdir(parents=True,exist_ok=True)
 if (out/'completion.json').exists():
  r=json.loads((out/'completion.json').read_text());require(r['script_sha256']==sha(__file__) and r['protocol_sha256']==sha(ROOT/'config/protocol.json'),'Cached source/config drift')
  if r['status']=='COMPLETE':
   require(sha(out/'poses.tsv')==r['poses_sha256'],'Cached table drift')
   for name,h in r['source_sha256'].items():require(sha(Path(name))==h,'Cached context source changed')
  return r
 start=time.monotonic();rows=[];failures=[]
 try:
  context=model_context(pool,cfg)
 except Exception as e:
  atomic_json(out/'completion.json',{'status':'CONTEXT_REJECTED','pool':pool,'error':str(e),'script_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'config/protocol.json')});return {'status':'CONTEXT_REJECTED','pool':pool,'error':str(e)}
 for p in pool['poses']:
  try:rows.append(pose_row(Path(p),context,cfg,pool))
  except Exception as e:failures.append({'path':p,'error':str(e)})
 pd.DataFrame(rows).to_csv(out/'poses.tsv',sep='\t',index=False);atomic_json(out/'failures.json',failures);r={'status':'COMPLETE','pool':{k:v for k,v in pool.items() if k!='poses'},'expected':len(pool['poses']),'parsed':len(rows),'rejected':len(failures),'bridges':sum(x['bridge'] for x in rows),'quality_pass':sum(x['passes_full_screen'] for x in rows),'poses_sha256':sha(out/'poses.tsv'),'script_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'config/protocol.json'),'source_sha256':context['source_sha256'],'target_sequence_identity':context['identity'],'reference_head_contacts_A':context['reference_head_contacts_A'],'seconds':time.monotonic()-start};atomic_json(out/'completion.json',r);return r

def main():
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--limit',type=int);a=p.parse_args();pools=json.loads((ROOT/'config/pools.json').read_text());pools=pools[:a.limit] if a.limit else pools;rows=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as ex:
  for r in ex.map(work,pools):rows.append(r);print(json.dumps({'complete':len(rows),'expected':len(pools),'key':r['pool']['key'],'status':r['status'],'quality_pass':r.get('quality_pass')}),flush=True);atomic_json(ROOT/'STATUS.json',{'phase':'SCREENING','pools_complete':len(rows),'pools_expected':len(pools),'active_job':True})
 atomic_json(ROOT/('benchmark.json' if a.limit else 'completion.json'),{'pools':len(rows),'records':rows,'script_sha256':sha(__file__),'protocol_sha256':sha(ROOT/'config/protocol.json')})
if __name__=='__main__':main()
