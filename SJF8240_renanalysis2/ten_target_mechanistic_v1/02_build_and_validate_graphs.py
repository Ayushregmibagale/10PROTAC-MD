"""Validate matching extraction against archived controls; build unknown-only inputs."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT))
import importlib.util
s=importlib.util.spec_from_file_location('extraction',ROOT/'01_extract_matching_inputs.py');e=importlib.util.module_from_spec(s);s.loader.exec_module(e)
AA='ACDEFGHIKLMNPQRSTVWY';ALIASES={'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H','CYX':'C','MSE':'M'}
def base_graph(packets):
 first=packets[0]
 for z in packets:
  for k in ['lys_resids','lys_node_indices','resids_A','resids_C','resids_D','resnames_A','resnames_C','resnames_D']:e.require(np.array_equal(z[k],first[k]),'Pose correspondence changed '+k)
 nl=len(first['resids_C']);nr=len(first['resids_D']);base=np.zeros((nl+nr,29),np.float32);lyscount=len(first['lys_resids'])
 for offset,c,ri in [(0,'C',0),(nl,'D',1)]:
  ids=first['resids_'+c];base[offset:offset+len(ids),ri]=1;base[offset:offset+len(ids),2]=(ids-ids.min())/max(int(ids.max()-ids.min()),1)
  for j,name in enumerate(first['resnames_'+c]):
   a=e.seq1(str(name),custom_map=ALIASES,undef_code='X');e.require(a in AA,'Unknown amino acid');base[offset+j,3+AA.index(a)]=1
 arr={}
 for w in range(4):
  mask=(first['times']>=5*w)&((first['times']<5*(w+1)) if w<3 else first['times']<=20);e.require(mask.sum() in [10,11],'Window frame count');feats=[];joints=[];lefts=[];rights=[]
  for z in packets:
   complex=z['complex_sasa_A2'][mask];burial=z['burial_A2'][mask];dist=z['min_e2_A'][mask];reach=z['reach80_fraction'][mask];bridge=z['bridge'][mask,None];feats.append(np.stack([np.ones(lyscount),complex.mean(0)/150,burial.mean(0)/150,dist.mean(0)/100,reach.mean(0),((complex>5)&(dist<80)&bridge).mean(0)],axis=1));l=z['vhl_left'][mask];r=z['vhl_right'][mask];joints.append(np.einsum('ti,tj->ij',(l<5).astype(float),(r<5).astype(float))/len(l));lefts.append(l.mean(0));rights.append(r.mean(0))
  x=base.copy();x[nl+first['lys_node_indices'],23:]=np.mean(feats,axis=0);joint=np.mean(joints,axis=0);left=np.mean(lefts,axis=0);right=np.mean(rights,axis=0);i,j=np.nonzero(joint>=.01);ei=np.array([np.r_[i,nl+j],np.r_[nl+j,i]],dtype=np.int64);ea=np.column_stack([joint[i,j],left[i]/10,right[j]/10]).astype(np.float32);ea=np.concatenate([ea,ea]);e.require(np.isfinite(x).all() and np.isfinite(ea).all(),'Nonfinite graph');arr.update({f'x{w}':x,f'edge_index{w}':ei,f'edge_attr{w}':ea})
 return arr

def main():
 e.freeze();rows=json.loads((ROOT/'extraction_completion.json').read_text());labels=json.loads((e.NATIVE/'config/labels.json').read_text())['labels'];checks=[]
 for r in rows:
  t,pose=r['target'],r['pose']
  if t not in labels:continue
  base=np.load(ROOT/'extracted'/t/pose/'base_packet.npz');old=np.load(e.MECH/'analysis/lysines/poses'/t/pose/'packet.npz');errs={}
  for key in old.files:
   if old[key].dtype.kind in 'fc':
    error=float(np.max(np.abs(base[key]-old[key]))) if old[key].size else 0.;e.require(error<.002,'Base descriptor mismatch '+t+' '+key+' '+str(error));errs['base_'+key]=error
   else:e.require(np.array_equal(base[key],old[key]),'Base identity mismatch '+key)
  new=np.load(ROOT/'extracted'/t/pose/'native_packet.npz');old=np.load(e.NATIVE/'results/poses'/t/pose/'packet.npz')
  for key in ['times_ns','donor_distance_A','nz_exposed_fraction','target_scaffold_max_overlap_A','anchor_RMSD_A','sidechain_outward_cosine']:
   error=float(np.max(np.abs(new[key]-old[key])));e.require(error<.002,'Native descriptor mismatch '+key+' '+str(error));errs['native_'+key]=error
  e.require(np.array_equal(new['bridge'],old['bridge']),'Bridge mismatch');checks.append({'target':t,'pose':pose,'errors':errs})
 e.require(len(checks)==8,'Eight controls required');manifest=[]
 for t in ['SIK3_full','RPS6KB2_full']:
  selected=sorted([r for r in rows if r['target']==t],key=lambda x:x['pose']);basepackets=[dict(np.load(ROOT/'extracted'/t/r['pose']/'base_packet.npz')) for r in selected];nativepackets=[dict(np.load(ROOT/'extracted'/t/r['pose']/'native_packet.npz')) for r in selected];arr=base_graph(basepackets);nl=int(arr['x0'][:,0].sum());lys_indices=basepackets[0]['lys_node_indices'];variants={'native_features':{},'zero_padded_control':{}}
  for w in range(4):
   values=[]
   for p in nativepackets:
    tt=p['times_ns'];mask=(tt>=5*w)&((tt<5*(w+1)) if w<3 else tt<=20);e.require(mask.sum()==(50 if w<3 else 51),'Native window count');d=p['donor_distance_A'][0,mask,0,:];ex=p['nz_exposed_fraction'][mask];clear=p['target_scaffold_max_overlap_A'][0,mask]<=.8;bridge=p['bridge'][mask];joint=(d<=20)&(ex>=.1)&clear[:,None];values.append(np.stack([d.mean(0)/100,(d<=20).mean(0),ex.mean(0),joint.mean(0),(joint&bridge[:,None]).mean(0)],axis=1))
   addition=np.zeros((len(arr[f'x{w}']),5),np.float32);addition[nl+lys_indices]=np.mean(values,axis=0)
   for v,a in [('native_features',addition),('zero_padded_control',np.zeros_like(addition))]:
    variants[v][f'x{w}']=np.c_[arr[f'x{w}'],a].astype(np.float32)
    for k in ['edge_index','edge_attr']:variants[v][f'{k}{w}']=arr[f'{k}{w}'].copy()
  for v,arrays in variants.items():
   d=ROOT/'inputs'/v;d.mkdir(parents=True,exist_ok=True);path=d/(t+'.npz');np.savez_compressed(path,**arrays);manifest.append({'target':t,'label':None,'variant':v,'sha256':e.sha(path),'path':str(path.relative_to(ROOT)),'poses':len(selected),'frames_per_pose':201,'node_width':34})
 e.save(ROOT/'validation/feature_matching.json',{'status':'PASS','same_definitions':True,'eight_original_targets_reextracted':checks,'unknown_graphs':manifest,'training_performed':False});print('PASS:8 matching controls;4 unknown graphs prepared')
if __name__=='__main__':main()
