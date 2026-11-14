"""Model missing donorGly76 using only peptide and E2-thioester chemistry.
No target lysine coordinate or degradation label enters the optimization.
"""
from pathlib import Path
import sys,json,importlib.util
import numpy as np
from scipy.optimize import least_squares
from Bio.PDB import MMCIFParser
ROOT=Path(__file__).resolve().parents[1];PRIOR=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2/crl2_ubiquitination_v1');spec=importlib.util.spec_from_file_location('prep',ROOT/'scripts/01_prepare_and_validate_native.py');e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def angle(a,b,c):
 u=a-b;v=c-b;return np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1)))
def dihedral(a,b,c,d):
 b0=-(b-a);b1=c-b;b2=d-c;b1=b1/np.linalg.norm(b1);v=b0-np.dot(b0,b1)*b1;w=b2-np.dot(b2,b1)*b1;return np.degrees(np.arctan2(np.dot(np.cross(b1,v),w),np.dot(v,w)))
def main():
 parser=MMCIFParser(QUIET=True,auth_chains=False);r2=parser.get_structure('r2',str(ROOT/'references/8R5H.cif'))[0];r1=parser.get_structure('r1',str(PRIOR/'references/8RX0.cif'))[0];trans=np.load(ROOT/'validation/E2_endpoint_transfer.npz');lig=next(iter(r2['L']));inputs=[('8R5H',r2['D'],r2['C'][(' ',93,' ')]['SG'].coord,np.array([lig[n].coord for n in ['N','C4','C3']])),('8RX0',r1['I'],trans['mapped_E2_SG'],None)];results=[];rng=np.random.default_rng(42)
 for name,ub,sg,seed in inputs:
  g75=ub[(' ',75,' ')];prevca=g75['CA'].coord.astype(float);prevc=g75['C'].coord.astype(float);prevo=g75['O'].coord.astype(float);sg=np.array(sg,float)
  if seed is None:
   # Existing local fit supplies a starting conformation only; final covalent checks are independent.
   attempt=json.loads((PRIOR/'validation/scaffold_preparation.json').read_text())['endpoint_attempts'];seed=np.array(next(a['xyz'][:3] for a in attempt if a['name']=='ub_tail_fit_75'))
  n,ca,c=seed;u=(ca-c)/np.linalg.norm(ca-c);v=(sg-c)/np.linalg.norm(sg-c);o=c-1.23*(u+v)/np.linalg.norm(u+v);initial=np.vstack([seed,o]);trials=[]
  def residual(flat):
   n,ca,c,o=flat.reshape(4,3);lengths=[np.linalg.norm(prevc-n)-1.33,np.linalg.norm(n-ca)-1.45,np.linalg.norm(ca-c)-1.52,np.linalg.norm(c-o)-1.23,np.linalg.norm(c-sg)-1.80];angles=[angle(prevca,prevc,n)-116,angle(prevc,n,ca)-122,angle(n,ca,c)-111,angle(ca,c,o)-121,angle(ca,c,sg)-116,angle(o,c,sg)-123];omega=dihedral(prevca,prevc,n,ca);omega_err=abs(abs(omega)-180);plane=np.dot(o-c,np.cross(ca-c,sg-c));previous_plane=np.dot(n-prevc,np.cross(prevca-prevc,prevo-prevc));previous_angle=angle(prevo,prevc,n)-123;return np.r_[np.array(lengths)/.02,np.array(angles)/3,omega_err/5,plane/.02,previous_plane/.02,previous_angle/3,(flat-initial.ravel())*.01]
  best=None
  for attempt in range(8):
   start=initial.ravel().copy()
   if attempt:start+=rng.normal(0,.3,start.shape)
   fit=least_squares(residual,start,max_nfev=3000,ftol=1e-10,xtol=1e-10,gtol=1e-10);x=fit.x.reshape(4,3);n,ca,c,o=x;bond=np.array([np.linalg.norm(prevc-n),np.linalg.norm(n-ca),np.linalg.norm(ca-c),np.linalg.norm(c-o),np.linalg.norm(c-sg)]);ang=np.array([angle(prevca,prevc,n),angle(prevc,n,ca),angle(n,ca,c),angle(ca,c,o),angle(ca,c,sg),angle(o,c,sg)]);bd=float(np.max(abs(bond-[1.33,1.45,1.52,1.23,1.80])));ad=float(np.max(abs(ang-[116,122,111,121,116,123])));omega=dihedral(prevca,prevc,n,ca);valid=bool(fit.success and bd<=.04 and ad<=6 and abs(abs(omega)-180)<=10 and abs(angle(prevo,prevc,n)-123)<=6 and abs(np.dot(n-prevc,np.cross(prevca-prevc,prevo-prevc)))<=.1);record={'attempt':attempt,'optimizer_success':bool(fit.success),'bond_lengths_A':bond.tolist(),'angles_deg':ang.tolist(),'omega_deg':omega,'previous_O_C_N_angle_deg':angle(prevo,prevc,n),'previous_carbonyl_plane_triple_A3':float(np.dot(n-prevc,np.cross(prevca-prevc,prevo-prevc))),'max_bond_error_A':bd,'max_angle_error_deg':ad,'covalent_geometry_valid':valid,'objective':float(np.sum(residual(fit.x)**2))};trials.append(record);key=(not valid,record['objective'])
   if best is None or key<best[0]:best=(key,x,record)
  _,x,record=best;results.append({'reference':name,'G76_atoms_N_CA_C_O':x.tolist(),'E2_C93_SG':sg.tolist(),'chemistry_checks':record,'all_trials':trials,'target_coordinates_used':False,'target_labels_used':False,'endpoint_observed':False,'description':'Constrained modeled native Gly76 thioester with fixedG75backbone andSG. Sidechain transfer for8RX0SG is modeled. No steric/energetic or reaction validation inferred.'})
  np.savez_compressed(ROOT/'references'/f'{name}_covalent_G76_model.npz',G76_xyz=x,E2_SG=sg,G75_CA=prevca,G75_C=prevc)
 (ROOT/'validation/covalent_tail_models.json').write_text(json.dumps({'models':results,'label_independent':True,'chemical_geometry_only':True},indent=2)+'\n');print(json.dumps([{r['reference']:r['chemistry_checks']} for r in results],indent=2))
if __name__=='__main__':main()
