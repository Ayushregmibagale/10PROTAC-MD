"""Bounded donor-tail repair allowing onlyUbG75O plus missingG76 to move.
E2/Ub backbone andSG remain fixed; no target atoms/labels used.
"""
from pathlib import Path
import json,sys,importlib.util
import numpy as np
from scipy.optimize import least_squares
from Bio.PDB import MMCIFParser
ROOT=Path(__file__).resolve().parents[1];PRIOR=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2/crl2_ubiquitination_v1');spec=importlib.util.spec_from_file_location('tail',ROOT/'scripts/09_model_covalent_donor_tail.py');t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t);RAD={'C':1.7,'N':1.55,'O':1.52,'S':1.8}
def main():
 parser=MMCIFParser(QUIET=True,auth_chains=False);rng=np.random.default_rng(42);results=[]
 for code in ['8R5H','8RX0']:
  s=parser.get_structure(code,str(ROOT/'references/8R5H.cif' if code=='8R5H' else PRIOR/'references/8RX0.cif'))[0];ec,uc=('C','D') if code=='8R5H' else ('F','I');g75=s[uc][(' ',75,' ')];pc=g75['C'].coord.astype(float);pa=g75['CA'].coord.astype(float);po=g75['O'].coord.astype(float);z=np.load(ROOT/'references'/f'{code}_covalent_G76_model.npz');sg=z['E2_SG'];initial=np.vstack([po,z['G76_xyz']]);core=[]
  for ch in [ec,uc]:
   for r in s[ch]:
    if r.id[0]!=' ':continue
    for a in r:
     if a.element in ['H','D'] or (ch==uc and r.id[1]==75 and a.name=='O'):continue
     core.append((ch,r.id[1],a.name,a.element,a.coord.astype(float)))
  if code=='8RX0':core.append((ec,93,'SG','S',sg))
  # Fixed8A neighborhood plus subsequent global replay; optimizer cannot shift>3A peraxis.
  core=[r for r in core if np.linalg.norm(r[4]-initial.mean(0))<12.];cx=np.array([r[4] for r in core]);cr=np.array([RAD[r[3]] for r in core]);nr=np.array([1.52,1.55,1.7,1.7,1.52]);excluded={(0,uc,75,'C'),(0,uc,75,'CA'),(1,uc,75,'C'),(1,uc,75,'CA'),(2,uc,75,'C'),(3,ec,93,'SG'),(3,ec,93,'CB'),(2,ec,93,'SG'),(4,ec,93,'SG')};mask=np.array([[(i,ch,res,name) not in excluded for ch,res,name,_,_ in core] for i in range(5)]);pair_i=np.array([0,0,0,1]);pair_j=np.array([2,3,4,4])
  def terms(flat):
   op,n,ca,c,o=flat.reshape(5,3);length=np.array([np.linalg.norm(pc-op),np.linalg.norm(pc-n),np.linalg.norm(n-ca),np.linalg.norm(ca-c),np.linalg.norm(c-o),np.linalg.norm(c-sg)]);ang=np.array([t.angle(pa,pc,op),t.angle(pa,pc,n),t.angle(op,pc,n),t.angle(pc,n,ca),t.angle(n,ca,c),t.angle(ca,c,o),t.angle(ca,c,sg),t.angle(o,c,sg)]);omega=t.dihedral(pa,pc,n,ca);plane1=np.dot(n-pc,np.cross(pa-pc,op-pc));plane2=np.dot(o-c,np.cross(ca-c,sg-c));return length,ang,omega,plane1,plane2
  def residual(flat):
   points=flat.reshape(5,3);length,ang,omega,plane1,plane2=terms(flat);over=nr[:,None]+cr[None,:]-np.linalg.norm(points[:,None,:]-cx[None,:,:],axis=2);internal=nr[pair_i]+nr[pair_j]-np.linalg.norm(points[pair_i]-points[pair_j],axis=1);return np.r_[(length-[1.23,1.33,1.45,1.52,1.23,1.80])/.02,(ang-[121,116,123,122,111,121,116,123])/3,(abs(abs(omega)-180))/5,plane1/.02,plane2/.02,np.maximum(over-.6,0)[mask]/.15,np.maximum(internal-.6,0)/.15,(points[0]-po)*.2,(points[1:]-initial[1:]).ravel()*.01]
  trials=[];best=None
  for attempt in range(8):
   start=initial.copy()
   if attempt:start+=rng.normal(0,.25,start.shape)
   lower=(initial-3.).ravel();upper=(initial+3.).ravel();fit=least_squares(residual,np.clip(start.ravel(),lower+1e-5,upper-1e-5),bounds=(lower,upper),max_nfev=1500,ftol=1e-9,xtol=1e-9,gtol=1e-9);x=fit.x.reshape(5,3);length,ang,omega,plane1,plane2=terms(fit.x);be=float(max(abs(length-[1.23,1.33,1.45,1.52,1.23,1.80])));ae=float(max(abs(ang-[121,116,123,122,111,121,116,123])));over=nr[:,None]+cr[None,:]-np.linalg.norm(x[:,None,:]-cx[None,:,:],axis=2);maxov=float(max(over[mask].max(),(nr[pair_i]+nr[pair_j]-np.linalg.norm(x[pair_i]-x[pair_j],axis=1)).max()));valid=bool(fit.success and be<=.04 and ae<=6 and abs(abs(omega)-180)<=10 and maxov<=.8 and abs(plane1)<=.1 and abs(plane2)<=.1);record={'attempt':attempt,'optimizer_success':bool(fit.success),'bond_lengths_A':length.tolist(),'angles_deg':ang.tolist(),'omega_deg':float(omega),'carbonyl_planes_A3':[float(plane1),float(plane2)],'max_bond_error_A':be,'max_angle_error_deg':ae,'maximum_checked_nonbonded_overlap_A':maxov,'G75O_displacement_A':float(np.linalg.norm(x[0]-po)),'validated_geometry_only':valid,'objective':float(np.sum(residual(fit.x)**2))};trials.append(record);rank=(not valid,record['objective'])
   if best is None or rank<best[0]:best=(rank,x,record)
  _,x,record=best;np.savez_compressed(ROOT/'references'/f'{code}_locally_relaxed_G76.npz',G75_O=x[0],G76_xyz=x[1:],E2_SG=sg);results.append({'reference':code,'record':record,'G75O_G76_points':x.tolist(),'all_trials':trials,'target_atoms_or_labels_used':False,'scope':'OnlyG75O andmissingG76 atoms allowed3A perCartesianaxis;otherUb/E2core andSG fixed. Geometric repair,not energy-minimized reaction.'});print(code,record,flush=True)
 (ROOT/'validation/locally_relaxed_tail_models.json').write_text(json.dumps({'models':results,'no_target_or_label_fit':True,'full_atomic_model_complete':False},indent=2)+'\n')
if __name__=='__main__':main()
