"""Analytic controls and explicit MZ1 native-reference projection."""
from pathlib import Path
import sys,json,importlib.util
import numpy as np
from Bio.PDB import MMCIFParser
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('extract',Path(__file__).with_name('02_extract_md_compatibility.py'));e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def main():
 # Surface and threshold controls do not use target labels.
 assert np.allclose(e.exposed_fraction(np.array([[0.,0,0]]),np.array([0]),np.array([[0.,0,0]]),np.array([0]),np.array([1.55])),1)
 surrounding=e.SPHERE*2.5;xyz=np.vstack([np.zeros((1,3)),surrounding]);assert e.exposed_fraction(np.zeros((1,3)),np.array([0]),xyz,np.arange(len(xyz)),np.full(len(xyz),1.7))[0]==0
 a=np.array([[0.,0,0],[1.,0,0],[0,2.,0],[0,0,3.]]);angle=.73;r=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]]);b=a@r.T+np.array([20,-5,7]);point=np.array([[3.,4.,5.]]);mapped,err,rr=e.fit_reference_points(a,b,point);assert err<1e-10 and np.allclose(mapped,point@r.T+[20,-5,7]);assert np.allclose((mapped-b.mean(0))@rr+a.mean(0),point)
 p=MMCIFParser(QUIET=True,auth_chains=False);native=p.get_structure('n',str(ROOT/'references/5T35.cif'))[0];closed=p.get_structure('c',str(ROOT/'references/8RX0.cif'))[0];pack=np.load(ROOT/'references/closed_observed.npz');rows=[]
 for targetchain,vhlchain in [('A','D'),('E','H')]:
  common=sorted(set(r.id[1] for r in native[vhlchain] if r.id[0]==' ' and 'CA' in r)&set(r.id[1] for r in closed['B'] if r.id[0]==' ' and 'CA' in r));mobile=np.array([native[vhlchain][(' ',i,' ')]['CA'].coord for i in common]);ref=np.array([closed['B'][(' ',i,' ')]['CA'].coord for i in common]);lys=[r for r in native[targetchain] if r.id[0]==' ' and 'NZ' in r and 'CE' in r];nz=np.array([r['NZ'].coord for r in lys]);ce=np.array([r['CE'].coord for r in lys]);mapped,rmsd,rotation=e.fit_reference_points(mobile,ref,np.vstack([nz,ce]));z=mapped[:len(lys)];cc=mapped[len(lys):];vec=pack['donor_C'][:,None,:]-z[None,:,:];dd=np.linalg.norm(vec,axis=2);cos=np.einsum('eli,li->el',vec,z-cc)/np.maximum(dd*np.linalg.norm(z-cc,axis=1)[None,:],1e-9)
  for j,res in enumerate(lys):rows.append({'crystal_copy':targetchain,'resid':res.id[1],'VHL_fit_RMSD_A':rmsd,'distance_to_two_modeled_G76C_A':dd[:,j].tolist(),'outward_cosines':cos[:,j].tolist()})
  far=z+np.array([500,0,0]);assert np.linalg.norm(pack['donor_C'][:,None,:]-far[None,:,:],axis=2).min()>30
 result={'software_controls':'PASS','isolated_NZ_exposure':1,'buried_NZ_exposure':0,'proper_rotation_and_inverse':'PASS','translated_far_negative':'PASS','native_MZ1_projection':rows,'native_control_is_supplier_SJF8240':False,'endpoint_observed':False,'known_positive_primary_distance_orientation_pass':any(min(r['distance_to_two_modeled_G76C_A'])<20 and any(d<20 and co>=0 for d,co in zip(r['distance_to_two_modeled_G76C_A'],r['outward_cosines'])) for r in rows),'known_K456_primary_distance_pass':any(r['resid']==456 and min(r['distance_to_two_modeled_G76C_A'])<20 for r in rows),'note':'Native MZ1 structural projection tests geometry and endpoint sensitivity, not a calibration to SJF8240 labels. K456 is tracked explicitly; missing donorG76 means reported distances are modeled.'};(ROOT/'validation/controls.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps([r for r in rows if r['resid'] in [368,445,456]],indent=2))
if __name__=='__main__':main()
