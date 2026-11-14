"""Explicit observed/modeled CRL2 references. No outcome-driven scaffold fitting."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fit(a,b,points):
 ac=a.mean(0);bc=b.mean(0);u,s,vt=np.linalg.svd((a-ac).T@(b-bc));r=vt.T@np.diag([1,1,np.linalg.det(vt.T@u.T)])@u.T;return (points-ac)@r.T+bc,float(np.sqrt(np.mean(np.sum(((a-ac)@r.T+bc-b)**2,axis=1))))
def atoms(chain):return [a for r in chain if r.id[0]==' ' for a in r if a.element not in ['H','D']]
def named(chain,resids,names):return np.array([chain[(' ',i,' ')][n].coord for i in resids for n in names],float)
def main():
 p=MMCIFParser(QUIET=True,auth_chains=False);closed=p.get_structure('closed',str(ROOT/'references/8RX0.cif'))[0];opened=p.get_structure('open',str(ROOT/'references/8RWZ.cif'))[0];ub=p.get_structure('ub',str(ROOT/'references/1UBQ.cif'))[0]['A'];records=[];endpoint=[]
 for start in [70,72,73]:
  a=named(ub,range(start,76),['N','CA','C']);b=named(closed['I'],range(start,76),['N','CA','C']);full=np.array([ub[(' ',76,' ')][n].coord for n in ['N','CA','C','O']],float);coords,error=fit(a,b,full);peptide=float(np.linalg.norm(coords[0]-closed['I'][(' ',75,' ')]['C'].coord));endpoint.append({'name':'ub_tail_fit_'+str(start),'xyz':coords.tolist(),'fit_RMSD_A':error,'G75C_G76N_distance_A':peptide,'modeled':True})
 # Retain all modeled tails, but only peptide-consistent tails can enter geometric scores.
 valid=[x for x in endpoint if 1.1<=x['G75C_G76N_distance_A']<=1.6];assert valid,'No valid reconstructed terminal reference'
 for label,model in [('closed_observed',closed),('open_RBX1_transferred_module',opened)]:
  common=sorted(set(r.id[1] for r in closed['H'] if r.id[0]==' ')&set(r.id[1] for r in opened['F'] if r.id[0]==' '))
  if label=='closed_observed':transform=lambda xyz:xyz;modulefit=0.
  else:
   aa=named(closed['H'],common,['CA']);bb=named(opened['F'],common,['CA']);transform=lambda xyz:fit(aa,bb,np.asarray(xyz,float))[0];modulefit=fit(aa,bb,aa)[1]
  points=[];elements=[];roles=[]
  # New components only: shared VHL/ELOB/ELOC stay as actually simulated.
  chain_groups=[(model['E'],'CUL2'),(model['H'] if label=='closed_observed' else model['F'],'RBX1')]
  for chain,role in chain_groups:
   for a in atoms(chain):points.append(a.coord);elements.append(a.element);roles.append(role)
  for cc,role in [('F','UBE2R1_mutant'),('G','NEDD8'),('I','donor_Ub_incomplete')]:
   for a in atoms(closed[cc]):points.append(transform(a.coord));elements.append(a.element);roles.append(role)
  points=np.array(points,float);np.savez_compressed(ROOT/'references'/f'{label}.npz',xyz=points,elements=np.array(elements),roles=np.array(roles),VHL_CA=np.array([r['CA'].coord for r in model['B'] if r.id[0]==' ' and 'CA' in r]),VHL_sequence=np.array(''.join(seq1(r.resname) for r in model['B'] if r.id[0]==' ' and 'CA' in r)),ELOC_CA=np.array([r['CA'].coord for r in model['C' if label=='closed_observed' else 'D'] if r.id[0]==' ' and 'CA' in r]),ELOC_sequence=np.array(''.join(seq1(r.resname) for r in model['C' if label=='closed_observed' else 'D'] if r.id[0]==' ' and 'CA' in r)),donor_C=np.array([transform(np.array(x['xyz']))[2] for x in valid]),donor_tail_names=np.array([x['name'] for x in valid]));records.append({'name':label,'module_RBX1_fit_RMSD_A':modulefit,'model_based_module_transfer':label!='closed_observed','heavy_atoms':len(points),'endpoint_count':len(valid),'file_sha256':sha(ROOT/'references'/f'{label}.npz')})
 (ROOT/'validation/scaffold_preparation.json').write_text(json.dumps({'sources_sha256':{p.name:sha(p) for p in (ROOT/'references').glob('*.cif')},'observed_missing':'8RX0 donor ubiquitin modeled only through Gly75; E2 and several catalytic components lack sidechain heavy atoms. Open scaffold lacks resolved E2/Ub/NEDD8 module.','endpoint_model':'1UBQ terminal backbone fitted to8RX0 Ub locally, retained only when G75C-G76N1.1–1.6A; all terminal estimates are modeled, not observed catalytic atoms.','endpoint_attempts':endpoint,'scaffolds':records,'atomic_sterics_complete':False,'reaction_probability':False},indent=2)+'\n');print(json.dumps({'endpoints':endpoint,'scaffolds':records},indent=2))
if __name__=='__main__':main()
