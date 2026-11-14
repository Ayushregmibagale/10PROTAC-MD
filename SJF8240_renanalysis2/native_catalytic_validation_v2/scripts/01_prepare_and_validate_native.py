"""Native-backbone-aware controls and deposited catalytic-center surrogate references."""
from pathlib import Path
import sys,json,hashlib,importlib.util
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];PRIOR=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2/crl2_ubiquitination_v1');spec=importlib.util.spec_from_file_location('prior_extract',PRIOR/'scripts/02_extract_md_compatibility.py');e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def heavy(chain):return [a for r in chain if r.id[0]==' ' for a in r if a.element not in ['H','D']]
def xyz(chain):return np.array([a.coord for a in heavy(chain)],float)
def fit(a,b,pts):return e.fit_reference_points(np.array(a,float),np.array(b,float),np.array(pts,float))
def refpack(model,chains,endpoint,path):
 points=[];els=[];roles=[]
 for ch,role in chains:
  for a in heavy(model[ch]):points.append(a.coord);els.append(a.element);roles.append(role)
 vch='H' if path.stem.startswith('8R5H') else 'B';ech='E' if path.stem.startswith('8R5H') else 'C';vhl=[r for r in model[vch] if r.id[0]==' ' and 'CA' in r];eloc=[r for r in model[ech] if r.id[0]==' ' and 'CA' in r];np.savez_compressed(path,xyz=np.array(points,float),elements=np.array(els),roles=np.array(roles),VHL_CA=np.array([r['CA'].coord for r in vhl]),VHL_sequence=np.array(''.join(seq1(r.resname) for r in vhl)),ELOC_CA=np.array([r['CA'].coord for r in eloc]),ELOC_sequence=np.array(''.join(seq1(r.resname) for r in eloc)),donor_C=np.array([endpoint]),donor_tail_names=np.array(['deposited_SY8_C3_surrogate' if 'SY8' in path.stem else 'modeled_G76C_fit75']))
def rebuild_target(backbone,crystal):
 points=[];elements=[];names=[];resids=[];receipts=[]
 for r in backbone:
  if r.id[0]!=' ':continue
  observed={a.name:a for a in r if a.element not in ['H','D']}
  for a in observed.values():points.append(a.coord);elements.append(a.element);names.append(a.name);resids.append(r.id[1])
  if r.id not in crystal:continue
  c=crystal[r.id]
  if c.resname!=r.resname or any(n not in r or n not in c for n in ['N','CA','C']):continue
  missing=[a for a in c if a.element not in ['H','D'] and a.name not in observed];aa=np.array([c[n].coord for n in ['N','CA','C']]);bb=np.array([r[n].coord for n in ['N','CA','C']]);mapped,rmsd,rot=fit(aa,bb,np.array([a.coord for a in missing]).reshape(-1,3) if missing else aa)
  if rmsd>.25:continue
  for a,p in zip(missing,mapped):points.append(p);elements.append(a.element);names.append(a.name);resids.append(r.id[1])
  receipts.append({'resid':r.id[1],'local_backbone_RMSD_A':rmsd,'transferred_sidechain_atoms':[a.name for a in missing]})
 return np.array(points,float),np.array(elements),np.array(names),np.array(resids),receipts

def main():
 parser=MMCIFParser(QUIET=True,auth_chains=False);r1=parser.get_structure('r1',str(PRIOR/'references/8RX0.cif'))[0];r2=parser.get_structure('r2',str(ROOT/'references/8R5H.cif'))[0];xtal=parser.get_structure('xtal',str(PRIOR/'references/5T35.cif'))[0];lig=next(iter(r2['L']));ccd=MMCIF2Dict(str(ROOT/'references/SY8.cif'));bondlist=list(zip(ccd['_chem_comp_bond.atom_id_1'],ccd['_chem_comp_bond.atom_id_2'],ccd['_chem_comp_bond.value_order']));assert any({a,b}=={'N','C4'} for a,b,_ in bondlist) and any({a,b}=={'C4','C3'} for a,b,_ in bondlist);assert abs(np.linalg.norm(lig['C3'].coord-r2['C'][(' ',93,' ')]['SG'].coord)-1.7)<.15
 connections=MMCIF2Dict(str(ROOT/'references/8R5H.cif'));conn=[]
 for i,kind in enumerate(connections['_struct_conn.conn_type_id']):
  if kind=='covale':conn.append({k.split('.')[-1]:v[i] for k,v in connections.items() if k.startswith('_struct_conn.') and ('label' in k or 'dist_value' in k)})
 endpoint=lig['C3'].coord.astype(float);# stable-crosslink surrogate of the donor terminus,not a native carbonyl atom
 original=np.load(PRIOR/'references/closed_observed.npz');g76=original['donor_C'][1];localids=[i for i in range(85,101) if (' ',i,' ') in r1['F'] and (' ',i,' ') in r2['C']];aa=np.array([r2['C'][(' ',i,' ')][n].coord for i in localids for n in ['N','CA','C']]);bb=np.array([r1['F'][(' ',i,' ')][n].coord for i in localids for n in ['N','CA','C']]);mapped,efit,_=fit(aa,bb,np.array([endpoint,r2['C'][(' ',93,' ')]['SG'].coord]));assert efit<=1.
 names=['8R5H_observed_SY8','8RX0_G76_fit75','8RX0_E2_aligned_SY8'];refpack(r2,[('A','CUL2'),('B','RBX1'),('C','UBE2R2'),('D','donor_Ub')],endpoint,ROOT/'references'/f'{names[0]}.npz');refpack(r1,[('E','CUL2'),('H','RBX1'),('F','UBE2R1'),('G','NEDD8'),('I','donor_Ub')],g76,ROOT/'references'/f'{names[1]}.npz');refpack(r1,[('E','CUL2'),('H','RBX1'),('F','UBE2R1'),('G','NEDD8'),('I','donor_Ub')],mapped[0],ROOT/'references'/f'{names[2]}.npz');controls=[];transfers=[]
 for name,model,targetchain,shared,donor in [(names[0],r2,'F',['E','G','H'],endpoint),(names[1],r1,'A',['B','C','D'],g76),(names[2],r1,'A',['B','C','D'],mapped[0])]:
  ref=np.load(ROOT/'references'/f'{name}.npz');sr=np.array([e.RAD[x] for x in ref['elements']]);tree=cKDTree(ref['xyz'])
  for cry in ['A','E']:
   target,els,anames,resids,evidence=rebuild_target(model[targetchain],xtal[cry]);transfers.append({'template':name,'crystal_copy':cry,'residues':evidence});radii=np.array([e.RAD[x] for x in els]);oc=np.vstack([target,*[xyz(model[ch]) for ch in shared]]);oids=np.arange(len(oc));orr=np.concatenate([radii,*[np.array([e.RAD[a.element] for a in heavy(model[ch])]) for ch in shared]]);nzindices=np.flatnonzero(anames=='NZ');nz=target[nzindices];ce=np.array([target[np.flatnonzero((resids==resids[i])&(anames=='CE'))[0]] for i in nzindices]);ex=e.exposed_fraction(nz,nzindices,oc,oids,orr);dist=np.linalg.norm(donor-nz,axis=1);vec=nz-ce;cos=np.einsum('ij,ij->i',donor-nz,vec)/(dist*np.linalg.norm(vec,axis=1));ov,_=e.overlap_positive(target,radii,tree,ref['xyz'],sr)
   for j,i in enumerate(nzindices):controls.append({'template':name,'crystal_copy':cry,'lysine':int(resids[i]),'NZ_distance_A':float(dist[j]),'NZ_exposed_fraction':float(ex[j]),'outward_cosine':float(cos[j]),'whole_target_scaffold_overlap_A':ov,'within20A':bool(dist[j]<=20),'Bai16A_geometry':bool(dist[j]<=16 and ex[j]>=.1 and ov<=.8),'primary20A_geometry':bool(dist[j]<=20 and ex[j]>=.1 and ov<=.8),'prior_oriented20A_rule':bool(dist[j]<=20 and ex[j]>=.1 and ov<=.8 and cos[j]>=0)})
  # Negative control is moved out of reach by500A; binding geometry is deliberately destroyed.
  assert np.linalg.norm(target+np.array([500,0,0])-donor,axis=1).min()>30
 np.savez_compressed(ROOT/'validation/E2_endpoint_transfer.npz',E2R2_backbone=aa,E2R1_backbone=bb,SY8_C3=endpoint,mapped_SY8_C3=mapped[0],mapped_E2_SG=mapped[1]);import pandas as pd;pd.DataFrame(controls).to_csv(ROOT/'results/native_controls.tsv',sep='\t',index=False)
 native={'sources_sha256':{str(p):sha(p) for p in [ROOT/'references/8R5H.cif',ROOT/'references/SY8.cif',PRIOR/'references/8RX0.cif',PRIOR/'references/5T35.cif']},'SY8_chemical_bonds':bondlist,'SY8_covalent_links':conn,'E2_active_loop_fit_RMSD_A':efit,'deposited_endpoint':'8R5H ligandSY8 atomC3; N-C4-C3 connectivity and covalent C93SG-C3 verified; stable-crosslink surrogate,not native Gly76 carbonyl','endpoint_types':names,'missing_components':'8R5H has no modeled NEDD8 and lacks some density/sidechains;8RX0 catalytic module is backbone-only. Nothing is called complete atomic chemistry.','native_sidechain_method':'Observed catalytic-target backbone retained; missing heavy sidechains transferred through residue-specific N,CA,C fits from5T35 copiesA/E,matching amino acid only and RMSD<=0.25A. No donor-directed rotamer search.','earlier_control_issue':'Global VHL-superposition of crystal BRD4 retained the crystal distal backbone/rotamers instead of catalytic BRD4 geometry. New controls honor deposited target backbone.','mutated_Cys368':'Kept as depositedCYS,never scored as wild-type lysine.','sidechain_transfers':transfers,'controls':'native_controls.tsv','biological_validation':False};(ROOT/'validation/native_reference_audit.json').write_text(json.dumps(native,indent=2)+'\n');print(pd.DataFrame(controls).query('lysine in [368,445,456]').to_string(index=False))
if __name__=='__main__':main()
