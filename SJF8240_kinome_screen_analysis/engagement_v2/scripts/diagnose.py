from pathlib import Path
import numpy as np
from rdkit import Chem
PARAM=Path('/home/ayushregmibagale/kinome_scan/data/kinase_pRosettaC/Q16832/Patchdock_Results/PT_1_14.params')
def mol(path):
 rw=Chem.RWMol();names=[];bonds=[]
 for line in path.read_text().splitlines():
  s=line.split()
  if s and s[0]=='ATOM':
   name=s[1];el=''.join(c for c in name if c.isalpha());a=Chem.Atom(el);rw.AddAtom(a);names.append(name)
  if s and s[0]=='BOND_TYPE':bonds.append(s[1:])
 for a,b,o in bonds:rw.AddBond(names.index(a),names.index(b),{'1':Chem.BondType.SINGLE,'2':Chem.BondType.DOUBLE,'3':Chem.BondType.TRIPLE,'ar':Chem.BondType.AROMATIC}[o])
 m=rw.GetMol();Chem.SanitizeMol(m);return m,names,bonds
m,n,b=mol(PARAM);print('full',len(n),Chem.MolToSmiles(m))
for h in ['PT0','PT1']:
 p=Path('/mnt/d/kinome_results/P42681');hm,hn,hb=mol(p/(h+'.params'));matches=m.GetSubstructMatches(hm,uniquify=False);print(h,'maps',len(matches),'atom sets',len(set(tuple(sorted(x)) for x in matches)))
 coords={l[12:16].strip():np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]) for l in (p/(h+'_0001.pdb')).read_text().splitlines() if l.startswith(('ATOM','HETATM'))}
 lengths=[(a,b,float(np.linalg.norm(coords[a]-coords[b]))) for a,b,o in hb];print('reference bond range',min(x[2] for x in lengths),max(x[2] for x in lengths));print('bad',[(a,b,round(d,2)) for a,b,d in lengths if d>2.1 or d<1])
 print('map example',[(hn[i],n[j]) for i,j in enumerate(matches[0])][:8] if matches else '')
