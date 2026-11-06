from pathlib import Path
import json,sys,time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
r=Path('/tmp/10protac_finish/prosettac_robust');w=r/'workflow';sys.path[:0]=[str(w/'baseline_engine'),str(w/'scripts/lib')]
from constrained_linker import constrained_linker
from official_parameters import docked_heads
p=r/'validation_controls/sjf8240_global50_diagnostic_v2/engine_raw';hp=[p/'vhl_head_H.sdf',p/'poi_head_H.sdf'];heads=[Chem.SDMolSupplier(str(f),removeHs=True)[0] for f in hp];full=Chem.MolFromSmiles((p/'full.smi').read_text().split()[0]);constructor=AllChem.ETKDGv3;out=[]
for f in sorted((p/'Patchdock_Results').glob('docked_*.sdf'))[:5]:
 pos=docked_heads(f,hp,root=p).GetConformer().GetPositions();cut=heads[0].GetNumAtoms();center=pos.mean(0);targets=[pos[:cut]-center,pos[cut:]-center]
 for mode in ['basic_DG','DG_no_trans_constraint']:
  def variant():
   params=constructor();params.useExpTorsionAnglePrefs=False;params.useBasicKnowledge=False
   if mode=='DG_no_trans_constraint':params.forceTransAmides=False
   return params
  AllChem.ETKDGv3=variant;start=time.monotonic();mol,d=constrained_linker(full,heads,targets,155921,attempts=3)
  row={'pose':f.stem,'mode':mode,'accepted':mol is not None,'seconds':time.monotonic()-start,'diagnostics':d};out.append(row);print(json.dumps({k:v for k,v in row.items() if k!='diagnostics'}),flush=True)
  if mol is not None:
   xyz=mol.GetConformer().GetPositions()+center
   for i,point in enumerate(xyz):mol.GetConformer().SetAtomPosition(i,point)
   writer=Chem.SDWriter(str(r/'validation_controls'/('variant_'+f.stem+'_'+mode+'.sdf')));writer.write(mol);writer.close()
AllChem.ETKDGv3=constructor
(r/'validation_controls/sjf8240_DG_variant_control.json').write_text(json.dumps({'status':'EXPLORATORY_NOT_ADOPTED','poses':out},indent=2)+'\n')
