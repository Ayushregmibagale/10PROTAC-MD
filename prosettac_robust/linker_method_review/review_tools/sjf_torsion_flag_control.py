from pathlib import Path
import sys,json,time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
r=Path('/tmp/10protac_finish/prosettac_robust');w=r/'workflow';sys.path[:0]=[str(w/'baseline_engine'),str(w/'scripts/lib')];from constrained_linker import constrained_linker;from official_parameters import docked_heads
p=r/'validation_controls/sjf8240_global50_diagnostic_v2/engine_raw';hp=[p/'vhl_head_H.sdf',p/'poi_head_H.sdf'];heads=[Chem.SDMolSupplier(str(x),removeHs=True)[0] for x in hp];full=Chem.MolFromSmiles((p/'full.smi').read_text().split()[0]);pos=docked_heads(p/'Patchdock_Results/docked_10_1.sdf',hp,root=p).GetConformer().GetPositions();pos-=pos.mean(0);targets=[pos[:30],pos[30:]];constructor=AllChem.ETKDGv3;rows=[]
for exp,basic in [(False,True),(True,False)]:
 def params():
  o=constructor();o.useExpTorsionAnglePrefs=exp;o.useBasicKnowledge=basic;return o
 AllChem.ETKDGv3=params;start=time.monotonic();mol,d=constrained_linker(full,heads,targets,155921,attempts=1);row={'experimental_torsion_preferences':exp,'basic_knowledge':basic,'accepted':mol is not None,'seconds':time.monotonic()-start,'details':d};rows.append(row);print(json.dumps(row),flush=True)
AllChem.ETKDGv3=constructor
(r/'validation_controls/sjf8240_torsion_flag_control.json').write_text(json.dumps({'status':'EXPLORATORY_NOT_ADOPTED','poses':rows},indent=2)+'\n')
