from pathlib import Path
import json,sys,time
import numpy as np
from rdkit import Chem
r=Path('/tmp/10protac_finish/prosettac_robust');w=r/'workflow';sys.path[:0]=[str(w/'baseline_engine'),str(w/'scripts/lib')]
from constrained_linker import constrained_linker
from official_parameters import docked_heads
p=r/'validation_controls/sjf8240_global50_diagnostic_v2/engine_raw';heads=[Chem.SDMolSupplier(str(p/n),removeHs=True)[0] for n in ['vhl_head_H.sdf','poi_head_H.sdf']];full=Chem.MolFromSmiles((p/'full.smi').read_text().split()[0]);out=[]
for f in sorted((p/'Patchdock_Results').glob('docked_*.sdf'))[:5]:
 ref=docked_heads(f,[p/'vhl_head_H.sdf',p/'poi_head_H.sdf'],root=p);pos=ref.GetConformer().GetPositions();cut=heads[0].GetNumAtoms();targets=[pos[:cut],pos[cut:]];center=pos.mean(0)
 for label,shift in [('absolute',np.zeros(3)),('centered',center)]:
  start=time.monotonic();mol,d=constrained_linker(full,heads,[a-shift for a in targets],155921,attempts=3)
  row={'pose':f.stem,'mode':label,'center':center.tolist(),'accepted':mol is not None,'seconds':time.monotonic()-start,'diagnostics':d};out.append(row);print(json.dumps(row),flush=True)
(r/'validation_controls/sjf8240_translation_control.json').write_text(json.dumps(out,indent=2)+'\n')
