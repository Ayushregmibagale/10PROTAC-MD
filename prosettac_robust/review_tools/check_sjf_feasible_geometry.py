from pathlib import Path
import json,sys,time
import numpy as np
from rdkit import Chem
r=Path(__file__).resolve().parents[1];w=r/'workflow';sys.path.insert(0,str(w/'baseline_engine'));from constrained_linker import match_heads,constrained_linker
p=r/'validation_controls/sjf8240_global50_diagnostic_v2/engine_raw';full=Chem.MolFromSmiles((p/'full.smi').read_text().split()[0]);heads=[Chem.SDMolSupplier(str(p/n),removeHs=True)[0] for n in ['vhl_head_H.sdf','poi_head_H.sdf']];reference=next(m for m in Chem.SDMolSupplier(str(p/'random_sampling.sdf'),removeHs=True) if m is not None);mapping=match_heads(full,heads)[0];targets=[reference.GetConformer().GetPositions()[list(idx)] for idx in mapping];start=time.monotonic();mol,d=constrained_linker(full,heads,targets,155921);result={'status':'KNOWN_FEASIBLE_SJF8240_SELF_GEOMETRY_CONTROL','accepted':mol is not None,'seconds':time.monotonic()-start,'diagnostic':d,'limitations':['Head coordinates are taken from a feasible free full-molecule conformer, not experimental bound heads.','No protein complex, binding or native recovery claim.']}
if mol is not None:
 writer=Chem.SDWriter(str(r/'validation_controls/sjf8240_self_geometry_control.sdf'));writer.write(mol);writer.close()
(r/'validation_controls/sjf8240_self_geometry_control.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
