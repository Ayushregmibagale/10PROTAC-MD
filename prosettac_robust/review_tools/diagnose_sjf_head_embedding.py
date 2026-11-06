from pathlib import Path
import sys,json,time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
r=Path(__file__).resolve().parents[1];w=r/'workflow';sys.path[:0]=[str(w/'baseline_engine'),str(w/'scripts/lib')]
from constrained_linker import match_heads
from official_parameters import docked_heads
p=r/'validation_controls/sjf8240_global50_diagnostic_v2/engine_raw';heads=[Chem.SDMolSupplier(str(p/n),removeHs=True)[0] for n in ['vhl_head_H.sdf','poi_head_H.sdf']];full=Chem.MolFromSmiles((p/'full.smi').read_text().split()[0]);maps=match_heads(full,heads)[0];docked=p/'Patchdock_Results/docked_1_1.sdf';reference=docked_heads(docked,[p/'vhl_head_H.sdf',p/'poi_head_H.sdf'],root=p);pos=reference.GetConformer().GetPositions();cut=heads[0].GetNumAtoms();targets=[pos[:cut],pos[cut:]];out=[]
for label,which in [('free',[]),('VHL_only',[0]),('POI_only',[1]),('both',[0,1])]:
 fixed={i:xyz for k in which for i,xyz in zip(maps[k],targets[k])};success=[]
 for trial in range(3):
  m=Chem.AddHs(Chem.Mol(full));params=AllChem.ETKDGv3();params.randomSeed=155921+trial;params.maxIterations=20;params.numThreads=1;params.useRandomCoords=True;params.trackFailures=True;params.SetCoordMap({i:Point3D(*map(float,p)) for i,p in fixed.items()});start=time.monotonic();cid=AllChem.EmbedMolecule(m,params);success.append({'seed':params.randomSeed,'cid':cid,'failure_counts':list(params.GetFailureCounts()),'seconds':time.monotonic()-start})
 row={'mode':label,'fixed_atoms':len(fixed),'attempts':success};out.append(row);print(json.dumps(row),flush=True)
(r/'validation_controls/sjf8240_head_embedding_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
