"""Replay fixed-linker closure with a small rigid refit of the VHL-bound head."""
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

raw=Path("broader_orientation_attempt/engine_raw").resolve()
sys.path.insert(0,sys.argv[1])
from constrained_linker import constrained_linker,fit,stereo_key
from official_parameters import docked_heads

heads_paths=[raw/"vhl_head_H.sdf",raw/"poi_head_H.sdf"]
heads=[Chem.SDMolSupplier(str(p),removeHs=True)[0] for p in heads_paths]
refit=Chem.SDMolSupplier("vhl_refit_result/poi_head.sdf",removeHs=True)[0]
full=Chem.MolFromSmiles((raw/"full.smi").read_text().split()[0])
identity=Chem.MolToInchiKey(full)
original=Chem.SDMolSupplier("mapk13_vhl_binary/poi_head.sdf",removeHs=True)[0]
if Chem.MolToSmiles(original)!=Chem.MolToSmiles(refit):
    raise ValueError("Refit altered the VHL head graph")
ptable=Chem.GetPeriodicTable()
lr=np.array([ptable.GetRvdw(a.GetAtomicNum()) for a in full.GetAtoms()])
out=Path("vhl_refit_linkers")
out.mkdir()
results=[]
for suffix in ("19_2","5_2"):
    source=raw/"Patchdock_Results"
    docked=docked_heads(str(source/f"docked_{suffix}.pdb"),heads_paths,raw)
    pos=docked.GetConformer().GetPositions()
    count=original.GetNumAtoms()
    r,t=fit(original.GetConformer().GetPositions(),pos[:count])
    discrepancy=np.max(np.linalg.norm(original.GetConformer().GetPositions()@r+t-pos[:count],axis=1))
    if discrepancy>.01:
        raise ValueError(f"Cannot map original VHL head to saved pose: {discrepancy}")
    target=[refit.GetConformer().GetPositions()@r+t,pos[count:]]
    protein=source/f"combined_{suffix}_0001.pdb"
    atoms=[line for line in protein.read_text().splitlines() if line.startswith("ATOM  ")
           and line[21] in "ABCD" and line[76:78].strip() not in {"H","D"}]
    px=np.array([[float(line[i:i+8]) for i in (30,38,46)] for line in atoms])
    pr=np.array([ptable.GetRvdw(line[76:78].strip().title()) for line in atoms])
    for seed in (42,7,123):
        mol,detail=constrained_linker(full,heads,target,seed)
        row={"pose":suffix,"seed":seed,"feasible":mol is not None}
        if mol is not None:
            if stereo_key(mol)!=identity:
                raise ValueError("Incorrect linker identity")
            distances=cdist(mol.GetConformer().GetPositions(),px)
            vdw=lr[:,None]+pr-distances
            _,idx=np.unravel_index(vdw.argmax(),vdw.shape)
            row.update(max_overlap_A=float(vdw.max()),min_distance_A=float(distances.min()),
                       worst_chain=atoms[idx][21], worst_residue=atoms[idx][17:20]+atoms[idx][22:27].strip())
            with Chem.SDWriter(str(out/f"{suffix}_{seed}.sdf")) as w:
                w.write(mol)
        results.append(row)
(out/"RESULTS.json").write_text(json.dumps(results,indent=2)+"\n")
print(results)
