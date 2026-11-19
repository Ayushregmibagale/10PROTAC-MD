"""Locally minimize archived MAPK13 linkers with explicit fixed-protein sterics."""
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.optimize import minimize
from scipy.spatial.distance import cdist

root, adapter = map(Path, sys.argv[1:])
sys.path.insert(0,str(adapter))
from constrained_linker import match_heads, stereo_key

raw=root/"engine_raw"
out=Path("linker_relaxation")
out.mkdir()
full=Chem.MolFromSmiles((raw/"full.smi").read_text().split()[0])
heads=[Chem.SDMolSupplier(str(raw/f),removeHs=True)[0] for f in ("vhl_head_H.sdf","poi_head_H.sdf")]
maps=match_heads(full,heads)
ptable=Chem.GetPeriodicTable()
rows=[]
for sdf in sorted((raw/"Patchdock_Results").glob("confs_*.sdf")):
    if not sdf.stat().st_size:
        continue
    suffix=sdf.stem.removeprefix("confs_")
    sampling=json.loads(Path(f"{sdf}.sampling.json").read_text())
    accepted=[x for x in sampling["diagnostics"] if x["accepted"]]
    protein=raw/"Patchdock_Results"/f"combined_{suffix}_0001.pdb"
    lines=[l for l in protein.read_text().splitlines() if l.startswith("ATOM  ")
           and l[21] in "ABCD" and l[76:78].strip() not in {"H","D"}]
    protein_xyz=np.array([[float(l[i:i+8]) for i in (30,38,46)] for l in lines])
    protein_r=np.array([ptable.GetRvdw(l[76:78].strip().title()) for l in lines])
    molecules=list(Chem.SDMolSupplier(str(sdf),removeHs=True))
    if len(molecules)!=len(accepted):
        raise ValueError("Conformer records and fixed-head receipts differ")
    for index,(m,receipt) in enumerate(zip(molecules,accepted),1):
        if m is None or stereo_key(m)!=Chem.MolToInchiKey(full):
            raise ValueError("Starting ligand graph or stereochemistry error")
        if [a.GetAtomicNum() for a in m.GetAtoms()] != [a.GetAtomicNum() for a in full.GetAtoms()] or m.GetSubstructMatch(full,useChirality=True)!=tuple(range(full.GetNumAtoms())):
            raise ValueError("Starting ligand atom order differs from registered input")
        fixed=np.array(sum(receipt["head_mappings"],[]))
        if tuple(tuple(x) for x in receipt["head_mappings"]) not in maps:
            raise ValueError("Invalid fixed-head mapping in receipt")
        hmol=Chem.AddHs(m,addCoords=True)
        start=hmol.GetConformer().GetPositions()
        free=np.ones(len(start),dtype=bool)
        free[fixed]=False
        n=m.GetNumAtoms()
        ligand_r=np.array([ptable.GetRvdw(a.GetAtomicNum()) for a in m.GetAtoms()])
        neighbors=cdist(start[:n],protein_xyz).min(axis=0)<6
        px=protein_xyz[neighbors]
        pr=protein_r[neighbors]
        field=AllChem.UFFGetMoleculeForceField(hmol)

        def objective(x):
            coords=start.copy()
            coords[free]=x.reshape(-1,3)
            flat=tuple(coords.ravel())
            value=field.CalcEnergy(flat)
            grad=np.asarray(field.CalcGrad(flat)).reshape(-1,3)
            delta=coords[:n,None,:]-px[None,:,:]
            dist=np.maximum(np.linalg.norm(delta,axis=2),1e-8)
            excess=np.maximum(ligand_r[:,None]+pr-dist-.35,0)
            value+=20*np.sum(excess**2)
            grad[:n]-=40*np.sum((excess/dist)[:,:,None]*delta,axis=1)
            return float(value),grad[free].ravel()

        opt=minimize(objective,start[free].ravel(),jac=True,method="L-BFGS-B",options={"maxiter":700,"ftol":1e-10})
        end=start.copy()
        end[free]=opt.x.reshape(-1,3)
        for i,coord in enumerate(end):
            hmol.GetConformer().SetAtomPosition(i,coord)
        refined=Chem.RemoveHs(hmol)
        valid=stereo_key(refined)==Chem.MolToInchiKey(full) and all(
            .9<np.linalg.norm(end[b.GetBeginAtomIdx()]-end[b.GetEndAtomIdx()])<2.1 for b in refined.GetBonds())
        def overlap(x):
            return float((ligand_r[:,None]+protein_r-cdist(x[:n],protein_xyz)).max())
        rows.append({"pose":suffix,"conformer":index,"converged":bool(opt.success),
                     "identity_and_bonds_valid":bool(valid),"before_A":overlap(start),"after_A":overlap(end),
                     "fixed_head_max_displacement_A":float(np.linalg.norm(end[fixed]-start[fixed],axis=1).max())})
        if opt.success and valid:
            with Chem.SDWriter(str(out/f"{suffix}_{index}.sdf")) as writer:
                writer.write(refined)
(out/"RESULTS.json").write_text(json.dumps(rows,indent=2)+"\n")
print(sorted(rows,key=lambda x:x["after_A"])[:4])
