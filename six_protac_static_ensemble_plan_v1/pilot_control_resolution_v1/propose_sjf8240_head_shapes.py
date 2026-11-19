"""Test ligand-driven SJF8240 orientation proposals using archived free conformers."""
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

pilot=Path(sys.argv[1])
adapter=pilot.parent/"v3_conformer_repair_v1/runtime_bundle/baseline_engine/runtime/adapter"
sys.path.insert(0,str(adapter))
from constrained_linker import fit,match_heads,stereo_key

raw=next((pilot/"raw_results/attempts").glob("SJF8240__*/engine_raw"))
paths=[raw/"vhl_head_H.sdf",raw/"poi_head_H.sdf"]
heads=[Chem.SDMolSupplier(str(p),removeHs=True)[0] for p in paths]
full=Chem.MolFromSmiles((raw/"full.smi").read_text().split()[0])
identity=Chem.MolToInchiKey(full)
mappings=match_heads(full,heads)
coords=[head.GetConformer().GetPositions() for head in heads]
proteins=[[line for line in (raw/name).read_text().splitlines()
   if line.startswith("ATOM  ") and line[76:78].strip() not in {"H","D"}] for name in ("Init0.pdb","Init1.pdb")]
px=[np.array([[float(line[i:i+8]) for i in (30,38,46)] for line in letters]) for letters in proteins]
tree=cKDTree(px[0])
ptable=Chem.GetPeriodicTable()
pr=[np.array([ptable.GetRvdw(l[76:78].strip().title()) for l in lines]) for lines in proteins]
lr=np.array([ptable.GetRvdw(a.GetAtomicNum()) for a in full.GetAtoms()])
# Constrain the internal geometry of each head independently. Their relative
# rigid-body arrangement remains free during the conformer proposal.
proposal=Chem.AddHs(full)
bounds=rdDistGeom.GetMoleculeBoundsMatrix(proposal)
a,b=mappings[0]
for group,template in zip((a,b),coords):
    distances=cdist(template,template)
    for i,atom in enumerate(group):
        for j,other in enumerate(group[:i]):
            low,high=sorted((atom,other))
            distance=distances[i,j]
            bounds[high,low]=max(0,distance-.10)
            bounds[low,high]=distance+.10
parameters=AllChem.ETKDGv3()
parameters.randomSeed=42
parameters.useRandomCoords=True
parameters.numThreads=2
parameters.maxIterations=200
parameters.trackFailures=True
parameters.SetBoundsMat(bounds)
ids=AllChem.EmbedMultipleConfs(proposal,numConfs=40,params=parameters)
results=[]
for number,cid in enumerate(ids,1):
    copy=Chem.Mol(proposal)
    conf=Chem.Conformer(proposal.GetConformer(cid))
    copy.RemoveAllConformers()
    copy.AddConformer(conf)
    mol=Chem.RemoveHs(copy)
    if mol is None or stereo_key(mol)!=identity:
        continue
    canonical=mol.GetSubstructMatch(full,useChirality=True)
    if len(canonical)!=full.GetNumAtoms():
        raise ValueError("Cannot map an archived conformer to registered connectivity")
    ligand=mol.GetConformer().GetPositions()[list(canonical)]
    for a,b in mappings:
        first,t1=fit(ligand[list(a)],coords[0])
        placed=ligand@first+t1
        second,t2=fit(coords[1],placed[list(b)])
        errs=[np.linalg.norm(placed[list(a)]-coords[0],axis=1),
              np.linalg.norm(coords[1]@second+t2-placed[list(b)],axis=1)]
        rms=[float(np.sqrt(np.mean(e**2))) for e in errs]
        if max(rms)>.35 or max(e.max() for e in errs)>1:
            results.append({"conformer":number,"head_rmsd_A":rms,"status":"HEAD_SHAPE_REJECTED"})
            continue
        moved=px[1]@second+t2
        d0=cdist(placed,px[0])
        d1=cdist(placed,moved)
        e3_to_target_min=float(tree.query(moved)[0].min())
        results.append({"conformer":number,"head_rmsd_A":rms,"min_e3_to_target_A":e3_to_target_min,
                        "min_ligand_protein_A":float(min(d0.min(),d1.min())),
                        "max_ligand_protein_vdw_overlap_A":float(max((lr[:,None]+pr[0]-d0).max(),
                                                                    (lr[:,None]+pr[1]-d1).max())),
                        "rotation":second.tolist(),"translation":t2.tolist(),"status":"GEOMETRY_EVALUATED"})
Path("SJF8240_HEAD_SHAPE_PROPOSALS.json").write_text(json.dumps(results,indent=2)+"\n")
supported=[r for r in results if r["status"]=="GEOMETRY_EVALUATED"]
separated=[r for r in supported if min(r["min_e3_to_target_A"],r["min_ligand_protein_A"])>=1.5]
print(len(ids),len(results),len(supported),len(separated),tuple(parameters.GetFailureCounts()))
print(sorted(separated,key=lambda x:x["max_ligand_protein_vdw_overlap_A"])[:2])
