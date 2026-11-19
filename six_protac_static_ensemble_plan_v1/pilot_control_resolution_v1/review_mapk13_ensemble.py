"""Measure every saved linker conformer against the corresponding ternary protein."""
from pathlib import Path
import json
import sys

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

root=Path(sys.argv[1])
raw=root/"engine_raw/Patchdock_Results"
table=Chem.GetPeriodicTable()
full=Chem.MolFromSmiles((root/"full.smi").read_text().split()[0])
identity=Chem.MolToInchiKey(full)
rows=[]
for sdf in sorted(raw.glob("confs_*.sdf")):
    if not sdf.stat().st_size:
        continue
    suffix=sdf.stem.removeprefix("confs_")
    for state,protein in (("docked",raw/f"combined_{suffix}.pdb"),
                          ("repacked",raw/f"combined_{suffix}_0001.pdb")):
        if not protein.is_file():
            continue
        atoms=[line for line in protein.read_text().splitlines() if line.startswith("ATOM  ")
               and line[21] in "ABCD" and line[76:78].strip() not in {"H","D"}]
        px=np.array([[float(line[i:i+8]) for i in (30,38,46)] for line in atoms])
        pr=np.array([table.GetRvdw(x[76:78].strip().title()) for x in atoms])
        for index, conformer in enumerate(Chem.SDMolSupplier(str(sdf),removeHs=True),1):
            if conformer is None or Chem.MolToInchiKey(conformer)!=identity:
                raise ValueError(f"Invalid ligand identity: {sdf}:{index}")
            lr=np.array([table.GetRvdw(a.GetAtomicNum()) for a in conformer.GetAtoms()])
            dist=cdist(conformer.GetConformer().GetPositions(),px)
            overlap=lr[:,None]+pr-dist
            i,j=np.unravel_index(overlap.argmax(),overlap.shape)
            rows.append({"pose":suffix,"conformer":index,"protein_state":state,
                        "max_overlap_A":float(overlap[i,j]),"min_distance_A":float(dist.min()),
                        "worst_chain":atoms[j][21],"worst_residue":atoms[j][17:20]+atoms[j][22:27].strip()})
Path("MAPK13_SAVED_CONFORMERS.json").write_text(json.dumps(rows,indent=2)+"\n")
for state in ("docked","repacked"):
    print(state,sorted((r for r in rows if r["protein_state"]==state),key=lambda r:r["max_overlap_A"])[:3])
