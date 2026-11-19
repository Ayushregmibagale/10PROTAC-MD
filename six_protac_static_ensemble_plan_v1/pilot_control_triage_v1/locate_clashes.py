"""Locate major heavy-atom overlaps in retained GSK215 and MAPK13 poses."""
import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist


def atoms(path, condition):
    return [line for line in path.read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM")) and condition(line)
            and line[76:78].strip().upper() not in {"H", "D"}]


def xyz(lines):
    return np.array([[float(line[i:i+8]) for i in (30,38,46)] for line in lines])


def pose_hits(final, reference):
    ligand = atoms(final, lambda line: line[17:20] == "PTC" and line[21] == "X")
    proteins = atoms(final, lambda line: line.startswith("ATOM  ") and line[21] in "ABCD")
    heads = {key: xyz(atoms(reference, lambda line: line[17:20] == key)) for key in ("PT0", "PT1")}
    ligand_xyz, protein_xyz = xyz(ligand), xyz(proteins)
    distances = cdist(ligand_xyz, protein_xyz)
    head_distances = {key: cdist(ligand_xyz, coords).min(axis=1) for key, coords in heads.items()}
    pt = Chem.GetPeriodicTable()
    def radii(lines):
        return np.array([pt.GetRvdw(line[76:78].strip().title()) for line in lines])
    overlap = radii(ligand)[:,None] + radii(proteins) - distances
    rows=[]
    for index in np.argsort(overlap.ravel())[-10:][::-1]:
        a,b=np.unravel_index(index,overlap.shape)
        if overlap[a,b] <= 0:
            continue
        head=min(head_distances,key=lambda key:head_distances[key][a])
        prot=proteins[b]
        rows.append({"ligand_atom":ligand[a][12:16].strip(),
                     "ligand_region":head if head_distances[head][a] <= .05 else "outside_fixed_heads",
                     "protein_chain":prot[21],"residue":prot[17:20]+prot[22:27].strip(),
                     "protein_atom":prot[12:16].strip(),"distance_A":float(distances[a,b]),
                     "overlap_A":float(overlap[a,b])})
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    cases={
      "GSK215":next((args.project/"ernie_six_control_pilot_v3/raw_results/attempts").glob("GSK215__*")),
      "SJF_delta_MAPK13":args.project/"pilot_control_review_v1/pilot_attempt"
    }
    review={}
    for name,attempt in cases.items():
        rows={}
        for pose in sorted((attempt/"Patchdock_Results").glob("combined_*.pdb")):
            parts=pose.stem.split("_")
            docked=attempt/"engine_raw/Patchdock_Results"/f"docked_{parts[1]}_{parts[2]}.pdb"
            rows[pose.stem]=pose_hits(pose,docked)
        review[name]=rows
    args.output.write_text(json.dumps(review,indent=2)+"\n")
    for pose,hits in review["SJF_delta_MAPK13"].items():
        print(pose,hits[0])


if __name__=="__main__":
    main()
