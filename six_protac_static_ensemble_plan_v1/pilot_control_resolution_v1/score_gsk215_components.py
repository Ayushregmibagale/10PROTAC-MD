"""Decompose archived ref2015 scores into component and interaction diagnostics."""
import json
import os
from pathlib import Path
import sys

import pyrosetta

raw=Path(sys.argv[1]).resolve()
pyrosetta.init("-ignore_unrecognized_res false -load_PDB_components false -out:level 100")
score=pyrosetta.create_score_function("ref2015")
root=Path.cwd()
results=[]
flags={"VBC":"ABC","PTK2":"D","PROTAC":"X","protein_pair":"ABCD","all":"ABCDX"}
from pyrosetta import generate_nonstandard_residue_set

def component_energy(records, flags, params):
    selected=[line for line in records if line[21] in flags]
    block=[]
    for line in selected:
        if block and block[-1].startswith(("ATOM  ","HETATM")) and block[-1][21]!=line[21]:
            block.append("TER")
        block.append(line)
    file=root/"component_input.pdb"
    file.write_text("\n".join(block)+"\nTER\nEND\n")
    pose=pyrosetta.Pose()
    if "X" in flags:
        generate_nonstandard_residue_set(pose,[str(params)])
    pyrosetta.pose_from_file(pose,str(file))
    if "X" in flags and sum(r.name3()=="PTC" for r in pose.residues)!=1:
        raise ValueError("PROTAC did not load uniquely")
    return float(score(pose))

for path in sorted(raw.glob("combined_*_0001.pdb")):
    pose=path.stem
    parts=pose.split("_")
    params=raw/f"PT_{parts[1]}_{parts[2]}.params"
    records=[line for line in path.read_text().splitlines() if line.startswith(("ATOM  ","HETATM"))]
    values={name:component_energy(records,ids,params) for name,ids in flags.items()}
    values.update({"pose":pose,
                   "protein_pair_delta":values["protein_pair"]-values["VBC"]-values["PTK2"],
                   "protac_protein_delta":values["all"]-values["protein_pair"]-values["PROTAC"]})
    results.append(values)
(root/"COMPONENT_SCORES.json").write_text(json.dumps(results,indent=2)+"\n")
print(len(results))
