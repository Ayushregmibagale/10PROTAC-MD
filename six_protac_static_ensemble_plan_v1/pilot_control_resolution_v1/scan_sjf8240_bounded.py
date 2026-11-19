"""Compare archived SJF8240/MET geometry with a chemically consistent synthetic control."""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem

pilot=Path(sys.argv[1])
adapter=pilot.parent/"v3_conformer_repair_v1/runtime_bundle/baseline_engine/runtime/adapter"
sys.path.insert(0,str(adapter))
from constrained_linker import constrained_linker, match_heads, stereo_key
from official_parameters import docked_heads

raw=next((pilot/"raw_results/attempts").glob("SJF8240__*/engine_raw"))
head_paths=[raw/"vhl_head_H.sdf",raw/"poi_head_H.sdf"]
heads=[Chem.SDMolSupplier(str(p),removeHs=True)[0] for p in head_paths]
full=Chem.MolFromSmiles((raw/"full.smi").read_text().strip())
match=match_heads(full,heads)[0]
good=None
for candidate in Chem.SDMolSupplier(str(raw/"random_sampling.sdf"),removeHs=True):
    if candidate is not None and stereo_key(candidate)==Chem.MolToInchiKey(full):
        good=candidate
        break
if good is None:
    raise ValueError("No archived verified conformer for the synthetic geometry control")
canonical=good.GetSubstructMatch(full,useChirality=True)
if len(canonical)!=full.GetNumAtoms():
    raise ValueError("Unable to map synthetic control conformer to the registered graph")
synthetic=good.GetConformer().GetPositions()[list(canonical)]
sets={"synthetic_full_molecule_control":[synthetic[list(indices)] for indices in match]}
for index in range(1,51):
    path=raw/"Patchdock_Results"/f"pd.{index}.pdb"
    xyz=docked_heads(str(path),head_paths,raw).GetConformer().GetPositions()
    cut=heads[0].GetNumAtoms()
    sets[path.stem]=[xyz[:cut],xyz[cut:]]
def run_case(item):
    label, coordinates = item
    compound, stats = constrained_linker(full, heads, coordinates, 42, attempts=20, embed_iterations=20)
    return {"geometry": label, "feasible": compound is not None, "details": stats}

path=Path("SJF8240_GLOBAL_CLOSURE_BOUNDED.jsonl")
with path.open("x", buffering=1) as ledger, ProcessPoolExecutor(max_workers=4) as workers:
    for row in workers.map(run_case,sets.items()):
        ledger.write(json.dumps(row)+"\n")
        print(row["geometry"],row["feasible"],flush=True)
