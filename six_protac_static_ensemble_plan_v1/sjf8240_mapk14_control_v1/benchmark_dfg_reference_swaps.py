"""Measure the historical DFG placement routine with all six frozen head references."""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys
from Bio import PDB

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

legacy, project, study = map(Path, sys.argv[1:])
script = legacy / "scripts/07_prepare_kinome_pRosettaC.py"
source = script.read_text()
boundary = "# ── load DDR2 reference"
if source.count(boundary) != 1:
    raise ValueError("Cannot isolate the legacy helper definitions")
helpers = {}
exec(compile(source.split(boundary,1)[0],str(script),"exec"),helpers)
parser = PDB.PDBParser(QUIET=True)

def structure_window(path):
    residues = helpers["get_std_res"](parser.get_structure(path.stem,path))
    motif = helpers["find_dfg"](residues)
    if motif < 0 or len(residues)<80:
        return None
    window,offset = helpers["extract_window"](residues,motif)
    return (window,offset,helpers["ca_array"](window)) if len(window)>=40 else None

refs={}
table=Chem.GetPeriodicTable()
for folder in sorted((study/"ernie_six_control_pilot_v3/inputs").glob("*/*")):
    head = Chem.SDMolSupplier(str(folder/"poi_head.sdf"),removeHs=True)[0]
    full = Chem.MolFromSmiles((folder/"full.smi").read_text().split()[0])
    window = structure_window(folder/"protein.pdb")
    if head is None or not full.HasSubstructMatch(head,useChirality=True) or window is None:
        raise ValueError(f"Reference unavailable: {folder}")
    ligand_r = np.array([table.GetRvdw(a.GetAtomicNum()) for a in head.GetAtoms()])
    refs[folder.parent.name]=(window,head.GetConformer().GetPositions(),ligand_r)

with (project/"structure_provenance/domains/domain_inventory.tsv").open() as f:
    genes={row["uniprot_id"]:row["gene"] for row in csv.DictReader(f,delimiter="\t")}
rows=[]
for uid,gene in sorted(genes.items()):
    path=legacy/"data/structures"/f"{uid}.pdb"
    window=structure_window(path) if path.is_file() else None
    if window is None:
        reason="missing_legacy_AF2_file" if not path.is_file() else "no_eligible_DFG_window"
        for drug in refs:
            rows.append({"drug":drug,"uid":uid,"gene":gene,"status":reason})
        continue
    residues,body_dfg,body_ca=window
    protein_atoms=[a for r in residues for a in r if a.element not in {"H","D"}]
    px=np.array([a.coord for a in protein_atoms])
    radius=np.array([table.GetRvdw(a.element.title()) for a in protein_atoms])
    tree=cKDTree(px)
    for drug, ((_,ref_dfg,ref_ca),xyz,head_radius) in refs.items():
        before=min(ref_dfg,body_dfg)
        after=min(len(ref_ca)-ref_dfg-1,len(body_ca)-body_dfg-1)
        if before+after+1<20:
            rows.append({"drug":drug,"uid":uid,"gene":gene,"status":"short_overlap"})
            continue
        fitted,rms=helpers["superimpose_and_transform"](
             ref_ca[ref_dfg-before:ref_dfg+after+1],body_ca[body_dfg-before:body_dfg+after+1],xyz)
        distance=float(tree.query(fitted)[0].min())
        clashes=[float(head_radius[i]+radius[j]-np.linalg.norm(fitted[i]-px[j]))
                 for i,neighbors in enumerate(tree.query_ball_point(fitted,4.5)) for j in neighbors]
        rows.append({"drug":drug,"uid":uid,"gene":gene,"status":"placement_generated",
                     "fit_rmsd_A":rms,"min_distance_A":distance,
                     "max_vdw_overlap_A":max(clashes,default=None),
                     "passes_minimum_contact_screen":1.5<=distance<=4.5})
counts={}
for drug in refs:
    selected=[r for r in rows if r["drug"]==drug]
    counts[drug]={"counts":dict(Counter(r["status"] for r in selected)),
                  "contact_screen_pass":sum(r.get("passes_minimum_contact_screen",False) for r in selected),
                  "overlap_at_most_0p6_A":sum(r.get("max_vdw_overlap_A") is not None and
                       r["max_vdw_overlap_A"]<=.6 and r.get("passes_minimum_contact_screen",False) for r in selected)}
Path("DFG_SIX_REFERENCE_BENCHMARK.json").write_text(json.dumps({"legacy_code_sha256":hashlib.sha256(script.read_bytes()).hexdigest(),
                       "scope":"Geometry-only diagnostic; no input release","counts":counts,"rows":rows},indent=2)+"\n")
print(json.dumps(counts,indent=2))
