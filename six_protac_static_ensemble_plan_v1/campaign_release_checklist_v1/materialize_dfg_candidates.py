"""Materialize and validate the 18 new DFG benchmark placements."""
import argparse
import hashlib
import json
from pathlib import Path
import runpy
import shutil
import sys

import numpy as np
from Bio import PDB
from rdkit import Chem

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("legacy",type=Path)
p.add_argument("project",type=Path)
p.add_argument("study",type=Path)
p.add_argument("output",type=Path)
args=p.parse_args()
source=args.legacy/"scripts/07_prepare_kinome_pRosettaC.py"
script=source.read_text()
boundary="# ── load DDR2 reference"
if script.count(boundary)!=1:
    raise ValueError("Unexpected legacy source")
helpers={}
exec(compile(script.split(boundary,1)[0],str(source),"exec"),helpers)
parser=PDB.PDBParser(QUIET=True)
sterics=runpy.run_path(str(args.study/"scripts/04_audit_prepared_coverage.py"))["sterics"]
audit=runpy.run_path(str(args.study/"ernie_candidate_manifest_v1/build_manifest.py"))["audit"]
registry=json.loads((args.project/"config/molecules.json").read_text())["compounds"]
benchmark=args.study/"sjf8240_mapk14_control_v1/DFG_SIX_REFERENCE_SUMMARY.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def window(path):
    residues=helpers["get_std_res"](parser.get_structure(path.stem,path))
    motif=helpers["find_dfg"](residues)
    if motif<0:
        raise ValueError(f"Missing DFG/DLG motif: {path}")
    return helpers["extract_window"](residues,motif)

args.output.mkdir(parents=True,exist_ok=False)
rows=[]
for candidate in json.loads(benchmark.read_text())["previously_missing_candidates"]:
    drug,uid=candidate["drug"],candidate["uid"]
    ref,=sorted((args.study/"ernie_six_control_pilot_v3/inputs"/drug).iterdir())
    full=Chem.MolFromSmiles((ref/"full.smi").read_text().split()[0])
    full_key=Chem.MolToInchiKey(full)
    if full_key!=registry[drug]["inchikey"]:
        raise ValueError(f"Supplier mismatch in reference {drug}")
    target=args.legacy/"data/structures"/f"{uid}.pdb"
    target_res,tidx=window(target)
    ref_res,ridx=window(ref/"protein.pdb")
    before=min(tidx,ridx)
    after=min(len(target_res)-tidx-1,len(ref_res)-ridx-1)
    molecule=Chem.SDMolSupplier(str(ref/"poi_head.sdf"),removeHs=True)[0]
    moved,rms=helpers["superimpose_and_transform"](
        helpers["ca_array"](ref_res)[ridx-before:ridx+after+1],
        helpers["ca_array"](target_res)[tidx-before:tidx+after+1],
        molecule.GetConformer().GetPositions())
    if abs(rms-candidate["fit_rmsd_A"])>1e-6:
        raise ValueError(f"DFG fit changed from benchmark {drug}/{uid}")
    name=f"{candidate['gene']}_{uid}__legacy_DFG_AF2"
    out=args.output/"inputs"/drug/name
    out.mkdir(parents=True)
    for f in ("vbc.pdb","vhl_head.sdf","full.smi","Protac_params.txt","selection.json"):
        shutil.copy2(ref/f,out/f)
    helpers["write_pdb"](out/"protein.pdb",target_res)
    protein=out/"protein.pdb"
    protein.write_text("".join(line[:21]+"D"+line[22:] if line.startswith("ATOM  ") else line
                               for line in protein.read_text().splitlines(keepends=True)))
    for i,coords in enumerate(moved):
        molecule.GetConformer().SetAtomPosition(i,coords)
    with Chem.SDWriter(str(out/"poi_head.sdf")) as writer:
        writer.write(molecule)
    params=dict((k.strip(),v.strip()) for k,v in
                (l.split(":",1) for l in (out/"Protac_params.txt").read_text().splitlines()))
    links=[]
    for head_name,anchor in zip(params["Heads"].split(),map(int,params["Anchor atoms"].split())):
        head=Chem.SDMolSupplier(str(out/head_name),removeHs=True)[0]
        mappings=full.GetSubstructMatches(head,useChirality=True,uniquify=False)
        sites={(m[anchor-1]+1,tuple(n.GetIdx()+1 for n in full.GetAtomWithIdx(m[anchor-1]).GetNeighbors()
                                  if n.GetIdx() not in m)) for m in mappings}
        if len(sites)!=1 or len(next(iter(sites))[1])!=1:
            raise ValueError(f"Nonunique attachment {drug}/{uid}/{head_name}")
        links.append({"head":head_name,"one_based_anchor":anchor,
                      "full_ligand_attachment":next(iter(sites))})
    geometry=sterics(out/"poi_head.sdf",protein)
    if not 1.5<=geometry["min_distance_A"]<=4.5 or geometry["max_vdw_overlap_A"]>.6:
        raise ValueError(f"Serialized placement fails benchmark geometry: {drug}/{uid}")
    info={"status":"CHEMISTRY_AND_STERIC_CHECKED_CANDIDATE","production_ready":False,
          "method":"legacy_DFG_offset_transfer","compound":drug,"uniprot":uid,"gene":candidate["gene"],
          "source_sha256":sha(target),"legacy_code_sha256":sha(source),"reference_head_sha256":sha(ref/"poi_head.sdf"),
          "benchmark_sha256":sha(benchmark),"geometry":geometry,"dfg_fit_rmsd_A":rms,
          "anchors":links,"files":{f.name:sha(f) for f in out.iterdir()}}
    (out/"preparation.json").write_text(json.dumps(info,indent=2)+"\n")
    hashes=audit(out,full_key)
    rows.append({**info,"input":str(out.relative_to(args.output)),"input_sha256":hashes})
(args.output/"DFG_CANDIDATE_RECEIPTS.json").write_text(json.dumps(rows,indent=2)+"\n")
print(len(rows))
