"""Compare the historical DFG input with an explicit 5IA4 to MAPK14 transfer."""
import importlib.util
import json
from pathlib import Path
import runpy
import shutil
import sys

import numpy as np
from Bio.PDB import PDBParser
from rdkit import Chem

project=Path(sys.argv[1])
study=Path(sys.argv[2])
legacy=Path.home()/"kinome_scan/data/kinase_pRosettaC/Q16539"
scripts=project/"scripts/01_setup"
sys.path.insert(0,str(project/"scripts/lib"))
import protac_qc as qc
def load(name):
    spec=importlib.util.spec_from_file_location(name,scripts/name)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
sterics=runpy.run_path(str(study/"scripts/04_audit_prepared_coverage.py"))["sterics"]
delta=load("00j_prepare_delta_p38_control.py")
bound=load("00d_prepare_bound_references.py")
frag=load("00f_prepare_foretinib_heads.py")
placement=load("01_place_warhead_batch.py")
references=project/"reference_sources"
structure=PDBParser(QUIET=True).get_structure("EPHA2",references/"5IA4.pdb")
parent,_,_=bound.bound_molecule(structure,references/"88Z.cif","88Z","A")
matches=parent.GetSubstructMatches(Chem.MolFromSmiles(frag.CORE),useChirality=True)
if len(matches)!=1:
    raise ValueError("Foretinib template head did not map uniquely")
head,_=frag.fragment(parent,matches[0])
base=study/"ernie_six_control_pilot_v3/inputs/SJF8240/MET_P08581"
full=Chem.MolFromSmiles((base/"full.smi").read_text().strip())
mapped=full.GetSubstructMatches(head,useChirality=True,uniquify=False)
attachments={i:{m[i] for m in mapped if any(n.GetIdx() not in m for n in full.GetAtomWithIdx(m[i]).GetNeighbors())}
             for i in range(head.GetNumAtoms())}
exits=[i for i,hits in attachments.items() if hits]
if len(exits)!=1 or len(attachments[exits[0]])!=1:
    raise ValueError("Ambiguous supplier-specific exit atom")
targets=[a for a in qc.pdb_atoms("1W82.pdb") if a["record"]=="ATOM" and a["chain"]=="A" and 24<=a["resnum"]<=308]
if not targets:
    raise ValueError("MAPK14 target template empty")
def ca(path,chain):
    return {str(a["resnum"])+a["icode"].strip():a["xyz"] for a in qc.pdb_atoms(path)
            if a["record"]=="ATOM" and a["chain"]==chain and a["name"]=="CA"}
keys,x,y=delta.pocket_pairs(json.loads((references/"KLIFS_7189_pocket_residue_map.json").read_text()),
        json.loads(Path("1W82_pocket.json").read_text()),ca(references/"5IA4.pdb","A"),ca("1W82.pdb","A"))
xyz=head.GetConformer().GetPositions()
mask=np.linalg.norm(x[:,None,:]-xyz[None,:,:],axis=2).min(axis=1)<10
if sum(mask)<10:
    raise ValueError("Too few mapped head-neighborhood positions")
rows={"legacy_dfg":sterics(legacy/"warhead.sdf",legacy/"protein.pdb")}
for label,include in (("whole_pocket",np.ones(len(keys),dtype=bool)),("head_neighborhood",mask)):
    r,t=placement.kabsch(x[include],y[include])
    transformed=xyz@r.T+t
    error=np.linalg.norm(x@r.T+t-y,axis=1)
    out=Path(label)
    out.mkdir()
    h=Chem.Mol(head)
    for i,p in enumerate(transformed):
        h.GetConformer().SetAtomPosition(i,p)
    with Chem.SDWriter(str(out/"poi_head.sdf")) as w:
        w.write(h)
    for f in ("vhl_head.sdf","vbc.pdb","full.smi","selection.json"):
        shutil.copy2(base/f,out/f)
    (out/"protein.pdb").write_text("\n".join(a["line"][:21]+"D"+a["line"][22:] for a in targets)+"\nTER\nEND\n")
    fields=qc.read_params(base/"Protac_params.txt")
    fields["Anchor atoms"]=fields["Anchor atoms"].split()[0]+" "+str(exits[0]+1)
    (out/"Protac_params.txt").write_text("".join(f"{k}: {v}\n" for k,v in fields.items()))
    result={"production_ready":False,"status":"MODELED_MAPK14_CONTROL_CANDIDATE","pose_reference":"5IA4 foretinib overlaid onto 1W82",
            "whole_pocket_rmsd_A":float(np.sqrt(np.mean(error**2))),
            "local_pocket_rmsd_A":float(np.sqrt(np.mean(error[mask]**2))),
            "fit_pocket_positions":int(include.sum()),"full_ligand_inchikey":Chem.MolToInchiKey(full),
            "poi_anchor":exits[0]+1,"ligand_full_anchor_one_based":next(iter(attachments[exits[0]]))+1,
            "sterics":sterics(out/"poi_head.sdf",out/"protein.pdb"),
            "source_hashes":{str(p):qc.sha256(p) for p in (references/"5IA4.pdb",references/"KLIFS_7189_pocket_residue_map.json",Path("1W82.pdb"),Path("1W82_pocket.json"))},
            "files":{p.name:qc.sha256(p) for p in out.iterdir()}}
    (out/"preparation.json").write_text(json.dumps(result,indent=2)+"\n")
    rows[label]=result
Path("MAPK14_INPUT_COMPARISON.json").write_text(json.dumps(rows,indent=2)+"\n")
for name,r in rows.items():
 print(name,r.get("whole_pocket_rmsd_A"),r.get("local_pocket_rmsd_A"),r.get("sterics",r)["max_vdw_overlap_A"])
