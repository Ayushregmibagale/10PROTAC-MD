"""Reconcile the full six-compound panel with verified, held candidate evidence."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import shutil
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
for key in ("project", "study", "dfg", "output"):
    parser.add_argument(key, type=Path)
parser.add_argument("--held-bundle", type=Path)
args = parser.parse_args()

def load(path):
    return json.loads(path.read_text())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def verify(path, hashes):
    if not hashes:
        raise ValueError(f"Missing input hashes: {path}")
    for name, digest in hashes.items():
        if Path(name).name != name or sha(path / name) != digest:
            raise ValueError(f"Input mismatch: {path}/{name}")

inventory = args.project / "structure_provenance/domains/domain_inventory.tsv"
with inventory.open() as handle:
    genes = {row["uniprot_id"]:row["gene"] for row in csv.DictReader(handle,delimiter="\t")}
registry = load(args.project / "config/molecules.json")["compounds"]
draft = args.study / "ernie_bulk_submission_draft_v2"
draft_manifest = draft / "release.json"
if load(draft_manifest)["qualified"]:
    raise ValueError("Historical bulk draft is unexpectedly qualified")
all_pairs={(drug, uid) for drug in registry for uid in genes}
known=defaultdict(list)
locations={}
sources=[draft_manifest]

def add(key, label, path, hashes, detail, primary=False):
    if key not in all_pairs:
        raise ValueError(f"Unexpected panel pair {key}")
    path = path if path.is_absolute() else args.study / path
    verify(path,hashes)
    locations[(key,label,path.name)]=path
    known[key].append({"source":label,"input_id":path.name,
                       "input_sha256":hashes,"detail":detail,"primary_control":primary})

for task in load(draft_manifest)["tasks"]:
    add((task["compound"],task["uniprot"]),"held_ernie_draft_v2",draft/task["input"],
        task["input_sha256"],task["input_review"])
latest_gap=args.study/"ernie_candidate_manifest_v1/supplemental_af2_whole_v1/UPDATED_GAP_LEDGER.tsv"
sources.append(latest_gap)
with latest_gap.open() as handle:
    historical_gaps={(row["compound"],row["uniprot"]):row for row in csv.DictReader(handle,delimiter="\t")}
for folder, filename in (("pocket_local_fit_v1","COMBINED_RECOVERY_REVIEW.json"),
                         ("supplemental_af2_whole_v1","SUPPLEMENTAL_RECOVERY_REVIEW.json")):
    file=args.study/"ernie_candidate_manifest_v1"/folder/filename
    sources.append(file)
    for row in load(file)["candidates"]:
        add((row["compound"],row["uniprot"]),row["method"],Path(row["source"]),
            row["input_sha256"],"Pocket-transfer candidate under review")
dfg_receipt=args.dfg/"DFG_CANDIDATE_RECEIPTS.json"
sources.append(dfg_receipt)
for row in load(dfg_receipt):
    add((row["compound"],row["uniprot"]),"legacy_DFG_checked_v1",args.dfg.resolve()/row["input"],
        row["input_sha256"],"Supplier graph, both anchors, and geometry checked; structural review pending")
mapk_dir=args.study/"sjf8240_mapk14_control_v1"
mapk_input=mapk_dir/"local_archive/training_input"
mapk_preparation=mapk_input/"preparation.json"
mapk_result=mapk_dir/"MAPK14_PILOT_REVIEW.json"
sources.extend((mapk_preparation,mapk_result))
review=load(mapk_result)
if review["status"]!="COMPLETE" or len(review["poses"])!=12 or any(x["identity"]!="PASS" for x in review["poses"]):
    raise ValueError("MAPK14 control check is incomplete")
add(("SJF8240","Q16539"),"mapk14_modeled_control_v1",mapk_input,
    {**load(mapk_preparation)["files"],"preparation.json":sha(mapk_preparation)},
    "Primary SJF8240 technical control; 12 independently verified retained poses",primary=True)

notes={
    ("GSK215","Q05397"):"Native-like poses recovered; alternative score exploratory",
    ("SJF8240","P08581"):"Published MET target; zero feasible linker at tested budgets",
    ("SJF8240","Q16539"):"Primary SJF8240 modeled technical control; no native ternary reference",
    ("SJF_delta","O15264"):"Published p38δ target; modeled ternary remains under steric review",
    ("SJF1521","P00533"):"Modeled EGFR head; anchors verified; production_ready false",
    ("JNJ1013","P51617"):"Repaired IRAK1 control; three scored poses passed corrected QC",
    ("MS143","P31749"):"Repaired AKT1 control; 27 scored poses passed corrected QC with TPO308",
}
rows=[]
index=[]
for drug,uid in sorted(all_pairs):
    records=known.get((drug,uid),[])
    if any(r["primary_control"] for r in records):
        status="MODELED_CONTROL_VERIFIED"
    elif any(r["source"]=="held_ernie_draft_v2" for r in records):
        status="STAGED_CANDIDATE_HELD"
    elif records:
        status="PREPARED_CANDIDATE_UNDER_REVIEW"
    else:
        status="MISSING_PREPARED_CANDIDATE"
    detail=notes.get((drug,uid),"")
    if not records and (drug,uid) in historical_gaps:
        detail=(detail+"; " if detail else "")+"No usable prepared candidate in the current ledger"
    row={"compound":drug,"uniprot":uid,"gene":genes[uid],"input_status":status,"release":"HOLD",
         "candidate_count":len(records),"evidence_sources":"; ".join(r["source"] for r in records),
         "note":detail}
    rows.append(row)
    if records:
        index.append({**row,"candidates":records})

if len(rows)!=2610 or len({(r["compound"],r["uniprot"]) for r in rows})!=2610:
    raise ValueError("Expected six × 435 unique pair rows")
args.output.mkdir(parents=True,exist_ok=False)
with (args.output/"RELEASE_CHECKLIST.tsv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
(args.output/"HELD_CANDIDATE_INDEX.json").write_text(json.dumps({
    "dispatch_permitted":False,"production_ready":False,"scope":"Wider verified candidate inventory",
    "source_sha256":{str(f.relative_to(args.study)) if f.is_relative_to(args.study) else f.name:sha(f) for f in sources},
    "pairs":index},indent=2)+"\n")
summary={"full_panel_pairs":len(rows),"pairs_with_candidate":len(index),
         "states":dict(Counter(row["input_status"] for row in rows)),
         "by_compound":{name:dict(Counter(row["input_status"] for row in rows if row["compound"]==name))
                        for name in registry},
         "release_approved_pairs":0,"production_ready":False}
(args.output/"RELEASE_CHECKLIST_SUMMARY.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))

if args.held_bundle:
    bundle=args.held_bundle
    bundle.mkdir(parents=True, exist_ok=False)
    repair=args.study/"v3_conformer_repair_v1/runtime_bundle"
    plant=args.study/"ernie_bulk_submission_draft_v2"
    shutil.copytree(repair,bundle/"runtime_bundle",ignore=shutil.ignore_patterns("__pycache__"))
    audit_fix=args.study/"pilot_control_review_v1/audit_outputs_candidate.py"
    shutil.copy2(audit_fix,bundle/"runtime_bundle/baseline_engine/audit_outputs.py")
    for name in ("submit.py","run_queue.py"):
        shutil.copy2(plant/name,bundle/name)
    tasks=[]
    for row in index:
        key=(row["compound"],row["uniprot"])
        candidates=row["candidates"]
        selected=next((c for c in candidates if c["primary_control"]),candidates[0])
        source=locations[(key,selected["source"],selected["input_id"])]
        task_id="__".join(key)
        destination=Path("inputs")/key[0]/key[1]
        directory=bundle/destination
        directory.mkdir(parents=True)
        for name in selected["input_sha256"]:
            shutil.copy2(source/name,directory/name)
        tasks.append({"id":task_id,"input":str(destination),"input_sha256":selected["input_sha256"],
                      "seed":int.from_bytes(hashlib.sha256(task_id.encode()).digest()[:4],"big")%2147483647,
                      "input_status":row["input_status"],"source":selected["source"]})
    runtime_hashes={str(f.relative_to(bundle)):sha(f) for f in (bundle/"runtime_bundle").rglob("*") if f.is_file()}
    for name in ("submit.py","run_queue.py"):
        runtime_hashes[name]=sha(bundle/name)
    package={"qualified":False,"scope":"Held subset covering the available candidate inputs in the 2610-pair panel",
             "intended_full_panel_pairs":len(rows),"missing_pair_count":len(rows)-len(tasks),
             "checklist_sha256":sha(args.output/"RELEASE_CHECKLIST.tsv"),
             "budgets":load(draft_manifest)["budgets"],"tasks":tasks,"source_sha256":runtime_hashes}
    (bundle/"release.json").write_text(json.dumps(package,indent=2)+"\n")
    print(f"Frozen held package with {len(tasks)} tasks")
