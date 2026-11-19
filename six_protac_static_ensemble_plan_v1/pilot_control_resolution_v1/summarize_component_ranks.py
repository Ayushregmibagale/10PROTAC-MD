"""Join saved component scores to a frozen GSK215 reference recovery table."""
import argparse
import csv
import json
from pathlib import Path

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("scores",type=Path)
parser.add_argument("recovery",type=Path)
parser.add_argument("output",type=Path)
args=parser.parse_args()
with args.recovery.open() as handle:
    native={row["pose"]:row for row in csv.DictReader(handle,delimiter="\t")
            if row["reference"]=="7PI4_GSK215"}
scores=json.loads(args.scores.read_text())
selected={row["pose"]:row for row in scores if row["pose"] in native}
if set(native)!=set(selected):
    raise ValueError("Missing scores for retained reference poses")
rows=[{**row,"raw_score":float(native[name]["raw_total_score"]),
       "reference_gate_pass":native[name]["frozen_gate_pass"]=="True"}
      for name,row in selected.items()]
orders={}
for metric in ("raw_score","protein_pair_delta","protac_protein_delta"):
    ranked=sorted(rows,key=lambda row:(row[metric],row["pose"]))
    orders[metric]={"passing_ranks":[i for i,row in enumerate(ranked,1)
                                    if row["reference_gate_pass"]],"top_id":ranked[0]["pose"]}
summary={"raw_score_max_absolute_replay_difference":max(abs(r["raw_score"]-r["all"]) for r in rows),
         "orders":orders,"poses":rows}
args.output.write_text(json.dumps(summary,indent=2)+"\n")
print(orders)
