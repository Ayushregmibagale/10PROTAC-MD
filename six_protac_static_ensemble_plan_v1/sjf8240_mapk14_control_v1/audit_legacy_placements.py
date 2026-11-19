"""Measure geometry of the legacy DFG candidate placements without replacing files."""
from collections import Counter
import csv
import json
from pathlib import Path
import runpy
import sys

legacy,helper=map(Path,sys.argv[1:])
sterics=runpy.run_path(str(helper))["sterics"]
out=[]
with (legacy/"prep_log.tsv").open() as handle:
    for report in csv.DictReader(handle,delimiter="\t"):
        if report["status"]!="ok":
            continue
        folder=legacy/report["uid"]
        try:
            scores=sterics(folder/"warhead.sdf",folder/"protein.pdb")
        except (OSError,ValueError,RuntimeError) as error:
            out.append({"uid":report["uid"],"error":str(error)})
            continue
        out.append({"uid":report["uid"],"gene":report["gene"],
                    "dfg_fit_rmsd_A":float(report["superimpose_rms"]),**scores})
valid=[row for row in out if "error" not in row]
sample={"archived_prepared":len(out),"geometry_measured":len(valid),
        "minimum_contact_below_1p5_A":sum(row["min_distance_A"]<1.5 for row in valid),
        "maximum_vdw_overlap_above_0p6_A":sum(row["max_vdw_overlap_A"]>.6 for row in valid),
        "failure_counts":dict(Counter(row["error"] for row in out if "error" in row)),
        "rows":out}
Path("LEGACY_GEOMETRY_AUDIT.json").write_text(json.dumps(sample,indent=2)+"\n")
print({k:v for k,v in sample.items() if k!="rows"})
