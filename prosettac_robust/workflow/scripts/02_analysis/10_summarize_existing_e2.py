#!/usr/bin/env python3
"""Summarize hash-bound retained E2 audit frames under explicit time/image rules.

No MD, new structures or GNN training. Keeps legacy ELOB/POI mask explicit.
"""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,sha256,atomic_json
from e2_summary import summarize_frames


def summarize(frames,pose_receipt,metrics,output,start_ns,end_ns,duplicates,image_policy):
    if output.exists() or any(output.resolve().is_relative_to(Path(p)) for p in ['/home/ayushregmibagale/protac-md','/mnt/d/PROTAC degradation Prediction']):raise ValidationError('Fresh output outside original work required')
    tracked=[frames,pose_receipt,metrics,Path(__file__),Path(__file__).resolve().parents[1]/'lib/e2_summary.py'];hashes={str(p):sha256(p) for p in tracked};receipt=json.loads(pose_receipt.read_text())
    if receipt.get('frames_sha256')!=hashes[str(frames)] or receipt.get('source_input_sha256',{}).get(str(metrics))!=hashes[str(metrics)]:raise ValidationError('E2 frame packet or metric input differs from audit receipt')
    if receipt.get('duplicate_coordinate_conflicts_in_audited_selection',1):raise ValidationError('Conflicting duplicate coordinates in E2 audit')
    with metrics.open() as f:rows=list(csv.DictReader(line for line in f if not line.startswith('#')))
    data=np.genfromtxt(frames,delimiter='\t',names=True,ndmin=1)
    report=summarize_frames(data,[float(r['time_ns']) for r in rows],[float(r['d_protac_vhl_A']) for r in rows],[float(r['d_protac_kinase_A']) for r in rows],start_ns,end_ns,duplicates,image_policy)
    if any(sha256(p)!=hashes[str(p)] for p in tracked):raise ValidationError('E2 summary input changed')
    report.update(poi=receipt['poi'],pose_id=receipt['pose_id'],input_sha256=hashes,old_outputs_modified=False,new_MD_or_GNN_training=False)
    output.mkdir(parents=True);atomic_json(output/'summary.json',report);return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['frames','pose-receipt','metrics','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--start-ns',type=float,required=True);p.add_argument('--end-ns',type=float,required=True);p.add_argument('--duplicates',choices=['reject','identical'],required=True);p.add_argument('--image-policy',choices=['reject','audit-only'],required=True);a=p.parse_args()
    print(json.dumps(summarize(a.frames,a.pose_receipt,a.metrics,a.output,a.start_ns,a.end_ns,a.duplicates,a.image_policy),indent=2))
