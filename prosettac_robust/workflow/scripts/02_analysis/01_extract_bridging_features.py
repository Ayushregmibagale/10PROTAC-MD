#!/usr/bin/env python3
"""Extract geometry from immutable, receipt-validated local attempts.

No guessed score columns, missing-feature sentinels or automatic ML predictions.
An explicit index links protac,gene,uid,attempt_id,local_dir. Output geometry is
not a complete input table for the existing MD/E2-dependent classifier.
"""
import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import verify_receipt, ValidationError

def extract_one(run):
    run=Path(run); selection=json.loads((run/'selection.json').read_text())
    rows=verify_receipt(run,selection)
    engine=json.loads((run/'engine_status.json').read_text())
    attempted=engine['stages']['linker_modeling']['attempted']
    shared={'n_protein_poses':attempted,'n_valid_ternary_poses':len(rows),
            'feature_schema':'validated_geometry_with_rejection_denominator_v2',
            'legacy_trained_model_compatible':False,
            'linker_geometry_yield':len(rows)/attempted if attempted else None,
            'protocol':engine.get('protocol','unspecified'),
            'outcome':engine['status']}
    if not rows:
        return [],dict(shared,n_poses=0,bridge_n_poses=0,bridge_pct=None,prostac_bridging_frac=None,
            bridge_frac_5A=None,bridge_frac_8A=None,mean_vhl_dist=None,mean_kin_dist=None,
            d_vhl_min=None,d_kin_min=None,d_sum_min=None,min_pose_d_sum=None,r_dvhl_dkin=None,
            kin_when_vhl_near=None,rosetta_score_mean=None,rosetta_score_bridge_mean=None)
    dv=np.array([r['d_vhl'] for r in rows]); dk=np.array([r['d_kin'] for r in rows])
    bridge=np.array([r['bridging_5A'] for r in rows]); scores=np.array([r['total_score'] for r in rows])
    return rows,dict(shared,**{'n_poses':len(rows),'bridge_n_poses':int(bridge.sum()),'bridge_pct':float(100*bridge.mean()),
                 'prostac_bridging_frac':float(bridge.mean()),'bridge_frac_5A':float(bridge.mean()),
                 'bridge_frac_8A':float(np.mean([r['bridging_8A'] for r in rows])),
                 'mean_vhl_dist':float(dv.mean()),'mean_kin_dist':float(dk.mean()),
                 'd_vhl_min':float(dv.min()),'d_kin_min':float(dk.min()),
                 'd_sum_min':float(dv.min()+dk.min()),
                 'min_pose_d_sum':float((dv+dk).min()),
                 'r_dvhl_dkin':float(np.corrcoef(dv,dk)[0,1]) if len(rows)>1 and dv.std()>0 and dk.std()>0 else None,
                 'kin_when_vhl_near':float(dk[dv<8].mean()) if np.any(dv<8) else None,
                 'rosetta_score_mean':float(scores.mean()),
                 'rosetta_score_bridge_mean':float(scores[bridge].mean()) if bridge.any() else None})

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); summaries=[]; poses=[]; errors=[]; seen=set()
    for row in csv.DictReader(a.index.open(),delimiter='\t'):
        try:
            key=tuple(row[k] for k in ('protac','gene','uid','attempt_id'))
            if key in seen: raise ValidationError('Duplicate attempt in index')
            seen.add(key)
            data,summary=extract_one(row['local_dir'])
            identity={k:row[k] for k in ('protac','gene','uid','attempt_id')}
            summaries.append(dict(identity,**summary))
            poses.extend(dict(identity,**r) for r in data)
        except (ValidationError,OSError,ValueError,KeyError) as exc:
            errors.append({'run':row,'error':str(exc)})
    a.output.mkdir(parents=True,exist_ok=True)
    from protac_qc import atomic_json
    atomic_json(a.output/'extraction_errors.json',errors)
    # Never replace the authoritative table with a partial, success-only subset.
    if errors: raise SystemExit(f'{len(errors)} invalid attempts; no feature tables written. See extraction_errors.json')
    if not summaries: raise SystemExit('No validated attempts')
    for name,records in [('bridging_features_all.tsv',summaries),('pose_geometry.tsv',poses)]:
        if not records:continue
        with (a.output/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t'); w.writeheader(); w.writerows(records)
    print(f'{len(summaries)} attempts, {len(poses)} validated poses. ML schema validation remains separate.')

if __name__=='__main__': main()
