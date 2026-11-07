#!/usr/bin/env python3
"""Mapped heavy-atom contacts with triclinic minimum-image distances.

All distances are in Angstrom. Block summaries describe within-trajectory
variation, not independent biological replicates or proof of convergence.
"""
import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,atomic_json,sha256
from md_mapping import mapped_heavy_groups

def validate_analysis_times(discard_ps,block_ps):
    if not np.isfinite(discard_ps) or not np.isfinite(block_ps) or discard_ps<0 or block_ps<=0:
        raise ValidationError('Explicit finite nonnegative discard and positive block duration required')

def validate_analysis_stage(completed,diagnostic):
    if diagnostic:
        if completed['status'] not in {'SMOKE_COMPLETE','COMPLETE'}:
            raise ValidationError('Diagnostic contacts require a completed stage')
    elif completed['status']!='COMPLETE' or completed['intent']['stage']!='prod':
        raise ValidationError('Scientific contacts require completed production MD')


def analyze(preparation,topology,trajectory,receipt,output,discard_ps,block_ps,diagnostic=False):
    validate_analysis_times(discard_ps,block_ps)
    tracked=[preparation,topology,trajectory,receipt]
    if topology.with_suffix('.provenance.json').exists():tracked.append(topology.with_suffix('.provenance.json'))
    before={str(p):sha256(p) for p in tracked}
    meta=json.loads(preparation.read_text());completed=json.loads(receipt.read_text())
    validate_analysis_stage(completed,diagnostic)
    if completed['intent']['inputs'].get('md_preparation.json')!=sha256(preparation):
        raise ValidationError('Atom mapping/preparation differs from executed system')
    if completed['outputs'].get(trajectory.name)!=sha256(trajectory):raise ValidationError('Trajectory differs from stage receipt')
    if completed['outputs'].get(topology.name)!=sha256(topology):
        converted=json.loads(topology.with_suffix('.provenance.json').read_text())
        if (converted['group']!='SOLUTE' or converted['target_sha256']!=sha256(topology)
            or converted['source_sha256'] not in [v for k,v in completed['outputs'].items() if k.endswith('.tpr')]):
            raise ValidationError('Reduced topology does not descend from executed stage')
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array
    from read_only_trajectory import ReadOnlyXTCReader
    if trajectory.suffix.lower()!='.xtc':raise ValidationError('Qualified read-only reader requires an XTC trajectory')
    universe=mda.Universe(str(topology),str(trajectory),format=ReadOnlyXTCReader)
    groups=mapped_heavy_groups(universe,meta)
    rows=[]
    last_time=None
    for ts in universe.trajectory:
        if not np.isfinite(ts.time) or (last_time is not None and ts.time<=last_time):
            raise ValidationError('Trajectory times are nonfinite, repeated or decreasing')
        last_time=ts.time
        if ts.time<discard_ps:continue
        box=np.asarray(ts.dimensions,dtype=float)
        if box.shape!=(6,) or not np.isfinite(box).all() or np.any(box[:3]<=0) or np.any((box[3:]<=0)|(box[3:]>=180)):
            raise ValidationError('Periodic box absent or invalid')
        if any(not np.isfinite(g.positions).all() for g in groups.values()):
            raise ValidationError('Nonfinite mapped coordinates')
        dv=float(distance_array(groups['protac'].positions,groups['vhl'].positions,box=ts.dimensions).min())
        dp=float(distance_array(groups['protac'].positions,groups['poi'].positions,box=ts.dimensions).min())
        rows.append({'time_ps':float(ts.time),'d_vhl_A':dv,'d_poi_A':dp,
                     'bridging_5A':int(dv<5 and dp<5),'bridging_8A':int(dv<8 and dp<8)})
    if not rows:raise ValidationError('No trajectory frames remain after discard')
    if any(sha256(p)!=before[str(p)] for p in tracked):raise ValidationError('Contact-analysis input changed during extraction')
    times=np.array([r['time_ps'] for r in rows])
    if len(times)>1 and np.any(np.diff(times)<=0):raise ValidationError('Trajectory times repeat or decrease')
    # Complete physical-time blocks only; the final short tail is excluded.
    nblocks=int(np.floor((times[-1]-discard_ps)/block_ps));blocks=[]
    for i in range(nblocks):
        selected=[r for r in rows if discard_ps+i*block_ps<=r['time_ps']<discard_ps+(i+1)*block_ps]
        if selected:blocks.append({'block':i,'n_frames':len(selected),
                                  'bridge_fraction_5A':float(np.mean([r['bridging_5A'] for r in selected]))})
    output.mkdir(parents=True,exist_ok=False)
    with (output/'contacts.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    summary={'status':'DIAGNOSTIC_ONLY' if diagnostic else 'PRODUCTION_CONTACT_SUMMARY',
        'units':{'distance':'Angstrom','time':'ps'},'periodic_method':'triclinic minimum-image heavy-atom distance',
        'n_frames':len(rows),'discard_ps':discard_ps,'block_ps':block_ps,'blocks':blocks,
        'bridge_fraction_5A':float(np.mean([r['bridging_5A'] for r in rows])),
        'bridge_fraction_8A':float(np.mean([r['bridging_8A'] for r in rows])),
        'convergence_certified':False,'independent_replica_uncertainty':None,
        'input_sha256':before,
        'MDAnalysis_version':mda.__version__}
    atomic_json(output/'summary.json',summary);return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['preparation','topology','trajectory','receipt','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--discard-ps',type=float,required=True);p.add_argument('--block-ps',type=float,required=True)
    p.add_argument('--diagnostic',action='store_true');a=p.parse_args()
    print(json.dumps(analyze(a.preparation,a.topology,a.trajectory,a.receipt,a.output,a.discard_ps,a.block_ps,a.diagnostic),indent=2))
