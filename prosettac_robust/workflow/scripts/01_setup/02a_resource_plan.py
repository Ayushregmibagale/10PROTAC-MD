#!/usr/bin/env python3
"""Measured baseline storage projection; read-only except the requested report.

Peak memory/time require representative measurement; this tool does not invent
these values or authorize a production campaign.
"""
import argparse,json,gzip,shutil,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import atomic_json,sha256,ValidationError

def plan(inputs,pilot,n_dock,n_local,conformers,available_bytes=None,pose_archive=None):
    if min(n_dock,n_local,conformers)<1:raise ValidationError('Sampling budgets must be positive')
    pairs=[d for d in inputs.glob('*/*') if (d/'preparation.json').exists()]
    if not pairs:raise ValidationError('No prepared inputs')
    pilot_dir=pilot/'Patchdock_Results'
    local=list(pilot_dir.glob('pd.*_docking_????.pdb'))
    if not local:raise ValidationError('No measured local-docking output')
    parameters=[json.loads(p.read_text()) for p in pilot_dir.glob('PT_*_parameter_provenance.json')]
    sampling=list(pilot_dir.glob('confs_*.sdf.sampling.json'))
    if sampling:
        if len(sampling)!=len(local):raise ValidationError('Missing per-local sampling budget evidence')
        generated={};actual=[];absent_zero_ensembles=[]
        for p in sampling:
            ledger=json.loads(p.read_text());count=ledger['generated_conformers']
            if ledger['requested_conformers']!=conformers or not 0<=count<=conformers:
                raise ValidationError('Measured requested conformer budget differs from projection')
            actual.append(count)
            sdf=Path(str(p)[:-len('.sampling.json')])
            if not sdf.exists() and count==0:
                absent_zero_ensembles.append(sdf.name);continue
            if not sdf.exists():raise ValidationError('A generated nonempty ensemble is missing')
            digest=sha256(sdf)
            if digest in generated and generated[digest]!=count:raise ValidationError('Conflicting conformer counts')
            generated[digest]=count
        for record in parameters:
            if generated.get(record['source_sha256'])!=record['input_conformers']:
                raise ValidationError('Parameter ensemble differs from recorded generated conformers')
        if {p['source_sha256'] for p in parameters}!={k for k,v in generated.items() if v>0}:
            raise ValidationError('Missing parameter provenance for a feasible ensemble')
        budget_basis='Explicit requested budgets; valid generated counts may be lower'
    else:
        budgets={p['input_conformers'] for p in parameters}
        if budgets!={conformers}:raise ValidationError('Measured full-conformer budget differs from requested budget')
        actual=[p['input_conformers'] for p in parameters];budget_basis='Historical generated-count evidence only';absent_zero_ensembles=[]
    # Include every measured local-stage artifact, not just the final scored PDB.
    files=[p for p in pilot_dir.rglob('*') if p.is_file()]
    raw=sum(p.stat().st_size for p in files)
    compressed=sum(len(gzip.compress(p.read_bytes(),compresslevel=6,mtime=0)) for p in files)
    free=available_bytes if available_bytes is not None else shutil.disk_usage(inputs).free
    scale=len(pairs)*n_dock*n_local/len(local)
    total_raw=round(raw*scale);total_gzip=round(compressed*scale)
    report={'status':'STORAGE_PROJECTION_ONLY','pairs':len(pairs),'n_dock':n_dock,'n_local':n_local,
        'resource_planner_sha256':sha256(Path(__file__)),
        'conformer_budget_evidence':budget_basis,'actual_generated_conformers_by_local_pose':actual,
        'zero_conformer_files_absent_after_upstream_empty_file_cleanup':absent_zero_ensembles,
        'requested_full_conformers':conformers,'pilot_local_pose_count':len(local),
        'pilot_conformer_budget_must_match':True,
        'pilot_total_intermediate_bytes':raw,'pilot_gzip6_bytes':compressed,
        'estimated_all_pairs_raw_bytes':total_raw,'estimated_all_pairs_gzip6_bytes':total_gzip,
        'available_bytes':free,'available_source':'explicit observation' if available_bytes is not None else 'local filesystem observation',
        'reserve_fraction':0.30,'raw_fits_with_reserve':total_raw<=0.7*free,
        'gzip_estimate_fits_with_reserve':total_gzip<=0.7*free,
        'limitations':['pilot conformer budget must match requested budget',
          'compressed estimate excludes full raw-plus-archive temporary overlap',
          'protein-size/chemotype variability, MD trajectories, retries and backups require additional budget',
          'RAM, walltime and concurrent-task limits require a representative pilot',
          'no files deleted or archived by this projection'],
        'pilot_artifact_sha256':{str(p.relative_to(pilot)):sha256(p) for p in files}}
    if pose_archive is not None:
        archive_bytes=sum(p.stat().st_size for p in pose_archive.iterdir() if p.is_file())
        metadata=[p for p in files if p.suffix in {'.json','.params','.log','.fasc'} or p.name=='score.sc']
        metadata_bytes=sum(len(gzip.compress(p.read_bytes(),compresslevel=6,mtime=0)) for p in metadata)
        retained=round((archive_bytes+metadata_bytes)*scale)
        report.update(pilot_lossless_pose_archive_bytes=archive_bytes,
            estimated_retained_archive_plus_compressed_metadata_bytes=retained,
            retained_estimate_fits_with_reserve=retained<=0.7*free,
            retention_condition='All scored/local/global poses archived and restoration SHA256 verified; retain params, scores, logs, chemistry, exact code/environment/seeds; remove only reproducible temporary artifacts under a separately validated policy',
            retention_policy_executed=False)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['inputs','pilot','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--n-dock',type=int,required=True);p.add_argument('--n-local',type=int,required=True)
    p.add_argument('--conformers',type=int,required=True);p.add_argument('--available-bytes',type=int)
    p.add_argument('--pose-archive',type=Path)
    a=p.parse_args();report=plan(a.inputs,a.pilot,a.n_dock,a.n_local,a.conformers,a.available_bytes,a.pose_archive)
    atomic_json(a.output,report);print(json.dumps({k:v for k,v in report.items() if k!='pilot_artifact_sha256'},indent=2))
