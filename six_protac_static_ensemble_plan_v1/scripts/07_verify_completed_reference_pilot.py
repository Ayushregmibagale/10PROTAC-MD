"""Independently verify completed pilot source/input/output receipt integrity.

Does not certify physical pose quality or degradation prediction performance.
"""
import argparse,hashlib,json,math
from pathlib import Path

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def scores(p):
    rows={}
    for line in p.read_text().splitlines():
        x=line.split()
        if len(x)==3 and x[0]=='SCORE:' and x[1]!='total_score':
            value=float(x[1]);name=x[2]
            if not math.isfinite(value) or name in rows:raise ValueError('Invalid/duplicate score')
            rows[name]=value
    return rows
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pilot',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--completed-only',action='store_true',help='Verify immutable completed tasks while other tasks still run; never certify the whole pilot')
    a=p.parse_args()
    source=json.loads((a.pilot/'BUNDLE_MANIFEST.json').read_text());state=json.loads((a.pilot/'STATUS.json').read_text())
    complete=state['status']=='REFERENCE_PILOT_COMPLETE_NEEDS_RECOVERY_ASSESSMENT'
    if not complete and not a.completed_only:raise ValueError('Pilot not complete')
    tasks=[t for t in state['tasks'] if t['status']=='COMPLETE'] if a.completed_only else state['tasks']
    if not tasks:raise ValueError('No completed tasks to verify')
    receipts=[]
    for task in tasks:
        run=a.pilot/'runs'/task['name'];raw=run/'engine_raw';engine=json.loads((run/'engine_status.json').read_text())
        if sha(run/'engine_status.json')!=task['engine_status_sha256']:raise ValueError('Engine receipt changed')
        if task['exit_code']!=0 or task['status']!='COMPLETE':raise ValueError('Task was not successfully completed')
        for rel,digest in engine['engine_sha256'].items():
            if source['files'].get('workflow/'+rel)!=digest:raise ValueError('Engine source differs from frozen bundle: '+rel)
        ref=next(t['reference'] for t in state['configuration']['tasks'] if t['name']==task['name'])
        for filename in ('full.smi','selection.json','vhl_head.sdf','poi_head.sdf','protein.pdb','vbc.pdb','Protac_params.txt'):
            if sha(run/filename)!=source['files']['inputs/'+ref+'/'+filename]:raise ValueError('Frozen input changed: '+filename)
        if (run/'native_complex.pdb').exists():raise ValueError('Unexpected native ternary evaluation coordinates in generation inputs')
        if sha(raw/'Patchdock_Results/score.sc')!=engine['raw_score_sha256']:raise ValueError('Raw score changed')
        if sha(raw/'baseline_output_qc.json')!=engine['raw_output_qc_sha256']:raise ValueError('Raw QC changed')
        for name,digest in engine['sampling_sha256'].items():
            if Path(name).name!=name or sha(raw/'Patchdock_Results'/name)!=digest:raise ValueError('Sampling ledger changed')
        quality=json.loads((raw/'baseline_output_qc.json').read_text());original=scores(raw/'Patchdock_Results/score.sc');retained=scores(run/'Patchdock_Results/score.sc')
        qc={r['description']:r for r in quality['poses']}
        if set(qc)!=set(original):raise ValueError('Raw QC and score membership differ')
        for name,row in qc.items():
            if sha(raw/'Patchdock_Results'/f'{name}.pdb')!=row['pdb_sha256']:raise ValueError('Raw pose changed')
        passed={name for name,r in qc.items() if r['status']=='MECHANICAL_QC_PASSED'}
        if set(retained)!=passed or len(retained)!=engine['stages']['pose_qc']['succeeded']:raise ValueError('Retained pose membership inconsistent')
        for name,value in retained.items():
            row=qc[name]
            if abs(value-original[name])>1e-6:raise ValueError('Retained score shifted/changed')
            if sha(run/'Patchdock_Results'/f'{name}.pdb')!=row['pdb_sha256']:raise ValueError('Retained pose differs from raw')
            head=row['head_engagement']
            if head['status']!='MEASURED_NOT_A_BINDING_LABEL' or max(head['vhl_core_rmsd_A'],head['poi_core_rmsd_A'])>.05:raise ValueError('Bound-head preservation check failed')
        if len(original)!=engine['stages']['scoring']['succeeded']:raise ValueError('Raw scoring count differs')
        run_bytes=sum(f.stat().st_size for f in run.rglob('*') if f.is_file())
        receipts.append({'name':task['name'],'reference':ref,'status':'INTEGRITY_PASS','engine_status':engine['status'],
            'raw_scored_poses':len(original),'retained_operational_QC_poses':len(retained),'sampling_ledgers':len(engine['sampling_sha256']),
            'run_logical_bytes':run_bytes,'wall_seconds':task['finished_unix']-task['started_unix'],
            'cpu_affinity_limit':state['configuration']['cpu_limit'],'peak_summed_group_RSS_bytes':task['peak_group_RSS_bytes'],
            'observed_environment':engine['environment'],'engine_status_sha256':sha(run/'engine_status.json'),
            'physical_pose_quality_certified':False})
    a.output.mkdir(parents=True,exist_ok=True)
    report={'status':'SOURCE_INPUT_OUTPUT_INTEGRITY_PASS' if complete else 'COMPLETED_TASK_INTEGRITY_PASS','pilot_status':state['status'],
      'remaining_tasks':[t['name'] for t in state['tasks'] if t['status']!='COMPLETE'],'runs':receipts,'MD':False,
      'limitations':['Only receipt/source/input/output integrity assessed','Operational QC is coarse; structural recovery and independent atom-type-aware steric checks remain separate',
         'Peak output in the controller is cumulative across runs; per-run bytes above are independently counted after completion',
         'Both references share the same seed; this is not a replicated convergence benchmark']}
    (a.output/'INTEGRITY_AND_RESOURCES.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({r['name']:{k:r[k] for k in ('status','raw_scored_poses','retained_operational_QC_poses','run_logical_bytes','wall_seconds')} for r in receipts},indent=2))

if __name__=='__main__':main()
