"""Read-only compiled-stage and completion-log evidence, not new MD."""
from pathlib import Path
import collections,concurrent.futures,json,re,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,atomic_json,ValidationError
PROJECT=Path('/home/ayushregmibagale/protac-md');GMX=Path('/home/ayushregmibagale/.conda/envs/gromacs/bin/gmx')
STAGES={'em':'em','nvt':'nvt','npt_restrained':'npt_restrained','npt':'npt'}


def log_evidence(path):
    if not path.exists():return {'available':False}
    before=sha256(path);text=path.read_text(errors='replace');steps=re.findall(r'^\s*Step\s+Time\s*\n\s*(\d+)\s+([-+\d.eE]+)',text,re.M);force=re.findall(r'^\s*Maximum force\s*=\s*(\S+)',text,re.M)
    value={'available':True,'path':str(path),'sha256':before,'finished_marker':'Finished mdrun' in text,'last_logged_step':int(steps[-1][0]) if steps else None,'last_logged_time_ps':float(steps[-1][1]) if steps else None,'LINCS_warning_count':len(re.findall('LINCS WARNING',text,re.I)),'constraint_error_count':len(re.findall('Constraint error',text,re.I)),'fatal_error_count':len(re.findall('Fatal error:',text,re.I)),'positive_EM_convergence_statement':bool(re.search(r'^\s*(?:Steepest Descents|Conjugate Gradients|L-BFGS) converged to Fmax <',text,re.M)),'last_maximum_force':force[-1] if force else None}
    if sha256(path)!=before:raise ValidationError('Stage log changed while reading')
    return value


def inspect(case,output):
    poi=case['poi'];pose=case['pose_id'];base=PROJECT/'random_screen'/poi/pose;rows=[]
    for stage,prefix in STAGES.items():
        folder=base/stage;tpr=folder/(prefix+'.tpr');entry={'poi':poi,'pose_id':pose,'stage':stage,'TPR_available':tpr.exists(),'log':log_evidence(folder/(prefix+'.log'))}
        if tpr.exists():
            before=sha256(tpr);target=output/'inputrec'/poi/pose;target.mkdir(parents=True,exist_ok=True);mdp=target/(stage+'_inspection.mdp');log=target/(stage+'_inspection.log');command=[str(GMX),'dump','-s',str(tpr),'-om',str(mdp)]
            with log.open('w') as handle:result=subprocess.run(command,stdout=subprocess.DEVNULL,stderr=handle,text=True,timeout=60)
            if result.returncode:raise ValidationError('Input-record inspection failed: '+str(log))
            values={}
            for line in mdp.read_text().splitlines():
                if '=' in line:
                    k,v=line.split('=',1);values[k.strip().replace('_','-')]=' '.join(v.split())
            selected={k:values.get(k) for k in ['integrator','dt','nsteps','emtol','emstep','tcoupl','ref-t','tau-t','pcoupl','ref-p','tau-p','rcoulomb','rvdw','constraint-algorithm','continuation']}
            entry.update(tpr=str(tpr),tpr_sha256=before,compiled_settings=selected,inspection_sha256=sha256(mdp),inspection_log_sha256=sha256(log),nominal_duration_ps=float(values['dt'])*int(values['nsteps']) if values['integrator'] in {'md','sd','bd','md-vv','md-vv-avek'} else None)
            if sha256(tpr)!=before:raise ValidationError('TPR changed while inspecting')
        rows.append(entry)
    production=log_evidence(base/'prod/prod.log');production.update(poi=poi,pose_id=pose)
    return rows,production


def main(output):
    if output.exists():raise ValidationError('Fresh stage audit output required')
    roster=ROOT/'sjf8240_review/md_graph_pose_lineage.json';source=Path(__file__);hashes={str(p):sha256(p) for p in [roster,source,GMX]};cases=json.loads(roster.read_text())
    if len(cases)!=174:raise ValidationError('Expected 174 retained production cases')
    output.mkdir(parents=True);start=time.monotonic();stages=[];production=[];failures=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(inspect,c,output):c for c in cases}
        for future in concurrent.futures.as_completed(futures):
            try:a,b=future.result();stages.extend(a);production.append(b)
            except Exception as exc:failures.append({'poi':futures[future]['poi'],'pose_id':futures[future]['pose_id'],'error':str(exc)})
            if (len(production)+len(failures))%30==0:print(json.dumps({'cases':len(production),'failures':len(failures),'elapsed_s':round(time.monotonic()-start,1)}),flush=True)
    if any(sha256(Path(p))!=h for p,h in hashes.items()):raise ValidationError('Source changed during stage audit')
    stages.sort(key=lambda r:(r['poi'],r['pose_id'],r['stage']));production.sort(key=lambda r:(r['poi'],r['pose_id']));atomic_json(output/'per_stage.json',stages);atomic_json(output/'production_logs.json',production)
    per_stage={}
    for stage in STAGES:
        selected=[r for r in stages if r['stage']==stage];available=[r for r in selected if r['TPR_available']]
        per_stage[stage]={'requested':len(selected),'compiled_TPR_available':len(available),'stage_log_available':sum(r['log']['available'] for r in selected),'nominal_duration_ps_counts':dict(collections.Counter(str(r['nominal_duration_ps']) for r in available)),'ref_t_counts':dict(collections.Counter(str(r['compiled_settings']['ref-t']) for r in available)),'positive_EM_convergence_statements':sum(r['log'].get('positive_EM_convergence_statement',False) for r in selected) if stage=='em' else None}
    summary={'status':'COMPLETE_STAGE_EVIDENCE_AUDIT' if not failures else 'INCOMPLETE_STAGE_EVIDENCE_AUDIT','cases':len(production),'per_stage':per_stage,'production_log_available':sum(r['available'] for r in production),'production_finished_markers':sum(r.get('finished_marker',False) for r in production),'production_last_step_counts':dict(collections.Counter(str(r.get('last_logged_step')) for r in production)),'production_cases_with_LINCS_warnings':sum(r.get('LINCS_warning_count',0)>0 for r in production),'production_cases_with_fatal_error_marker':sum(r.get('fatal_error_count',0)>0 for r in production),'failures':failures,'elapsed_s':time.monotonic()-start,'source_sha256':hashes,'stage_records_sha256':sha256(output/'per_stage.json'),'production_records_sha256':sha256(output/'production_logs.json'),'original_files_modified':False,'new_MD_executed':False,'limitations':['Absent intermediate files in this retained layout do not establish that equilibration was omitted.','Requested TPR steps and log completion markers are evidence of execution, not conformational convergence or independent replicas.','Log warnings require case review; neither accepting nor discarding a biological label follows automatically.','Inputrec inspection files are not reusable grompp inputs.']};atomic_json(output/'summary.json',summary);print(json.dumps({k:v for k,v in summary.items() if k not in ['source_sha256','limitations']},indent=2))

if __name__=='__main__':
    p=__import__('argparse').ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
