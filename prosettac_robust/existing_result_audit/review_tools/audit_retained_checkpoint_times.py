"""Read final checkpoint times without simulating or changing original files."""
from pathlib import Path
import collections,concurrent.futures,json,re,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import ValidationError,sha256,atomic_json
PROJECT=Path('/home/ayushregmibagale/protac-md');GMX=Path('/home/ayushregmibagale/.conda/envs/gromacs/bin/gmx')


def inspect(case,output):
    records=[]
    for stage,name,expected in [('nvt','nvt',200.),('npt_restrained','npt_restrained',500.),('npt','npt',500.),('prod','prod',20000.)]:
        path=PROJECT/'random_screen'/case['poi']/case['pose_id']/stage/(name+'.cpt');row={'poi':case['poi'],'pose_id':case['pose_id'],'stage':stage,'checkpoint_available':path.exists(),'nominal_requested_time_ps':expected}
        if path.exists():
            before=sha256(path);command=[str(GMX),'check','-f',str(path)];result=subprocess.run(command,capture_output=True,text=True,timeout=60);text=result.stdout+result.stderr;target=output/'logs'/case['poi']/case['pose_id'];target.mkdir(parents=True,exist_ok=True);log=target/(stage+'_checkpoint.txt');log.write_text(text)
            times=re.findall(r'Last frame\s+-?\d+\s+time\s+([-+\d.eE]+)',text);atoms=re.findall(r'# Atoms\s+(\d+)',text)
            if result.returncode or len(times)!=1 or len(atoms)!=1:raise ValidationError('Checkpoint read failed: '+str(log))
            row.update(checkpoint=str(path),checkpoint_sha256=before,checkpoint_time_ps=float(times[0]),checkpoint_atoms=int(atoms[0]),reaches_requested_time=abs(float(times[0])-expected)<.001,command=command,log_sha256=sha256(log))
            if sha256(path)!=before:raise ValidationError('Checkpoint changed during read')
        records.append(row)
    return records


def main(output):
    if output.exists():raise ValidationError('Fresh checkpoint audit output required')
    roster=ROOT/'sjf8240_review/md_graph_pose_lineage.json';hashes={str(p):sha256(p) for p in [roster,GMX,Path(__file__)]};cases=json.loads(roster.read_text())
    if len(cases)!=174:raise ValidationError('Expected 174 retained production cases')
    output.mkdir(parents=True);rows=[];failures=[];start=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(inspect,c,output):c for c in cases}
        for future in concurrent.futures.as_completed(futures):
            try:rows.extend(future.result())
            except Exception as exc:failures.append(dict(futures[future],error=str(exc)))
    if any(sha256(Path(p))!=h for p,h in hashes.items()):raise ValidationError('Audit source changed')
    rows.sort(key=lambda r:(r['poi'],r['pose_id'],r['stage']));atomic_json(output/'per_stage.json',rows)
    summary={'status':'COMPLETE_RETAINED_CHECKPOINT_TIME_AUDIT' if not failures else 'INCOMPLETE_CHECKPOINT_TIME_AUDIT','cases_requested':174,'stages':{stage:{'checkpoint_available':sum(r['checkpoint_available'] for r in rows if r['stage']==stage),'checkpoint_reaches_requested_time':sum(r.get('reaches_requested_time',False) for r in rows if r['stage']==stage),'observed_time_ps_counts':dict(collections.Counter(str(r['checkpoint_time_ps']) for r in rows if r['stage']==stage and r['checkpoint_available']))} for stage in ['nvt','npt_restrained','npt','prod']},'failures':failures,'elapsed_s':time.monotonic()-start,'source_sha256':hashes,'records_sha256':sha256(output/'per_stage.json'),'original_files_modified':False,'new_MD_executed':False,'limitations':['Checkpoint time/parse evidence is not a conformational-convergence certificate or independent biological replicate.','File paths and hashes describe retained artifacts; no missing historical generation receipt is fabricated.','Absent intermediate checkpoint files do not prove a stage was omitted.']};atomic_json(output/'summary.json',summary);print(json.dumps({k:v for k,v in summary.items() if k not in ['source_sha256','limitations']},indent=2))

if __name__=='__main__':
    p=__import__('argparse').ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
