#!/usr/bin/env python3
"""Receipt-compatible entry point for the repaired existing pRosettaC baseline.

Keeps all raw engine outputs. A separate retained table contains only poses that
pass molecular identity, overlap and fixed-head checks. Native-frame mode is an
integration control, not an independent recovery benchmark.
"""
import argparse,json,os,shutil,subprocess,sys,time,importlib.metadata as metadata
import platform
from pathlib import Path
BASE=Path(__file__).resolve().parent;ROOT=BASE.parent
sys.path.insert(0,str(ROOT/'scripts/lib'))
from protac_qc import ValidationError,atomic_json,sha256,parse_scores,engine_source_hashes
from audit_outputs import audit
from constrained_linker import DEFAULT_EMBED_ITERATIONS,DEFAULT_OUTER_ATTEMPTS

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',type=Path,required=True);p.add_argument('--seed',type=int,required=True)
    p.add_argument('--protocol',choices=['global','native_frame'],default='global')
    p.add_argument('--n-dock',type=int,default=200);p.add_argument('--n-local',type=int,default=50)
    p.add_argument('--conformers',type=int,default=100);p.add_argument('--workers',type=int,default=2)
    p.add_argument('--max-global-trials',type=int,default=10000000);a=p.parse_args()
    run=a.run_dir.resolve();raw=run/'engine_raw';output=run/'Patchdock_Results'
    if raw.exists() or output.exists():raise ValidationError('Engine attempt already used')
    if min(a.n_dock,a.n_local,a.conformers,a.workers,a.max_global_trials)<1 or a.n_local>50 or a.conformers>100:
        raise ValidationError('Invalid baseline sampling budget')
    report={'schema_version':2,'status':'RUNNING','score_role':'raw','seed':a.seed,
        'engine_kind':'existing_pRosettaC_centered_fixed_head_ETKDG_no_basic_knowledge','protocol':a.protocol,
        'linker_method':'CENTERED_TORSION_DG_FIXED_HEADS_UFF',
        'docking_microstate':'supplier_neutral_reference','started_unix':time.time(),'stages':{},
        'engine_sha256':engine_source_hashes(ROOT),
        'settings':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}}
    report['settings'].update(linker_embedding_outer_attempts=DEFAULT_OUTER_ATTEMPTS,
                              linker_embedding_internal_iterations=DEFAULT_EMBED_ITERATIONS)
    report['engine_sha256']['scripts/lib/protac_qc.py']=sha256(ROOT/'scripts/lib/protac_qc.py')
    atomic_json(run/'engine_status.json',report)
    try:
        required={'numpy':'2.2.6','scipy':'1.15.2','rdkit':'2025.9.5','pyrosetta':'2026.3+releasequarterly.5e498f1409'}
        observed={name:metadata.version(name) for name in required}
        report['environment']={'python':sys.version,'versions':observed}
        from openbabel import openbabel
        report['environment'].update(openbabel_library_version=openbabel.OBReleaseVersion(),
            scikit_learn_version=metadata.version('scikit-learn'),platform=platform.platform(),
            available_cpu_affinity=len(os.sched_getaffinity(0)))
        if sys.version_info[:2]!=(3,10) or observed!=required:
            raise ValidationError('Use the validated Python 3.10 protac environment; package versions differ')
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        command=[sys.executable,str(BASE/'execute.py'),'--input',str(run),'--run',str(raw),
            '--seed',str(a.seed),'--protocol',a.protocol,'--n-dock',str(a.n_dock),
            '--local-poses',str(a.n_local),'--conformers',str(a.conformers),'--workers',str(a.workers),
            '--max-global-trials',str(a.max_global_trials)]
        subprocess.run(command,check=True,env=env)
        if engine_source_hashes(ROOT)!=report['engine_sha256']:
            raise ValidationError('Execution snapshot changed during the run')
        execution=json.loads((raw/'diagnostic.json').read_text())
        local=list((raw/'Patchdock_Results').glob('pd.*_docking_????.pdb'))
        sampling=list((raw/'Patchdock_Results').glob('confs_*.sdf.sampling.json'))
        if len(local)!=len(sampling):raise ValidationError('Missing linker accounting for a local pose')
        feasible=sum(json.loads(f.read_text())['generated_conformers']>0 for f in sampling)
        global_file=raw/'global_sampling.json'
        if global_file.exists():global_counts=json.loads(global_file.read_text())
        elif a.protocol=='native_frame':global_counts={'attempted':1,'succeeded':1,'rejected':0,'failed':0}
        else:raise ValidationError('Missing global sampling ledger')
        if len(local)!=global_counts['succeeded']*a.n_local:raise ValidationError('Local sampling output count mismatch')
        quality=audit(raw)
        rows=parse_scores(raw/'Patchdock_Results/score.sc',allow_empty=feasible==0)
        if len(rows)!=feasible:raise ValidationError('Feasible-linker/scoring count mismatch')
        accepted=[]
        for item in quality['poses']:
            if item['status']!='MECHANICAL_QC_PASSED':continue
            binding=item.get('head_engagement',{})
            if binding.get('status')!='MEASURED_NOT_A_BINDING_LABEL':raise ValidationError('Missing mapped bound-head reference')
            if max(binding['vhl_core_rmsd_A'],binding['poi_core_rmsd_A'])>0.05:
                raise ValidationError('Fixed-head preservation contract violated after scoring')
            accepted.append(item)
        output.mkdir()
        with (output/'score.sc').open('w') as f:
            f.write('SEQUENCE: \nSCORE: total_score description\n')
            for row in accepted:
                name=row['description'];shutil.copy2(raw/'Patchdock_Results'/(name+'.pdb'),output/(name+'.pdb'))
                f.write(f'SCORE: {row["total_score"]:.6f} {name}\n')
        n=len(local);scored=len(rows);kept=len(accepted)
        report['stages']={'global_sampling':global_counts,
            'local_docking':{'attempted':n,'succeeded':n,'rejected':0,'failed':0},
            'linker_modeling':{'attempted':n,'succeeded':scored,'rejected':n-scored,'failed':0},
            'scoring':{'attempted':scored,'succeeded':scored,'rejected':0,'failed':0},
            'pose_qc':{'attempted':scored,'succeeded':kept,'rejected':scored-kept,'failed':0}}
        report.update(status='COMPLETE' if accepted else 'COMPLETE_NO_VALID_POSES',
            raw_score_sha256=sha256(raw/'Patchdock_Results/score.sc'),
            raw_output_qc_sha256=sha256(raw/'baseline_output_qc.json'),
            sampling_sha256={f.name:sha256(f) for f in sampling},
            interpretation='Counts describe this sampling budget; no biological degrader/nondegrader label is inferred')
    except BaseException as exc:
        report.update(status='FAILED',error=repr(exc));raise
    finally:
        report['finished_unix']=time.time();atomic_json(run/'engine_status.json',report)

if __name__=='__main__':main()
