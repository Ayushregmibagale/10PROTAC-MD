"""Foreground LOCAL controller for the repaired original baseline.

Inputs, attempts and transitions are durable. Pilot mode has enforced small
budgets; campaign mode requires a versioned release receipt and storage budget.
No SSH, detached submission, automatic deletion or experiment-engine invocation.
"""
import argparse,csv,fcntl,hashlib,importlib.util,json,os,shutil,sys,tempfile,time,uuid,subprocess
from pathlib import Path
from protac_qc import ValidationError,sha256,verify_receipt,engine_source_hashes,atomic_json

def module(root,name,relative):
    spec=importlib.util.spec_from_file_location(name,root/relative)
    loaded=importlib.util.module_from_spec(spec);spec.loader.exec_module(loaded);return loaded

def source_hashes(root):
    return engine_source_hashes(root)


def snapshot_sources(root,dest,expected):
    """Freeze executable code/config before starting a worker or recorder."""
    dest.mkdir(parents=True,exist_ok=False)
    shutil.copytree(root/'baseline_engine',dest/'baseline_engine',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for relative in ['scripts/lib/protac_qc.py','scripts/lib/campaign.py',
                     'scripts/01_setup/02_batch_prosettac.py','scripts/01_setup/03_record_attempt.py']:
        target=dest/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/relative,target)
    if engine_source_hashes(dest)!=expected:
        raise ValidationError('Sources changed while creating execution snapshot; preserve and retry')
    atomic_json(dest/'source_manifest.json',{'status':'VERIFIED_EXECUTION_SNAPSHOT','source_sha256':expected})
    return dest

def transition(manifest,pair,attempt,status,key,reason=''):
    manifest.parent.mkdir(parents=True,exist_ok=True)
    with manifest.with_suffix('.tsv.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        rows=list(csv.DictReader(manifest.open(),delimiter='\t')) if manifest.exists() else []
        matches=[row for row in rows if (row['protac'],row['run_name'])==(pair['protac'],pair['run_name'])]
        if len(matches)>1:raise ValidationError('Duplicate manifest identities')
        old=matches[0] if matches else {}
        if status=='PREPARING' and old.get('status') in {'RUNNING','PREPARING','SUBMITTED','SUBMITTING'}:
            raise ValidationError('An active/stale attempt requires reconciliation before retry: '+old.get('local_attempt_dir',''))
        row=dict(old,**pair);row.update(status=status,attempt_id=attempt.name,local_attempt_dir=str(attempt),
            attempt_key=key,reason=reason,last_updated_unix=str(time.time()),scheduler_id='LOCAL_FOREGROUND')
        if matches:rows[rows.index(old)]=row
        else:rows.append(row)
        fields=sorted({k for r in rows for k in r});fd,name=tempfile.mkstemp(dir=manifest.parent,prefix='.manifest-')
        try:
            with os.fdopen(fd,'w') as out:
                writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows);out.flush();os.fsync(out.fileno())
            os.replace(name,manifest)
        finally:
            if os.path.exists(name):os.unlink(name)
        with (manifest.parent/'attempt_events.jsonl').open('a') as out:out.write(json.dumps(row)+'\n')

def validate_plan(root,pairs,args):
    if not pairs:raise ValidationError('No selected inputs')
    if args.protocol=='native_frame' and not args.pilot:raise ValidationError('Native-frame controls are restricted to explicit pilot mode')
    if args.pilot and (len(pairs)>5 or args.n_dock>5 or args.n_local>5 or args.conformers>10 or args.workers>2 or args.max_global_trials>20000 or args.timeout_seconds>3600):
        raise ValidationError('Pilot exceeds bounds: 5 pairs/global/local, 10 conformers, 2 workers, 20,000 trials or one hour')
    if min(args.n_dock,args.n_local,args.conformers,args.workers,args.max_global_trials)<1 or args.timeout_seconds<=0:
        raise ValidationError('Resource/sampling budgets must be positive')
    if args.workers>len(os.sched_getaffinity(0)):raise ValidationError('Workers exceed available CPU allocation')
    preflight=module(root,'campaign_preflight','scripts/01_setup/01b_preflight.py')
    registry=json.loads((root/'config/molecules.json').read_text())['compounds']
    for pair in pairs:
        result=preflight.validate_pair(Path(pair['local_dir']),registry)
        invalid=[issue for issue in result['issues'] if issue not in {'BENCHMARK_RELEASE_PENDING','MODELED_HEAD_REVIEW_PENDING'}]
        if invalid:raise ValidationError('Invalid input '+pair['run_name']+': '+str(invalid))
    if args.pilot:return {'mode':'BOUNDED_PILOT','minimum_free_bytes':5_000_000_000}
    path=root/'audits/baseline_release.json'
    if not path.exists():raise ValidationError('No baseline campaign release; use bounded pilot mode for validation')
    release=json.loads(path.read_text())
    if release.get('status')!='RELEASED' or release.get('engine_sha256')!=source_hashes(root):raise ValidationError('Baseline release is absent, changed or obsolete')
    if not release.get('evidence_sha256'):raise ValidationError('Scientific release evidence is missing')
    for relative,digest in release['evidence_sha256'].items():
        file=(root/relative).resolve()
        if not file.is_relative_to(root.resolve()) or sha256(file)!=digest:raise ValidationError('Release evidence changed')
    for pair in pairs:
        if release['input_sha256'].get(pair['protac']+'/'+pair['run_name'])!=pair['input_hash']:
            raise ValidationError('Input is not part of this release: '+pair['run_name'])
        preparation=json.loads((Path(pair['local_dir'])/'preparation.json').read_text())
        if preparation.get('contains_modeled_head') and preparation['head_bundle_sha256'] not in release.get('reviewed_modeled_head_bundles',[]):
            raise ValidationError('Modeled-head evidence is outside this release')
    for field in ['n_dock','n_local','conformers','workers','max_global_trials','timeout_seconds']:
        if getattr(args,field)>release['budget'][field]:raise ValidationError('Budget exceeds release: '+field)
    return {'mode':'RELEASED_CAMPAIGN','minimum_free_bytes':int(release['budget']['bytes_per_pair'])*len(pairs),
            'release_sha256':sha256(path),'retention':'full_raw_outputs'}

def execute(root,pairs,args):
    root=root.resolve();policy=validate_plan(root,pairs,args)
    if args.python is None or args.seed is None:raise ValidationError('Explicit --python and --seed are required')
    runs=args.runs.resolve();runs.mkdir(parents=True,exist_ok=True);manifest=runs/'manifest.tsv'
    protocol={k:getattr(args,k) for k in ['seed','protocol','n_dock','n_local','conformers','workers','max_global_trials','timeout_seconds']}
    protocol['engine_sha256']=source_hashes(root)
    with (runs/'campaign.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if shutil.disk_usage(runs).free<policy['minimum_free_bytes']:raise ValidationError('Insufficient free space for selected full-retention plan')
        for pair in pairs:
            key=hashlib.sha256(json.dumps({'input':pair['input_hash'],'protocol':protocol},sort_keys=True).encode()).hexdigest()
            previous=list(csv.DictReader(manifest.open(),delimiter='\t')) if manifest.exists() else []
            old=next((r for r in previous if (r['protac'],r['run_name'])==(pair['protac'],pair['run_name'])),{})
            if old.get('attempt_key')==key and old.get('status') in {'COMPLETE','COMPLETE_NO_VALID_POSES'} and not args.rerun:
                old_run=Path(old['local_attempt_dir']);verify_receipt(old_run,json.loads((old_run/'selection.json').read_text()))
                print('Already validated: '+pair['protac']+'/'+pair['run_name'],flush=True);continue
            if shutil.disk_usage(runs).free<5_000_000_000:raise ValidationError('Free-space floor reached; campaign stopped before next attempt')
            attempt=runs/'attempts'/pair['protac']/pair['run_name']/str(uuid.uuid4())
            transition(manifest,pair,attempt,'PREPARING',key)
            terminal_recorded=False
            try:
                attempt.mkdir(parents=True,exist_ok=False)
                for file in Path(pair['local_dir']).iterdir():
                    if file.is_file():shutil.copy2(file,attempt/file.name)
                snapshot=snapshot_sources(root,attempt/'source_snapshot',protocol['engine_sha256'])
                seed=1+int.from_bytes(hashlib.sha256((str(args.seed)+'|'+pair['protac']+'/'+pair['run_name']).encode()).digest()[:4],'big')%2147483646
                atomic_json(attempt/'campaign_intent.json',{'pair':pair,'policy':policy,'protocol':protocol,'derived_seed':seed,
                    'source_snapshot_manifest_sha256':sha256(snapshot/'source_manifest.json')})
                command=[str(args.python),str(snapshot/'baseline_engine/run.py'),'--run-dir',str(attempt),'--seed',str(seed),'--protocol',args.protocol,
                    '--n-dock',str(args.n_dock),'--n-local',str(args.n_local),'--conformers',str(args.conformers),
                    '--workers',str(args.workers),'--max-global-trials',str(args.max_global_trials)]
                transition(manifest,pair,attempt,'RUNNING',key)
                recorded=subprocess.run([str(args.python),str(snapshot/'scripts/01_setup/03_record_attempt.py'),
                    '--run-dir',str(attempt),'--timeout-seconds',str(args.timeout_seconds),'--',*command],
                    text=True,capture_output=True)
                (attempt/'recorder.stdout.log').write_text(recorded.stdout)
                (attempt/'recorder.stderr.log').write_text(recorded.stderr)
                result=json.loads((attempt/'completion.json').read_text())
                transition(manifest,pair,attempt,result['status'],key,result.get('reason',''))
                terminal_recorded=True
                print(json.dumps({'attempt':str(attempt),'status':result['status']}),flush=True)
                if result['status'] not in {'COMPLETE','COMPLETE_NO_VALID_POSES'}:raise ValidationError('Stopped after failed attempt')
            except BaseException as exc:
                if not terminal_recorded:transition(manifest,pair,attempt,'INTERRUPTED' if isinstance(exc,KeyboardInterrupt) else 'FAILED',key,str(exc))
                raise
