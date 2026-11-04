#!/usr/bin/env python3
"""Plan or execute immutable baseline-workflow inputs using the LOCAL backend.

--plan is read-only. --run --pilot is bounded; full --run needs a versioned scientific release.
"""
import argparse
import csv
import fcntl
import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import sha256, ValidationError

def pair_hash(folder):
    # Hash config and every local dependency, including the original symlink target content.
    import hashlib
    files=sorted(p for p in folder.iterdir() if p.is_file())
    return hashlib.sha256(json.dumps([(p.name,sha256(p)) for p in files]).encode()).hexdigest()

def get_all_inputs(input_dir,protac_filter=None,gene_filter=None):
    pairs=[]
    for d in sorted(Path(input_dir).glob('*/*')):
        if not d.is_dir(): continue
        gene,_,uid=d.name.rpartition('_')
        if (d/'preparation.json').exists():
            info=json.loads((d/'preparation.json').read_text())
            gene,uid=info.get('gene') or gene,info.get('uniprot_id') or uid
        if protac_filter and d.parent.name!=protac_filter: continue
        if gene_filter and gene!=gene_filter: continue
        pairs.append({'protac':d.parent.name,'gene':gene,'uid':uid,'run_name':d.name,'local_dir':str(d.resolve()),
                      'input_hash':pair_hash(d)})
    return pairs

def build_manifest(manifest,pairs):
    """Locked atomic union; filtered plans never delete unrelated manifest rows."""
    manifest=Path(manifest); manifest.parent.mkdir(parents=True,exist_ok=True)
    with (manifest.parent/(manifest.name+'.lock')).open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        existing=list(csv.DictReader(manifest.open(),delimiter='\t')) if manifest.exists() else []
        key=lambda r:(r['protac'],r.get('run_name') or r['gene']+'_'+r['uid'])
        merged={key(r):r for r in existing}
        if len(merged)!=len(existing): raise ValidationError('Duplicate keys in existing manifest')
        for pair in pairs:
            old=merged.get(key(pair),{})
            changed=old.get('input_hash')!=pair['input_hash']
            if changed and old.get('status') in {'RUNNING','SUBMITTED','SUBMITTING'}:
                raise ValidationError('Inputs changed during active attempt; inspect before continuing')
            state=dict(old)
            state.update(pair)
            if not old or changed:
                state.update(status='BLOCKED',attempt_id='',scheduler_id='',reason='PREFLIGHT_AND_ENGINE_RELEASE_REQUIRED')
            merged[key(pair)]=state
        rows=[merged[k] for k in sorted(merged)]
        fields=sorted({k for r in rows for k in r})
        fd,name=tempfile.mkstemp(dir=manifest.parent,prefix='.manifest-')
        try:
            with os.fdopen(fd,'w') as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
                f.flush(); os.fsync(f.fileno())
            os.replace(name,manifest)
        finally:
            if os.path.exists(name): os.unlink(name)
        return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True)
    for mode in ['plan','build-manifest','status','run']: g.add_argument('--'+mode,action='store_true')
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--inputs',type=Path); p.add_argument('--protac'); p.add_argument('--gene')
    p.add_argument('--python',type=Path);p.add_argument('--seed',type=int)
    p.add_argument('--workers',type=int);p.add_argument('--n-dock',type=int)
    p.add_argument('--n-local',type=int);p.add_argument('--conformers',type=int)
    p.add_argument('--max-global-trials',type=int);p.add_argument('--timeout-seconds',type=float)
    p.add_argument('--runs',type=Path);p.add_argument('--pilot',action='store_true')
    p.add_argument('--protocol',choices=['global','native_frame'],default='global')
    p.add_argument('--max-pairs',type=int);p.add_argument('--rerun',action='store_true')
    a=p.parse_args()
    defaults={'workers':1 if a.pilot else 2,'n_dock':2 if a.pilot else 200,
              'n_local':3 if a.pilot else 50,'conformers':2 if a.pilot else 100,
              'max_global_trials':10000 if a.pilot else 10000000,
              'timeout_seconds':600 if a.pilot else 86400}
    for key,value in defaults.items():
        if getattr(a,key) is None:setattr(a,key,value)
    if a.runs is None:a.runs=a.root/'runs'
    manifest=a.runs/'manifest.tsv'
    if a.run:
        from campaign import execute
        pairs=get_all_inputs(a.inputs or a.root/('prepared_inputs' if a.pilot else 'prepared_screen'),a.protac,a.gene)
        limit=a.max_pairs if a.max_pairs is not None else (1 if a.pilot else len(pairs))
        if limit<1:raise ValidationError('max-pairs must be positive')
        execute(a.root,pairs[:limit],a)
        return
    if a.status:
        if not manifest.exists(): print('No manifest.'); return
        rows=list(csv.DictReader(manifest.open(),delimiter='\t'))
    else:
        rows=get_all_inputs(a.inputs or a.root/('prepared_inputs' if a.pilot else 'prepared_screen'),a.protac,a.gene)
        if a.max_pairs is not None or a.pilot:rows=rows[:a.max_pairs if a.max_pairs is not None else 1]
        if a.build_manifest: rows=build_manifest(manifest,rows)
    from collections import Counter
    print(json.dumps({'pairs':len(rows),'by_protac':dict(Counter(r['protac'] for r in rows)),
                      'status':dict(Counter(r.get('status','UNVALIDATED') for r in rows)),
                      'mode':'BOUNDED_PILOT' if a.pilot else 'CAMPAIGN_REQUIRES_VERSIONED_RELEASE',
                      'budget':defaults | {k:getattr(a,k) for k in defaults},'retention':'FULL_RAW_OUTPUTS',
                      'execution_started':False},indent=2))

if __name__=='__main__': main()
