#!/usr/bin/env python3
"""Run one explicit MD stage; separate bounded smoke runs from scientific stages.

Each stage has immutable inputs and hashed completion receipts. Interrupted
stages resume with GROMACS checkpoint checksum validation and append semantics.
"""
import argparse
import json
import os
import math
import re
import fcntl
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,atomic_json,sha256

STAGES={'em':(None,'system.gro'),'nvt':('em','em.gro'),
        'npt_restrained':('nvt','nvt.gro'),'npt':('npt_restrained','npt_restrained.gro'),
        'prod':('npt','npt.gro')}

def validate_stage_log(log,stage,smoke_steps,mdp_text):
    if 'Finished mdrun' not in log:raise ValidationError('No GROMACS completion marker')
    if any(x.lower() in log.lower() for x in ['LINCS WARNING','Constraint error','not finite']):
        raise ValidationError('Dynamics instability found; stage cannot be accepted')
    settings={}
    for line in mdp_text.splitlines():
        line=line.split(';',1)[0]
        if '=' in line:
            key,value=line.split('=',1);settings[key.strip().lower().replace('_','-')]=value.strip()
    requested=int(settings['nsteps']);initial=int(settings.get('init-step','0'))
    if requested<0:raise ValidationError('A finite requested stage length is required')
    steps=re.findall(r'^\s*Step\s+Time\s*\n\s*(\d+)\s+[-+\d.eE]+',log,re.MULTILINE)
    if not steps:raise ValidationError('No actual GROMACS step count in stage log')
    final=int(steps[-1]);converged=False
    if stage=='em':
        # A negative convergence statement must never satisfy this test.
        converged=bool(re.search(r'^\s*(?:Steepest Descents|Conjugate Gradients|L-BFGS) converged to Fmax <',log,re.MULTILINE))
        force=re.search(r'^Maximum force\s*=\s*(\S+)',log,re.MULTILINE)
        if not force or not math.isfinite(float(force[1])):raise ValidationError('Missing/nonfinite minimized maximum force')
        threshold=float(settings.get('emtol','10'))
        if not smoke_steps and (not converged or float(force[1])>=threshold):
            raise ValidationError('Minimization did not meet its force convergence criterion')
    if final<initial+requested and not (stage=='em' and converged):
        raise ValidationError(f'Stage stopped before requested final step: {final} < {initial+requested}; checkpoint may be resumed')
    return {'last_logged_step':final,'requested_final_step':initial+requested,'force_converged':converged if stage=='em' else None}


def verify_tpr_binding(work,prefix,payload):
    tpr=work/(prefix+'.tpr');binding=work/(prefix+'.tpr_generation.json')
    if not binding.exists():raise ValidationError('Existing TPR lacks a generation receipt; preserve this directory and prepare a fresh system')
    recorded=json.loads(binding.read_text())
    if recorded['intent']!=payload or sha256(tpr)!=recorded['tpr_sha256']:
        raise ValidationError('TPR differs from its input/protocol generation receipt')


def verify_preparation(work,stage):
    meta=json.loads((work/'md_preparation.json').read_text())
    if meta.get('status')!='TOPOLOGY_PREPARED':raise ValidationError('System topology has not passed preparation')
    ff=work/'amber99sb-ildn.ff'
    if not (ff/'forcefield.itp').is_file():raise ValidationError('Portable, hash-bound force field is required; rebuild the older system')
    files=[work/'topol.top',work/'index.ndx',work/'system.gro',work/'mdp'/(stage+'.mdp')]
    files+=list(work.glob('*.itp'))+[p for p in ff.rglob('*') if p.is_file()]
    for path in files:
        relative=str(path.relative_to(work))
        if meta.get('files',{}).get(relative)!=sha256(path):
            raise ValidationError('MD preparation file missing from receipt or altered: '+relative)
    return meta

def run(work,stage,gmx,threads,smoke_steps=None,device='auto',gpu_id=None):
    # All stages in one system share outputs; independent replicas use separate systems.
    with (work/'.md_stage.lock').open('w') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValidationError('Another MD stage is using this system directory')
        return _run(work,stage,gmx,threads,smoke_steps,device,gpu_id)


def _run(work,stage,gmx,threads,smoke_steps=None,device='auto',gpu_id=None):
    work=work.resolve(); gmx=gmx.resolve()
    verify_preparation(work,stage)
    if threads<1:raise ValidationError('threads must be positive')
    if smoke_steps is not None and (smoke_steps<1 or stage not in {'em','nvt'}):
        raise ValidationError('Smoke runs require positive steps and stage em or nvt')
    previous,coordinates=STAGES[stage]
    prefix=('smoke_' if smoke_steps else '')+stage
    receipt=work/(prefix+'.receipt.json');intent=work/(prefix+'.intent.json')
    if previous and not smoke_steps:
        parent=json.loads((work/(previous+'.receipt.json')).read_text())
        if parent['status']!='COMPLETE':raise ValidationError('Preceding stage incomplete')
        for name,digest in parent['outputs'].items():
            if sha256(work/name)!=digest:raise ValidationError('Preceding stage output changed: '+name)
        validate_stage_log((work/(previous+'.log')).read_text(),previous,None,(work/'mdp'/(previous+'.mdp')).read_text())
    if smoke_steps and stage=='nvt':
        coordinates='smoke_em.gro'
        parent=json.loads((work/'smoke_em.receipt.json').read_text())
        if parent['status']!='SMOKE_COMPLETE':raise ValidationError('Smoke EM is incomplete')
        for name,digest in parent['outputs'].items():
            if sha256(work/name)!=digest:raise ValidationError('Smoke EM output changed: '+name)
        validate_stage_log((work/'smoke_em.log').read_text(),'em',parent['intent']['smoke_steps'],(work/'smoke_em.mdp').read_text())
    mdp=work/'mdp'/(stage+'.mdp')
    if smoke_steps:
        text=mdp.read_text(); lines=[]
        for line in text.splitlines():
            key=line.split('=')[0].strip()
            if key=='nsteps':line='nsteps = '+str(smoke_steps)
            if key in {'nstenergy','nstlog','nstxout-compressed'}:line=key+' = 1'
            lines.append(line)
        mdp=work/(prefix+'.mdp');expected='\n'.join(lines)+'\n'
        if mdp.exists() and mdp.read_text()!=expected:raise ValidationError('Smoke protocol changed')
        mdp.write_text(expected)
    inputs=[mdp,work/coordinates,work/'topol.top',work/'index.ndx',work/'md_preparation.json']
    inputs+=sorted(work.glob('*.itp'))
    inputs+=sorted(p for p in (work/'amber99sb-ildn.ff').rglob('*') if p.is_file())
    if previous and stage!='nvt' and not smoke_steps:inputs.append(work/(previous+'.cpt'))
    if device not in {'auto','cpu','gpu'}:raise ValidationError('Unknown GROMACS device mode')
    payload={'stage':stage,'smoke_steps':smoke_steps,'threads':threads,'device':device,'gpu_id':gpu_id,
             'gmx_sha256':sha256(gmx),'inputs':{str(p.relative_to(work)):sha256(p) for p in inputs}}
    if intent.exists() and json.loads(intent.read_text())!=payload:
        raise ValidationError('Stage inputs changed; use a new system directory')
    if receipt.exists():
        old=json.loads(receipt.read_text())
        if old['intent']!=payload:raise ValidationError('Completion receipt does not match inputs')
        for name,digest in old['outputs'].items():
            if sha256(work/name)!=digest:raise ValidationError('Completed output changed: '+name)
        validate_stage_log((work/(prefix+'.log')).read_text(),stage,smoke_steps,mdp.read_text())
        return old
    atomic_json(intent,payload)
    env=os.environ.copy();env['OMP_NUM_THREADS']=str(threads)
    def command(args,label):
        result=subprocess.run([str(gmx),*args],cwd=work,env=env,text=True,capture_output=True)
        with (work/(prefix+'.'+label+'.log')).open('a') as f:f.write(result.stdout+result.stderr)
        if result.returncode:raise ValidationError('GROMACS '+label+' failed; inspect stage log')
    if (work/(prefix+'.tpr')).exists():
        verify_tpr_binding(work,prefix,payload)
    else:
        if (work/(prefix+'.tpr_generation.json')).exists():raise ValidationError('Previously generated TPR is missing')
        args=['grompp','-f',str(mdp),'-c',coordinates,'-r',coordinates,'-p','topol.top',
              '-n','index.ndx','-o',prefix+'.tpr']
        if previous and stage!='nvt' and not smoke_steps:args+=['-t',previous+'.cpt']
        command(args,'grompp')
        atomic_json(work/(prefix+'.tpr_generation.json'),{'intent':payload,'tpr_sha256':sha256(work/(prefix+'.tpr')),
                    'grompp_command':[str(gmx),*args]})
    args=['mdrun','-deffnm',prefix,'-ntmpi','1','-ntomp',str(threads),'-nb',device,'-pme',device]
    if gpu_id is not None:args+=['-gpu_id',gpu_id]
    if (work/(prefix+'.cpt')).exists():args+=['-cpi',prefix+'.cpt','-append']
    command(args,'mdrun')
    log=(work/(prefix+'.log')).read_text()
    completion=validate_stage_log(log,stage,smoke_steps,mdp.read_text())
    outputs={p.name:sha256(p) for suffix in ['.gro','.cpt','.edr','.xtc','.log','.tpr']
             if (p:=work/(prefix+suffix)).exists()}
    if prefix+'.gro' not in outputs:raise ValidationError('Final coordinates absent')
    outputs[prefix+'.tpr_generation.json']=sha256(work/(prefix+'.tpr_generation.json'))
    result={'status':'SMOKE_COMPLETE' if smoke_steps else 'COMPLETE','intent':payload,'outputs':outputs,'completion_check':completion,
            'scientific_equilibration_complete':False if smoke_steps or stage!='prod' else None}
    atomic_json(receipt,result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--stage',choices=STAGES,required=True);p.add_argument('--gmx',type=Path,required=True)
    p.add_argument('--threads',type=int,default=2);p.add_argument('--smoke-steps',type=int)
    p.add_argument('--device',choices=['auto','cpu','gpu'],default='auto');p.add_argument('--gpu-id')
    a=p.parse_args();print(json.dumps(run(a.work,a.stage,a.gmx,a.threads,a.smoke_steps,a.device,a.gpu_id),indent=2))
