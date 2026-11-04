#!/usr/bin/env python3
"""Run one explicit MD stage; separate bounded smoke runs from scientific stages.

Each stage has immutable inputs and hashed completion receipts. Interrupted
stages resume with GROMACS checkpoint checksum validation and append semantics.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,atomic_json,sha256

STAGES={'em':(None,'system.gro'),'nvt':('em','em.gro'),
        'npt_restrained':('nvt','nvt.gro'),'npt':('npt_restrained','npt_restrained.gro'),
        'prod':('npt','npt.gro')}

def run(work,stage,gmx,threads,smoke_steps=None,device='auto',gpu_id=None):
    work=work.resolve(); gmx=gmx.resolve()
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
    if smoke_steps and stage=='nvt':coordinates='smoke_em.gro'
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
        return old
    atomic_json(intent,payload)
    env=os.environ.copy();env['OMP_NUM_THREADS']=str(threads)
    def command(args,label):
        result=subprocess.run([str(gmx),*args],cwd=work,env=env,text=True,capture_output=True)
        with (work/(prefix+'.'+label+'.log')).open('a') as f:f.write(result.stdout+result.stderr)
        if result.returncode:raise ValidationError('GROMACS '+label+' failed; inspect stage log')
    if not (work/(prefix+'.tpr')).exists():
        args=['grompp','-f',str(mdp),'-c',coordinates,'-r',coordinates,'-p','topol.top',
              '-n','index.ndx','-o',prefix+'.tpr']
        if previous and stage!='nvt' and not smoke_steps:args+=['-t',previous+'.cpt']
        command(args,'grompp')
    args=['mdrun','-deffnm',prefix,'-ntmpi','1','-ntomp',str(threads),'-nb',device,'-pme',device]
    if gpu_id is not None:args+=['-gpu_id',gpu_id]
    if (work/(prefix+'.cpt')).exists():args+=['-cpi',prefix+'.cpt','-append']
    command(args,'mdrun')
    log=(work/(prefix+'.log')).read_text()
    if 'Finished mdrun' not in log:raise ValidationError('No GROMACS completion marker')
    if stage=='em' and not smoke_steps and 'converged to Fmax <' not in log:
        raise ValidationError('Minimization reached its step limit without the force convergence criterion')
    if any(x in log for x in ['LINCS WARNING','Constraint error','not finite']):
        raise ValidationError('Dynamics instability found; stage cannot be accepted')
    outputs={p.name:sha256(p) for suffix in ['.gro','.cpt','.edr','.xtc','.log','.tpr']
             if (p:=work/(prefix+suffix)).exists()}
    if prefix+'.gro' not in outputs:raise ValidationError('Final coordinates absent')
    result={'status':'SMOKE_COMPLETE' if smoke_steps else 'COMPLETE','intent':payload,'outputs':outputs,
            'scientific_equilibration_complete':False if smoke_steps or stage!='prod' else None}
    atomic_json(receipt,result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--stage',choices=STAGES,required=True);p.add_argument('--gmx',type=Path,required=True)
    p.add_argument('--threads',type=int,default=2);p.add_argument('--smoke-steps',type=int)
    p.add_argument('--device',choices=['auto','cpu','gpu'],default='auto');p.add_argument('--gpu-id')
    a=p.parse_args();print(json.dumps(run(a.work,a.stage,a.gmx,a.threads,a.smoke_steps,a.device,a.gpu_id),indent=2))
