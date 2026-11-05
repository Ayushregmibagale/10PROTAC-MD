#!/usr/bin/env python3
"""Receipt-bound thermodynamic diagnostics, never a convergence certificate."""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,atomic_json,sha256


def summarize(time,values,discard,block):
    time=np.asarray(time,dtype=float);values=np.asarray(values,dtype=float)
    if not np.isfinite(discard) or not np.isfinite(block) or discard<0 or block<=0:
        raise ValidationError('Explicit finite discard/block times are required')
    if time.ndim!=1 or values.ndim!=2 or len(time)!=len(values) or len(time)<2:
        raise ValidationError('Invalid thermodynamic time series')
    if not np.isfinite(time).all() or not np.isfinite(values).all() or not (np.diff(time)>0).all():
        raise ValidationError('Nonfinite or non-increasing thermodynamic time series')
    keep=time>=discard
    if keep.sum()<2:raise ValidationError('Insufficient retained thermodynamic samples')
    t=time[keep];v=values[keep];blocks=[]
    for i in range(int(np.floor((t[-1]-discard)/block))):
        low=discard+i*block;high=low+block
        mask=(t>=low)&((t<=high) if high==t[-1] else (t<high))
        if mask.sum()<2:raise ValidationError('A complete time block has too few samples')
        blocks.append({'start_ps':low,'end_ps':high,'samples':int(mask.sum()),'means':v[mask].mean(0).tolist()})
    return {'retained_samples':int(keep.sum()),'retained_time_ps':[float(t[0]),float(t[-1])],
        'sample_means':v.mean(0).tolist(),'sample_standard_deviations':v.std(0,ddof=1).tolist(),
        'linear_slopes_per_ps':np.polyfit(t-t[0],v,1)[0].tolist(),
        'block_means':blocks,'independent_replicates':1,'convergence_certified':False,
        'interpretation':'Descriptive diagnostics; time samples/blocks are not independent simulation replicas and no confidence interval is inferred'}


def analyze(work,stage,gmx,discard,block,output):
    work=work.resolve();receipt_path=work/(stage+'.receipt.json')
    receipt=json.loads(receipt_path.read_text());edr=work/(stage+'.edr')
    if receipt['status']!='COMPLETE' or receipt['outputs'].get(edr.name)!=sha256(edr):
        raise ValidationError('Energy file is not bound to a completed stage')
    xvg=work/(stage+'_thermodynamics_checked.xvg')
    command=[str(gmx),'energy','-f',str(edr),'-o',str(xvg)]
    requested=['Temperature','Pressure','Potential']+([] if stage=='nvt' else ['Density'])
    result=subprocess.run(command,input='\n'.join(requested+['0'])+'\n',text=True,capture_output=True)
    log=work/(stage+'_thermodynamics_checked.log');log.write_text(result.stdout+result.stderr)
    if result.returncode:raise ValidationError('Gromacs energy extraction failed')
    text=xvg.read_text();legends={int(i):name for i,name in re.findall(r'@\s+s(\d+)\s+legend\s+"([^"]+)"',text)}
    names=[legends[i] for i in range(len(legends))]
    if set(names)!=set(requested):raise ValidationError('Extracted energy fields differ from requested fields')
    data=np.loadtxt([line for line in text.splitlines() if line and not line.startswith(('#','@'))],ndmin=2)
    if data.shape[1]!=len(names)+1:raise ValidationError('Energy field count mismatch')
    summary=summarize(data[:,0],data[:,1:],discard,block)
    report={'status':'THERMODYNAMIC_DIAGNOSTIC_ONLY','stage':stage,'discard_ps':discard,'block_ps':block,
        'fields':names,'units':{'Temperature':'K','Pressure':'bar','Potential':'kJ/mol','Density':'kg/m^3'},
        'summary':summary,'command':command,'stdin_fields':requested,
        'receipt_sha256':sha256(receipt_path),'edr_sha256':sha256(edr),'gmx_sha256':sha256(gmx),
        'extracted_series_sha256':sha256(xvg),'extraction_log_sha256':sha256(log),
        'limitations':['Chosen time window is a post-hoc validation diagnostic, not a prespecified production analysis',
            'Thermodynamic stability does not establish conformational convergence or binding stability',
            'Statistics use the extracted recorded samples; Gromacs energy-accumulator statistics can differ']}
    atomic_json(output,report);print(json.dumps({'stage':stage,'fields':names,'means':summary['sample_means'],
        'blocks':len(summary['block_means']),'convergence_certified':False},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['work','gmx','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--stage',choices=['nvt','npt_restrained','npt','prod'],required=True)
    p.add_argument('--discard-ps',type=float,required=True);p.add_argument('--block-ps',type=float,required=True)
    a=p.parse_args();analyze(a.work,a.stage,a.gmx,a.discard_ps,a.block_ps,a.output)
