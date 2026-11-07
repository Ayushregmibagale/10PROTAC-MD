"""Summarize completed control stages; descriptive diagnostics, not convergence."""
from pathlib import Path
import json,subprocess,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,atomic_json
from gromacs_energy import named_energy_series
WORK=ROOT/'linker_method_review/md_handoff_DG_v1/system';OUT=ROOT/'linker_method_review/md_handoff_DG_v1/energy_audits_v2';OUT.mkdir(exist_ok=True)
gmx=Path('/home/ayushregmibagale/software/gmx-plumed/bin/gmx');rows=[]
for stage,begin,names in [('nvt',100,['Temperature','Pressure','Potential']),('npt_restrained',250,['Temperature','Pressure','Density','Potential']),('npt',250,['Temperature','Pressure','Density','Potential'])]:
    receipt=WORK/(stage+'.receipt.json')
    if not receipt.exists():continue
    dest=OUT/(stage+'.json')
    if dest.exists():rows.append(json.loads(dest.read_text()));continue
    r=json.loads(receipt.read_text());edr=WORK/(stage+'.edr');assert r['status']=='COMPLETE' and r['outputs'][edr.name]==sha256(edr)
    xvg=OUT/(stage+'.xvg');assert not xvg.exists()
    cmd=[str(gmx),'energy','-f',str(edr),'-o',str(xvg),'-b',str(begin),'-xvg','xmgrace']
    result=subprocess.run(cmd,input='\n'.join(names+['0'])+'\n',text=True,capture_output=True,timeout=120)
    (OUT/(stage+'.log')).write_text(result.stdout+result.stderr);assert result.returncode==0
    time,series=named_energy_series(xvg.read_text(),names);assert len(time)>=3
    assert r['outputs'][edr.name]==sha256(edr)
    report={'stage':stage,'window_ps':[float(time[0]),float(time[-1])],'samples':len(time),'statistics':{name:{'mean':float(series[name].mean()),'sample_std':float(series[name].std(ddof=1))} for i,name in enumerate(names)},'units':{'Temperature':'K','Pressure':'bar','Density':'kg/m^3','Potential':'kJ/mol'},'source_edr_sha256':sha256(edr),'parent_receipt_sha256':sha256(receipt),'xvg_sha256':sha256(xvg),'gmx_sha256':sha256(gmx),'purpose':'Descriptive last-half stage diagnostic; samples are correlated and no equilibrium/convergence claim is made'}
    atomic_json(dest,report);rows.append(report)
atomic_json(ROOT/'linker_method_review/evidence/md_completed_stage_thermodynamics.json',{'status':'DESCRIPTIVE_CONTROL_STATISTICS_NOT_CONVERGENCE','stages':rows,'source_sha256':sha256(Path(__file__))})
print(json.dumps([{'stage':r['stage'],'statistics':r['statistics']} for r in rows],indent=2))
