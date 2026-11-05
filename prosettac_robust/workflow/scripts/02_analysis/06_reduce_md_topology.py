#!/usr/bin/env python3
"""Create a receipt-bound topology for the ordered SOLUTE-only trajectory."""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,sha256,atomic_json


def solute_indices(text,count):
    selected=[];active=False;found=0
    for line in text.splitlines():
        line=line.split(';',1)[0].strip()
        if line.startswith('['):
            match=re.fullmatch(r'\[\s*(.*?)\s*\]',line)
            if not match:raise ValidationError('Malformed index group')
            active=match[1].strip().upper()=='SOLUTE'
            found+=int(active)
        elif active and line:selected.extend(map(int,line.split()))
    if found!=1 or selected!=list(range(1,count+1)):
        raise ValidationError('One contiguous ordered SOLUTE group is required')
    return selected


def reduce(work,stage,gmx,output):
    work=work.resolve();output=output.resolve();parent=work/(stage+'.receipt.json')
    receipt=json.loads(parent.read_text());topology=work/(stage+'.tpr');index=work/'index.ndx'
    preparation=work/'md_preparation.json'
    if receipt['status'] not in {'COMPLETE','SMOKE_COMPLETE'}:raise ValidationError('Stage is incomplete')
    if receipt['outputs'].get(topology.name)!=sha256(topology):raise ValidationError('Parent TPR changed')
    for p in [index,preparation]:
        if receipt['intent']['inputs'].get(p.name)!=sha256(p):raise ValidationError('Executed preparation/index changed')
    count=json.loads(preparation.read_text())['solute_atoms'];solute_indices(index.read_text(),count)
    provenance=output.with_suffix('.provenance.json')
    if output.exists() or provenance.exists():raise ValidationError('Output exists; preserve it and use a fresh name')
    output.parent.mkdir(parents=True,exist_ok=True)
    command=[str(gmx),'convert-tpr','-s',str(topology),'-n',str(index),'-o',str(output)]
    result=subprocess.run(command,input='SOLUTE\n',text=True,capture_output=True)
    log=output.with_suffix('.conversion.log');log.write_text(result.stdout+result.stderr)
    if result.returncode:raise ValidationError('Topology reduction failed')
    import MDAnalysis as mda
    if mda.Universe(str(output)).atoms.n_atoms!=count:raise ValidationError('Reduced atom count differs from complete ordered solute')
    report={'group':'SOLUTE','atom_count':count,'source_sha256':sha256(topology),'target_sha256':sha256(output),
        'index_sha256':sha256(index),'preparation_sha256':sha256(preparation),'parent_receipt_sha256':sha256(parent),
        'gmx_sha256':sha256(gmx),'command':command,'stdin':'SOLUTE','log_sha256':sha256(log),
        'atom_order':'Original indices 1 through solute_atoms, without permutation'}
    atomic_json(provenance,report);print(json.dumps({'status':'ORDERED_SOLUTE_TOPOLOGY_VERIFIED','atom_count':count,'output':str(output)},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['work','gmx','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--stage',required=True)
    a=p.parse_args();reduce(a.work,a.stage,a.gmx,a.output)
