#!/usr/bin/env python3
"""Reparameterize/repack one retained conformer ensemble with an invariant frame."""
import argparse
import json
import os
import shutil
import sys
from pathlib import Path
from rdkit import Chem
from constrained_linker import fixed_head_neighbor,match_heads
from official_parameters import generate
from audit_outputs import head_engagement
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/lib'))
from protac_qc import atomic_json,sha256,validate_ligand_identity,pose_features,ValidationError


def control(source,identifier,dest,seed):
    source=source.resolve();dest=dest.resolve();dest.mkdir(parents=True,exist_ok=False)
    results=dest/'Patchdock_Results';results.mkdir()
    for name in ['full.smi','selection.json','vhl_head.sdf','poi_head.sdf','vbc.pdb','protein.pdb',
                 'PT0_parameter_provenance.json','PT1_parameter_provenance.json']:
        shutil.copy2(source/name,dest/name)
    for name in ['confs_'+identifier+'.sdf','docked_'+identifier+'.sdf','docked_'+identifier+'.pdb']:
        shutil.copy2(source/'Patchdock_Results'/name,results/name)
    for source_head in source.glob('*_head_H.sdf'):
        shutil.copy2(source_head,dest/source_head.name)
    ensemble=results/('confs_'+identifier+'.sdf');molecules=list(Chem.SDMolSupplier(str(ensemble),removeHs=True))
    if len(molecules)<2:raise ValidationError('Frame regression control requires multiple retained conformers')
    heads=[Chem.MolFromMolFile(str(dest/(role+'_head.sdf'))) for role in ['vhl','poi']]
    mapping=match_heads(molecules[0],heads)[0];neighbor=fixed_head_neighbor(molecules[0],mapping[0])
    os.environ['PROTAC_ATTEMPT_ROOT']=str(dest);os.environ['PROTAC_SEED']=str(seed)
    prefix=results/('PT_'+identifier)
    pdb,params=generate(str(ensemble),'PTC',str(prefix),conformers=True,nbr=neighbor)
    original=source/'Patchdock_Results'/('combined_'+identifier+'.pdb')
    protein=[l for l in original.read_text().splitlines() if l.startswith(('ATOM  ','TER'))]
    combined=results/('combined_'+identifier+'.pdb')
    combined.write_text('\n'.join(protein)+'\n'+Path(pdb).read_text()+'\nEND\n')
    os.chdir(results)
    from pyrosetta_adapter import relax
    relax(str(combined),params)
    scored=combined.with_name(combined.stem+'_0001.pdb');selection=json.loads((dest/'selection.json').read_text())
    validate_ligand_identity(scored,selection['protac'],dest/'full.smi')
    engagement=head_engagement(dest,scored,selection)
    if engagement['status']!='MEASURED_NOT_A_BINDING_LABEL' or max(engagement['vhl_core_rmsd_A'],engagement['poi_core_rmsd_A'])>0.05:
        raise ValidationError('Invariant-frame control still moves bound heads')
    report={'status':'FIXED_FRAME_REPACK_CONTROL_PASSED','seed':seed,'input_conformers':len(molecules),
        'original_input_sha256':sha256(original),'ensemble_sha256':sha256(ensemble),
        'neighbor_atom_zero_based':neighbor,'head_engagement':engagement,
        'geometry':pose_features(scored,selection),'scored_pose_sha256':sha256(scored),
        'scope':'Actual retained multi-conformer regression; not an independent binding benchmark'}
    atomic_json(dest/'validation.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--local-id',required=True);p.add_argument('--seed',type=int,required=True)
    a=p.parse_args();control(a.source,a.local_id,a.output,a.seed)
