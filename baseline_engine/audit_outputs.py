#!/usr/bin/env python3
"""Report baseline raw output chemistry/geometry without relabeling failed poses.

Writes a separate QC report. Original scores and coordinates are never replaced.
An empty valid subset is a sampling outcome, not an experimental negative.
"""
import argparse,json,sys
import numpy as np
from scipy.spatial import cKDTree
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/lib'))
from protac_qc import parse_scores,pose_features,validate_ligand_identity,sha256,atomic_json,ValidationError,pdb_atoms

def head_engagement(run,pose,selection):
    """Direct bound-core deviations in the actual docked protein frame, no fitting."""
    parts=pose.stem.split('_')
    target=run/'Patchdock_Results'/('docked_'+parts[1]+'_'+parts[2]+'.pdb')
    if not target.exists():return {'status':'BOUND_HEAD_REFERENCE_UNAVAILABLE'}
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    from official_parameters import docked_heads
    paths=[];heads=[];core_masks=[]
    for role,protein in enumerate(['vbc.pdb','protein.pdb']):
        provenance=run/('PT'+str(role)+'_parameter_provenance.json')
        if not provenance.exists():return {'status':'HEAD_PARAMETER_PROVENANCE_UNAVAILABLE'}
        meta=json.loads(provenance.read_text())
        matches=[p for p in run.glob('*.sdf') if sha256(p)==meta['source_sha256']]
        if len(matches)!=1:raise ValidationError('Bound head input cannot be identified uniquely by hash')
        paths.append(matches[0]);head=Chem.SDMolSupplier(str(matches[0]),removeHs=True)[0];heads.append(head)
        protein_xyz=np.array([a['xyz'] for a in pdb_atoms(run/protein) if a['record']=='ATOM' and a['element'] not in {'H','D'}])
        core_masks.append(cKDTree(protein_xyz).query(head.GetConformer().GetPositions())[0]<=4.5)
    reference=docked_heads(str(target),paths,run)
    offset=heads[0].GetNumAtoms();targets=reference.GetConformer().GetPositions()
    atoms=[a for a in pdb_atoms(pose) if a['record']==selection['protac']['record'] and a['chain'] in selection['protac']['chains'] and a['element'] not in {'H','D'}]
    ligand=Chem.MolFromPDBBlock('\n'.join(a['line'] for a in atoms)+'\nEND\n',sanitize=False,removeHs=True)
    ligand=AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles((run/'full.smi').read_text().strip()),ligand)
    Chem.AssignStereochemistryFrom3D(ligand,replaceExistingTags=True)
    mappings=[ligand.GetSubstructMatches(h,useChirality=True) for h in heads]
    results=[]
    for a in mappings[0]:
        for b in mappings[1]:
            if set(a)&set(b):continue
            errors=[np.linalg.norm(ligand.GetConformer().GetPositions()[list(a)]-targets[:offset],axis=1),
                    np.linalg.norm(ligand.GetConformer().GetPositions()[list(b)]-targets[offset:],axis=1)]
            if not all(mask.any() for mask in core_masks):raise ValidationError('Bound contact core is empty')
            rmsd=[float(np.sqrt(np.mean(error[mask]**2))) for error,mask in zip(errors,core_masks)]
            results.append((max(rmsd),rmsd))
    if not results:raise ValidationError('No disjoint stereochemical bound-head mapping')
    best=min(results,key=lambda r:r[0])
    return {'status':'MEASURED_NOT_A_BINDING_LABEL','contact_core_cutoff_A':4.5,
        'comparison':'No superposition; final PROTAC heads versus their docked bound-head coordinates',
        'core_atoms':[int(mask.sum()) for mask in core_masks],
        'vhl_core_rmsd_A':best[1][0],'poi_core_rmsd_A':best[1][1],
        'head_mapping_alternatives':len(results),'docked_head_reference_sha256':sha256(target)}

def audit(run,min_distance=None):
    run=run.resolve();selection=json.loads((run/'selection.json').read_text())
    diagnostic=json.loads((run/'diagnostic.json').read_text()) if (run/'diagnostic.json').exists() else {}
    empty=diagnostic.get('status')=='NO_FEASIBLE_LINKER_POSES_IN_DIAGNOSTIC'
    rows=parse_scores(run/'Patchdock_Results/score.sc',allow_empty=empty)
    if min_distance is not None:selection['qc']={'min_nonbonded_heavy_distance_A':min_distance}
    outcomes=[]
    for row in rows:
        pose=run/'Patchdock_Results'/(row['description']+'.pdb');entry=dict(row,pdb_sha256=sha256(pose))
        try:
            validate_ligand_identity(pose,selection['protac'],run/'full.smi')
            geometry=pose_features(pose,selection);entry.update(geometry)
            entry['head_engagement']=head_engagement(run,pose,selection)
            entry.update(status='GEOMETRY_REJECTED' if geometry['severe_overlap'] else 'MECHANICAL_QC_PASSED',
                         reason=f'Severe protein-ligand overlap (<{geometry["clash_cutoff_A"]} Angstrom)' if geometry['severe_overlap'] else '')
        except (ValidationError,ValueError,OSError) as exc:entry.update(status='CHEMISTRY_OR_OUTPUT_REJECTED',reason=str(exc))
        outcomes.append(entry)
    good=[r for r in outcomes if r['status']=='MECHANICAL_QC_PASSED']
    report={'status':'OUTPUT_QC_PASSED_REQUIRES_SCIENTIFIC_VALIDATION' if good else 'NO_QC_VALID_POSES_IN_RECORDED_DIAGNOSTIC',
        'raw_scored_poses':len(rows),'qc_valid_poses':len(good),'qc_rejected_poses':len(rows)-len(good),
        'raw_pose_geometry_acceptance':len(good)/len(rows) if rows else None,
        'local_pose_count':diagnostic.get('local_pose_count'),
        'local_poses_with_feasible_linkers':diagnostic.get('local_poses_with_feasible_linkers'),
        'conditional_valid_bridge_fraction_5A':sum(r['bridging_5A'] for r in good)/len(good) if good else None,
        'biological_degrader_label':None,'production_release':False,
        'scope':'Chemistry and gross-clash checks only; head engagement, sampling and experimental validation remain required',
        'selection_sha256':sha256(run/'selection.json'),'overlap_cutoff_A':outcomes[0].get('clash_cutoff_A') if outcomes else None,
        'explicit_cutoff_override':min_distance,'score_sha256':sha256(run/'Patchdock_Results/score.sc'),'poses':outcomes}
    filename='baseline_output_qc.json' if min_distance is None else 'baseline_output_qc_min'+str(min_distance).replace('.','p')+'.json'
    atomic_json(run/filename,report);return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--min-distance',type=float)
    a=p.parse_args();print(json.dumps(audit(a.run,a.min_distance),indent=2))
