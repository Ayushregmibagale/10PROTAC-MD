"""Read-only E2 timestamp/replay/imaging audit over the complete core roster.

CA image checks and virtual-point MIC are diagnostics, not replacement physics.
"""
from pathlib import Path
import argparse,ast,collections,concurrent.futures,csv,json,sys,time,shutil
import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations
from MDAnalysis.lib.distances import distance_array
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import ValidationError,sha256,atomic_json,pdb_atoms
from read_only_trajectory import ReadOnlyXTCReader
from e2_geometry import validate_time_alignment,fit_reference_points,WholeProteinGuard
PROJECT=Path('/home/ayushregmibagale/protac-md')


def audit_case(case,output,reference,ensemble,rotation,saved,role_record):
    poi=case['poi'];pose=case['pose_id'];base=PROJECT/'random_screen'/poi/pose
    paths=[base/'prod/prod.tpr',PROJECT/'analysis/random_screen'/poi/pose/'prod_nojump.xtc',PROJECT/'analysis/random_screen'/poi/pose/'metrics.csv']
    hashes={str(p):sha256(p) for p in paths}
    if hashes[str(paths[0])]!=role_record['source_sha256']['tpr']:raise ValidationError('TPR differs from sequence-role audit')
    with paths[2].open() as f:metrics=list(csv.DictReader(line for line in f if not line.startswith('#')))
    mt=np.asarray([float(r['time_ns']) for r in metrics]);dv=np.asarray([float(r['d_protac_vhl_A']) for r in metrics]);dk=np.asarray([float(r['d_protac_kinase_A']) for r in metrics])
    if not np.isfinite(np.column_stack([mt,dv,dk])).all() or np.any(dv<0) or np.any(dk<0):raise ValidationError('Invalid retained metrics')
    u=mda.Universe(str(paths[0]),str(paths[1]),format=ReadOnlyXTCReader)
    whole_guard=WholeProteinGuard(u.select_atoms('protein'));unwrap_calls=0
    anchor=u.select_atoms('moltype Protein_chain_B and name CA');vhl=u.select_atoms('moltype Protein_chain_C and name CA');target=u.select_atoms('moltype Protein_chain_D and name CA');lys=u.select_atoms('moltype Protein_chain_D and resname LYS and name NZ')
    if len(anchor)!=85 or not len(lys):raise ValidationError('ELOC or target lysine selection invalid')
    if not np.array_equal(anchor.indices,u.select_atoms('segid *Protein_chain_B and name CA').indices) or not np.array_equal(lys.indices,u.select_atoms('segid *Protein_chain_D and resname LYS and name NZ').indices):raise ValidationError('Legacy E2 selections differ from identified molecules')
    pairs=[('ELOC_VHL_CA',anchor,vhl),('VHL_POI_CA',vhl,target),('ELOC_POI_CA',anchor,target)]
    rows=[];pdb_c=reference.mean(0);last=None;duplicate_conflicts=0
    for ts in u.trajectory:
        unwrap_calls+=int(whole_guard(ts))
        mobile=anchor.positions;center=mobile.mean(0);R=rotation(reference-pdb_c,mobile-center);mapped=(ensemble-pdb_c)@R.T+center
        transformed,rmsd,_=fit_reference_points(reference,mobile,ensemble)
        nz=lys.positions
        direct=np.sqrt(((nz[:,None]-mapped[None])**2).sum(2)).min(0)
        box=np.asarray(ts.dimensions,dtype=float)
        if box.shape!=(6,) or not np.isfinite(box).all() or np.any(box[:3]<=0):raise ValidationError('Invalid trajectory box')
        periodic=distance_array(mapped,nz,box=box).min(1)
        row={'frame':int(ts.frame),'time_ns':float(ts.time)/1000,'min_E2_direct_A':float(direct.min()),'min_E2_MIC_sensitivity_A':float(periodic.min()),'ELOC_anchor_fit_RMSD_A':rmsd,'rotation_determinant':float(np.linalg.det(R)),'old_vs_checked_mapping_max_A':float(np.max(np.linalg.norm(mapped-transformed,axis=1)))}
        for cutoff in [50,80,100]:
            row[f'reach_{cutoff}A_direct']=float(np.mean(direct<cutoff));row[f'reach_{cutoff}A_MIC_sensitivity']=float(np.mean(periodic<cutoff))
        for name,a,b in pairs:
            d=float(distance_array(a.positions,b.positions).min());p=float(distance_array(a.positions,b.positions,box=box).min());row[name+'_direct_A']=d;row[name+'_periodic_A']=p;row[name+'_gap_A']=d-p
        selected=np.vstack([anchor.positions,vhl.positions,target.positions,nz])
        if last is not None and ts.time==last[0] and not np.array_equal(selected,last[1]):duplicate_conflicts+=1
        last=(float(ts.time),selected.copy());rows.append(row)
    u.trajectory.close()
    tt=np.asarray([r['time_ns'] for r in rows]);alignment_error=None;alignment_failure=None
    try:alignment_error=validate_time_alignment(mt,tt,5.1e-5)
    except ValidationError as exc:alignment_failure=str(exc)
    # Preserve the old positional estimator for investigation; never publish it as aligned when it is not.
    n=min(len(metrics),len(rows));bridge=(dv[:n]<5)&(dk[:n]<5)
    values={'bridge_pct':float(np.mean((dv<5)&(dk<5))*100),'n_kinase_lys':len(lys),'min_e2_dist_mean':float(np.mean([r['min_E2_direct_A'] for r in rows[:n]]))}
    for cutoff in [50,80,100]:
        array=np.asarray([r[f'reach_{cutoff}A_direct'] for r in rows[:n]])
        for label,mask in [('bridge',bridge),('nonbridge',~bridge)]:values[f'reach_{cutoff}A_{label}']=float(array[mask].mean()) if mask.any() else None
    differences={};missing_saved=[]
    for field,value in values.items():
        text=saved.get(field,'')
        original=float(text) if text else None
        if value is None or original is None:
            if value is not None or original is not None:missing_saved.append(field)
        else:differences[field]=float(value-original)
    folder=output/'poses'/poi/pose;folder.mkdir(parents=True)
    with (folder/'frames.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    if any(sha256(p)!=hashes[str(p)] for p in paths):raise ValidationError('Original E2 input changed during audit')
    summary={'poi':poi,'pose_id':pose,'trajectory_frames':len(rows),'original_unwrap_fallback_calls':unwrap_calls,'metric_rows':len(metrics),'time_alignment_error_ns':alignment_error,'time_alignment_failure':alignment_failure,'metric_duplicate_times':int(np.sum(np.diff(mt)==0)),'trajectory_duplicate_times':int(np.sum(np.diff(tt)==0)),'duplicate_coordinate_conflicts_in_audited_selection':duplicate_conflicts,'trajectory_window_ns':[float(tt[0]),float(tt[-1])],'old_positional_estimator_replay':values,'saved_minus_replay_direction':'reported differences are replay minus saved','replay_minus_saved':differences,'missing_value_disagreements':missing_saved,'all_saved_values_reproduced_at_0_001_precision':not missing_saved and all(abs(x)<.001 for x in differences.values()),'maximum_CA_direct_minus_periodic_A':{name:max(r[name+'_gap_A'] for r in rows) for name,_,_ in pairs},'frames_with_CA_image_gap_over_0_001A':{name:sum(r[name+'_gap_A']>.001 for r in rows) for name,_,_ in pairs},'maximum_ELOC_anchor_fit_RMSD_A':max(r['ELOC_anchor_fit_RMSD_A'] for r in rows),'rotation_determinant_range':[min(r['rotation_determinant'] for r in rows),max(r['rotation_determinant'] for r in rows)],'maximum_old_vs_checked_mapping_A':max(r['old_vs_checked_mapping_max_A'] for r in rows),'maximum_virtual_MIC_minimum_change_A':max(abs(r['min_E2_direct_A']-r['min_E2_MIC_sensitivity_A']) for r in rows),'source_input_sha256':hashes,'frames_sha256':sha256(folder/'frames.tsv')}
    atomic_json(folder/'summary.json',summary);return summary


def main(output):
    if output.exists() or output.resolve().is_relative_to(PROJECT):raise ValidationError('Fresh audit output outside original required')
    paths=[Path(__file__),ROOT/'workflow/scripts/lib/e2_geometry.py',ROOT/'workflow/scripts/lib/read_only_trajectory.py',ROOT/'existing_result_audit/contact_core_v1/plan.json',ROOT/'existing_result_audit/protein_roles_core_complete/per_pose_roles.json',ROOT/'existing_result_audit/e2_reference_frame_v1/summary.json',PROJECT/'scripts/02_analysis/16_all_poi_e2_reachability.py',PROJECT/'reference_structures/e2_c85_ensemble.json',PROJECT/'analysis/e2_reachability_all_pois.tsv',PROJECT/'simulations/DDR2/structure/rank1/complex_protein.pdb']
    hashes={str(p):sha256(p) for p in paths};cases=json.loads(paths[3].read_text())['cases'];roles={(r['poi'],r['pose_id']):r for r in json.loads(paths[4].read_text())}
    if len(cases)!=141:raise ValidationError('Expected frozen 141-case core roster')
    static=json.loads(paths[5].read_text())
    if static['max_random_pose_ELOC_CA_RMSD_to_common_reference_A']>1e-6:raise ValidationError('Reference frames differ')
    tree=ast.parse(paths[6].read_text());function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='kabsch_rotation');namespace={'np':np};exec(compile(ast.Module(body=[function],type_ignores=[]),str(paths[6]),'exec'),namespace)
    arrays=[np.asarray(v,dtype=float) for v in json.loads(paths[7].read_text()).values()]
    if any(not np.array_equal(a,arrays[0]) for a in arrays):raise ValidationError('E2 prior arrays differ; do not deduplicate silently')
    reference=np.asarray([a['xyz'] for a in pdb_atoms(paths[9]) if a['record']=='ATOM' and a['chain']=='B' and a['name']=='CA'],dtype=np.float32)
    if reference.shape!=(85,3):raise ValidationError('Incorrect Elongin-C reference')
    with paths[8].open() as f:saved_rows=list(csv.DictReader(f,delimiter='\t'))
    saved={(r['poi'],r['pose_id']):r for r in saved_rows}
    if len(saved)!=len(saved_rows):raise ValidationError('Duplicate saved E2 rows')
    output.mkdir(parents=True);source=output/'sources';source.mkdir()
    for i,p in enumerate(paths):shutil.copyfile(p,source/(str(i)+'_'+p.name))
    atomic_json(output/'intent.json',{'cases':cases,'workers':2,'timestamp_tolerance_ns':5.1e-5,'virtual_point_MIC_policy':'sensitivity only, never substituted as physical correction','source_sha256':hashes,'new_MD_or_training':False})
    start=time.monotonic();outcomes=[];failures=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(audit_case,c,output,reference,arrays[0],namespace['kabsch_rotation'],saved[(c['poi'],c['pose_id'])],roles[(c['poi'],c['pose_id'])]):c for c in cases}
        for f in concurrent.futures.as_completed(futures):
            try:outcomes.append(f.result())
            except Exception as exc:failures.append(dict(futures[f],error=str(exc)))
            if (len(outcomes)+len(failures))%10==0:print(json.dumps({'completed':len(outcomes),'failures':len(failures),'elapsed_s':round(time.monotonic()-start,1)}),flush=True)
    if any(sha256(p)!=hashes[str(p)] for p in paths):raise ValidationError('Source changed during full E2 audit')
    outcomes.sort(key=lambda r:(r['poi'],r['pose_id']));atomic_json(output/'per_pose_summary.json',outcomes)
    report={'status':'COMPLETE_RETAINED_E2_FRAME_AUDIT' if not failures else 'INCOMPLETE_RETAINED_E2_FRAME_AUDIT','cases_requested':len(cases),'cases_completed':len(outcomes),'frames_checked':sum(r['trajectory_frames'] for r in outcomes),'timestamp_aligned_cases':sum(r['time_alignment_failure'] is None for r in outcomes),'saved_E2_reproduced_cases':sum(r['all_saved_values_reproduced_at_0_001_precision'] for r in outcomes),'failures':failures,'cases_with_CA_image_gap':sum(any(v for v in r['frames_with_CA_image_gap_over_0_001A'].values()) for r in outcomes),'maximum_CA_image_gap_A':max((v for r in outcomes for v in r['maximum_CA_direct_minus_periodic_A'].values()),default=None),'maximum_virtual_MIC_minimum_change_A':max((r['maximum_virtual_MIC_minimum_change_A'] for r in outcomes),default=None),'elapsed_s':time.monotonic()-start,'source_sha256':hashes,'summaries_sha256':sha256(output/'per_pose_summary.json'),'original_files_modified':False,'new_MD_or_GNN_training':False,'limitations':['All available frames of the frozen core E2 cohort only; two prediction targets did not feed this eight-target E2 producer.','CA minimum-distance image checks are diagnostics, not an exhaustive heavy-atom interaction or complete coherent assembly proof.','MIC distances for virtual Gaussian points are a sensitivity stress test, not a physical replacement rule.','Reproducing retained E2 values is not validation of the positional prior, experimental labels or GNN biological accuracy.']}
    atomic_json(output/'summary.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['source_sha256','limitations']},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);main(parser.parse_args().output)
