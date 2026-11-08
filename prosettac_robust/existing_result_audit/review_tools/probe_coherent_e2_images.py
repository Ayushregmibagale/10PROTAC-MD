"""Bounded retained-data comparison of reference/contact imaging conventions.

Includes all 12 previously flagged core cases plus a prespecified MET control.
No virtual-point wrapping, new simulation, original edit or GNN training.
"""
from pathlib import Path
import ast,collections,concurrent.futures,csv,json,os,shutil,sys,time
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import ValidationError,sha256,atomic_json,pdb_atoms
from read_only_trajectory import ReadOnlyXTCReader
from e2_geometry import WholeProteinTree,fit_reference_points,validate_time_alignment
from coherent_images import assemble,cell_vectors
PROJECT=Path('/home/ayushregmibagale/protac-md')
ALIASES={'HID':'H','HIE':'H','HIP':'H','HSD':'H','HSE':'H','HSP':'H','CYX':'C','MSE':'M'}


def seq(resnames):return ''.join(seq1(str(name),custom_map=ALIASES,undef_code='X') for name in resnames)

def descriptor(points,lysines):
    distances=np.sqrt(np.sum((points[:,None]-lysines[None])**2,axis=2)).min(1)
    return distances,float(distances.min()),[float(np.mean(distances<c)) for c in [50,80,100]]


def work(case,ensemble,rotation,output):
    poi=case['poi'];pose=case['pose_id'];base=PROJECT/'random_screen'/poi/pose;analysis=PROJECT/'analysis/random_screen'/poi/pose
    files=[base/'prod/prod.tpr',analysis/'prod_nojump.xtc',base/'structure/complex_protein.pdb',analysis/'metrics.csv',ROOT/'existing_result_audit/e2_full_frames_v3/poses'/poi/pose/'frames.tsv',ROOT/'existing_result_audit/contact_core_v1/poses'/poi/pose/'contacts_A.tsv']
    hashes={str(p):sha256(p) for p in files};atoms=pdb_atoms(files[2]);refs={chain:[a for a in atoms if a['record']=='ATOM' and a['chain']==chain and a['name']=='CA'] for chain in ['B','C','D']}
    reference={chain:np.asarray([a['xyz'] for a in rows],dtype=np.float32) for chain,rows in refs.items()}
    u=mda.Universe(str(files[0]),str(files[1]),format=ReadOnlyXTCReader);whole=WholeProteinTree(u.select_atoms('protein'));groups={chain:u.select_atoms('moltype Protein_chain_'+chain+' and name CA') for chain in ['B','C','D']};lys=u.select_atoms('moltype Protein_chain_D and resname LYS and name NZ')
    for chain in groups:
        if len(groups[chain])!=len(refs[chain]) or seq(groups[chain].resnames)!=seq([a['resname'] for a in refs[chain]]):raise ValidationError('Reference/TPR CA correspondence differs: '+chain)
    with files[3].open() as f:metrics=list(csv.DictReader(line for line in f if not line.startswith('#')))
    previous=np.genfromtxt(files[4],delimiter='\t',names=True);contacts=np.genfromtxt(files[5],delimiter='\t',names=True)
    focus=case['focus_frame'];rows=[]
    for ts in u.trajectory:
        if ts.time/1000>20+5.1e-5:break
        whole(ts);a=groups['B'].positions.copy();v=groups['C'].positions.copy();p=groups['D'].positions.copy();nz=lys.positions.copy();box=ts.dimensions.copy()
        im=assemble(a,v,p,reference['B'],reference['C'],reference['D'],box)
        ref_center=reference['B'].mean(0);center=a.mean(0);R=rotation(reference['B']-ref_center,a-center);points=(ensemble-ref_center)@R.T+center
        variants={'old':nz,'contact_image':nz+im['poi_shift_A'],'reference_image':nz+im['reference_poi_shift_A']};row={'frame':int(ts.frame),'time_ns':float(ts.time)/1000,'image_methods_agree':int(im['POI_image_methods_agree']),'ELOC_fit_RMSD_A':im['ELOC_fit_RMSD_A'],'VHL_reference_fit_RMSD_A':im['VHL_reference_fit_RMSD_A']}
        for name,position in variants.items():
            _,minimum,reach=descriptor(points,position);row[name+'_min_E2_A']=minimum
            for cutoff,value in zip([50,80,100],reach):row[f'{name}_reach_{cutoff}A']=value
        periodic_errors=[]
        for name,left,right,lifted_left,lifted_right in [('ELOC_VHL',a,v,a,im['vhl']),('VHL_POI',v,p,im['vhl'],im['poi']),('ELOC_POI',a,p,a,im['poi'])]:
            before=float(distance_array(left,right,box=box).min());after=float(distance_array(lifted_left,lifted_right,box=box).min());direct=float(distance_array(lifted_left,lifted_right).min());row[name+'_candidate_direct_minus_periodic_A']=direct-after;periodic_errors.append(abs(before-after))
        row['periodic_distance_preservation_error_A']=max(periodic_errors)
        for name,values in [('VHL',im['vhl_lattice_shift']),('POI_contact',im['poi_lattice_shift']),('POI_reference',im['reference_poi_lattice_shift'])]:
            for axis,value in zip('abc',values):row[name+'_lattice_'+axis]=int(value)
        tested=ts.frame in {0,50,100,150,200,focus};row['box_translation_tested']=int(tested);row['translation_invariance_error_A']='';row['reach_translation_changed']=''
        if tested:
            vectors=cell_vectors(box);sa=np.array([1,-1,1])@vectors;sv=np.array([-1,1,0])@vectors;sp=np.array([0,-1,1])@vectors
            shifted=assemble(a+sa,v+sv,p+sp,reference['B'],reference['C'],reference['D'],box)
            canonical,_,_=fit_reference_points(reference['B'],a,ensemble);translated,_,_=fit_reference_points(reference['B'],a+sa,ensemble)
            d0,_,r0=descriptor(canonical,nz+im['poi_shift_A']);d1,_,r1=descriptor(translated,nz+sp+shifted['poi_shift_A'])
            row['translation_invariance_error_A']=max(float(np.max(np.abs(d0-d1))),float(np.max(np.linalg.norm(shifted['vhl']-sa-im['vhl'],axis=1))),float(np.max(np.linalg.norm(shifted['poi']-sa-im['poi'],axis=1))))
            row['reach_translation_changed']=int(r0!=r1)
        rows.append(row)
    u.trajectory.close();times=[r['time_ns'] for r in rows];mt=[float(r['time_ns']) for r in metrics];validate_time_alignment(mt,times,5.1e-5);validate_time_alignment(contacts['time_ns'],times,5.1e-5)
    if len(rows)!=201:raise ValidationError('Expected complete 201-frame core window')
    old_error=float(np.max(np.abs(previous['min_E2_direct_A']-np.asarray([r['old_min_E2_A'] for r in rows]))))
    if old_error>.001:raise ValidationError('Original descriptor replay changed')
    if any(r['periodic_distance_preservation_error_A']>.001 for r in rows):raise ValidationError('Candidate changed a periodic protein distance')
    tested_rows=[r for r in rows if r['box_translation_tested']]
    if any(r['translation_invariance_error_A']>.001 for r in tested_rows):raise ValidationError('Candidate is not invariant to box translations')
    states={'legacy_ELOB_POI':np.array([float(r['d_protac_vhl_A'])<5 and float(r['d_protac_kinase_A'])<5 for r in metrics]),'true_VHL_POI_heavy':(contacts['true_VHL_heavy_A']<5)&(contacts['POI_heavy_A']<5)}
    comparisons={}
    for variant in variants:
        comparisons[variant]={'mean_min_E2_A':float(np.mean([r[variant+'_min_E2_A'] for r in rows])),'conditional':{}}
        for name,mask in states.items():
            sub={'bridge_frames':int(mask.sum()),'nonbridge_frames':int((~mask).sum())}
            for cutoff in [50,80,100]:
                values=np.array([r[f'{variant}_reach_{cutoff}A'] for r in rows]);b=float(values[mask].mean()) if mask.any() else None;n=float(values[~mask].mean()) if (~mask).any() else None;sub[f'reach_{cutoff}A_bridge']=b;sub[f'reach_{cutoff}A_nonbridge']=n;sub[f'delta_{cutoff}A_paired']=b-n if b is not None and n is not None else None
            comparisons[variant]['conditional'][name]=sub
    folder=output/'poses'/poi/pose;folder.mkdir(parents=True)
    with (folder/'frames.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    if any(sha256(p)!=hashes[str(p)] for p in files):raise ValidationError('Original source changed')
    report={'status':'CANDIDATE_IMAGING_SENSITIVITY_NOT_A_TRAINING_RELEASE','poi':poi,'pose_id':pose,'role':case['role'],'frames':len(rows),'POI_image_conventions_agree_frames':sum(r['image_methods_agree'] for r in rows),'max_VHL_reference_fit_RMSD_A':max(r['VHL_reference_fit_RMSD_A'] for r in rows),'max_candidate_CA_direct_minus_periodic_A':max(r[n+'_candidate_direct_minus_periodic_A'] for r in rows for n in ['ELOC_VHL','VHL_POI','ELOC_POI']),'max_periodic_distance_preservation_error_A':max(r['periodic_distance_preservation_error_A'] for r in rows),'box_translation_test_frames':len(tested_rows),'max_translation_invariance_error_A':max(r['translation_invariance_error_A'] for r in tested_rows),'reachability_translation_change_frames':sum(r['reach_translation_changed'] for r in tested_rows),'old_descriptor_replay_max_error_A':old_error,'comparisons':comparisons,'source_input_sha256':hashes,'frames_sha256':sha256(folder/'frames.tsv'),'automatic_training_release':False,'new_MD_or_GNN_training':False};atomic_json(folder/'summary.json',report);return report


def main(output):
    if output.exists() or output.resolve().is_relative_to(PROJECT):raise ValidationError('Fresh isolated candidate output required')
    paths=[Path(__file__),ROOT/'workflow/scripts/lib/coherent_images.py',ROOT/'workflow/scripts/lib/e2_geometry.py',ROOT/'workflow/scripts/lib/read_only_trajectory.py',ROOT/'existing_result_audit/e2_full_frames_v3/per_pose_summary.json',PROJECT/'scripts/02_analysis/16_all_poi_e2_reachability.py',PROJECT/'reference_structures/e2_c85_ensemble.json']
    hashes={str(p):sha256(p) for p in paths};records=json.loads(paths[4].read_text());cases=[]
    for r in records:
        if not any(r['frames_with_CA_image_gap_over_0_001A'].values()):continue
        frames=np.genfromtxt(ROOT/'existing_result_audit/e2_full_frames_v3/poses'/r['poi']/r['pose_id']/'frames.tsv',delimiter='\t',names=True);gap=np.maximum.reduce([frames[k+'_gap_A'] for k in ['ELOC_VHL_CA','VHL_POI_CA','ELOC_POI_CA']]);cases.append({'poi':r['poi'],'pose_id':r['pose_id'],'role':'all_previously_flagged_cases','focus_frame':int(frames['frame'][np.argmax(gap)])})
    if len(cases)!=12:raise ValidationError('Frozen flagged set differs from 12 cases')
    cases.append({'poi':'MET','pose_id':'pose002','role':'prespecified_clear_control_not_independent_benchmark','focus_frame':0})
    tree=ast.parse(paths[5].read_text());function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='kabsch_rotation');namespace={'np':np};exec(compile(ast.Module(body=[function],type_ignores=[]),str(paths[5]),'exec'),namespace)
    arrays=[np.asarray(v,dtype=float) for v in json.loads(paths[6].read_text()).values()]
    if any(not np.array_equal(a,arrays[0]) for a in arrays):raise ValidationError('E2 prior multiplicity changed')
    output.mkdir(parents=True);sources=output/'sources';sources.mkdir()
    for i,p in enumerate(paths):shutil.copyfile(p,sources/(str(i)+'_'+p.name))
    atomic_json(output/'intent.json',{'status':'GEOMETRIC_CONVENTION_SENSITIVITY_ONLY','cases':cases,'PID':os.getpid(),'workers':2,'input_sha256':hashes,'rules':['ELOC fit anchors E3 reference image','VHL shifted to nearest reference-predicted centroid image','POI shifted by closest periodic VHL/POI CA pair','Independent initial-reference-pose POI image is compared, not silently adopted','Gaussian virtual points are not wrapped','No experimental labels used to select images']})
    start=time.monotonic();rows=[];failures=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(work,c,arrays[0],namespace['kabsch_rotation'],output):c for c in cases}
        for f in concurrent.futures.as_completed(futures):
            try:rows.append(f.result())
            except Exception as exc:failures.append(dict(futures[f],error=str(exc)))
            print(json.dumps({'completed':len(rows),'failures':len(failures),'elapsed_s':round(time.monotonic()-start,1)}),flush=True)
    if any(sha256(p)!=hashes[str(p)] for p in paths):raise ValidationError('Candidate code/input changed')
    rows.sort(key=lambda r:(r['poi'],r['pose_id']));atomic_json(output/'per_case.json',rows)
    summary={'status':'COMPLETE_CANDIDATE_SENSITIVITY_NOT_SCIENTIFIC_RELEASE' if not failures else 'INCOMPLETE_CANDIDATE_SENSITIVITY','cases_requested':13,'cases_completed':len(rows),'frames':sum(r['frames'] for r in rows),'failures':failures,'cases_with_complete_POI_image_method_agreement':sum(r['POI_image_conventions_agree_frames']==r['frames'] for r in rows),'frames_with_POI_image_method_disagreement':sum(r['frames']-r['POI_image_conventions_agree_frames'] for r in rows),'maximum_periodic_distance_preservation_error_A':max((r['max_periodic_distance_preservation_error_A'] for r in rows),default=None),'maximum_translation_invariance_error_A':max((r['max_translation_invariance_error_A'] for r in rows),default=None),'box_translation_frames_tested':sum(r['box_translation_test_frames'] for r in rows),'elapsed_s':time.monotonic()-start,'original_files_modified':False,'new_MD_or_GNN_training':False,'automatic_training_release':False,'source_sha256':hashes,'records_sha256':sha256(output/'per_case.json'),'limitations':['Retrospective all-flagged-case diagnostic plus one previously inspected control, not independent validation.','Reference and nearest-contact conventions are explicit geometric choices; agreement does not establish a productive E2 or catalytic probability.','A dissociated or image-ambiguous assembly requires separate interpretation; no model scores are changed here.']};atomic_json(output/'summary.json',summary);print(json.dumps({k:v for k,v in summary.items() if k not in ['source_sha256','limitations']},indent=2))

if __name__=='__main__':
    p=__import__('argparse').ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
