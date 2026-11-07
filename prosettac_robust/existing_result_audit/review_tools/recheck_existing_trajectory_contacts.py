"""Re-analyse one existing trajectory with sequence-identified roles; no new MD."""
from pathlib import Path
import argparse,importlib.util,json,sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import ValidationError,sha256,atomic_json
from read_only_trajectory import ReadOnlyXTCReader
spec=importlib.util.spec_from_file_location('sequence_role_audit',Path(__file__).with_name('audit_existing_protein_roles.py'));roles=importlib.util.module_from_spec(spec);spec.loader.exec_module(roles)


def main(a):
    project=a.project.resolve();out=a.output.resolve()
    if out.exists() or out.is_relative_to(project) or out.is_relative_to(Path('/mnt/d/PROTAC degradation Prediction')):raise ValidationError('Fresh output outside originals required')
    base=project/'random_screen'/a.poi/a.pose_id;analysis=project/'analysis/random_screen'/a.poi/a.pose_id
    tpr=base/'prod/prod.tpr';xtc=analysis/'prod_nojump.xtc';index=analysis/'index.ndx';dv=analysis/'dist_protac_vhl.xvg';dk=analysis/'dist_protac_kinase.xvg'
    paths=[tpr,xtc,index,dv,dk];hashes={str(p):sha256(p) for p in paths}
    cache=ROOT/'workflow/sequence_provenance/uniprot_batches/sjf8240_reference_cohort.json';entries={x['primaryAccession']:x['sequence']['value'] for x in json.loads(cache.read_text())['results']}
    cohort=json.loads((ROOT/'workflow/config/sjf8240_publication_cohort.json').read_text())['protein_groups'];uid=cohort[a.poi] if a.poi in cohort else cohort[a.poi.removesuffix('_full')]
    candidates={k:entries[k] for k in ['P40337','Q15370','Q15369',uid]}
    u=mda.Universe(str(tpr),str(xtc),format=ReadOnlyXTCReader);identified={}
    for name in sorted(set(u.atoms.moltypes)):
        if not str(name).startswith('Protein'):continue
        protein=u.select_atoms('moltype '+str(name));sequence=''.join(seq1(x.resname,undef_code='X') for x in protein.select_atoms('name CA'))
        found=roles.identify(sequence,candidates);identified[found['uniprot']]={'moltype':str(name),'sequence_match':found,'atoms':protein}
    true_vhl=identified['P40337']['atoms'];poi=identified[uid]['atoms'];groups=roles.reader.index_groups(index)
    ligand=u.atoms[np.asarray(groups['Other'])-1];legacy=u.atoms[np.asarray(groups['VHL'])-1];legacy_poi=u.atoms[np.asarray(groups['Kinase'])-1]
    if set(poi.indices)!=set(legacy_poi.indices):raise ValidationError('Legacy POI differs from sequence-identified target; review separately')
    heavy=lambda g:g[g.masses>1.5]
    tv,storedv=roles.reader.distance_series(dv,'Other-VHL',allow_identical_duplicates=True);tk,storedk=roles.reader.distance_series(dk,'Other-Kinase',allow_identical_duplicates=True)
    selected=tv<=a.time_end_ns+1e-8;tv=tv[selected];storedv=storedv[selected];storedk=storedk[selected]
    rows=[]
    for ts in u.trajectory:
        time=float(ts.time)/1000
        if time>a.time_end_ns+1e-8:break
        pairs=[(ligand,legacy),(ligand,true_vhl),(ligand,poi),(heavy(ligand),heavy(true_vhl)),(heavy(ligand),heavy(poi))]
        values=[float(distance_array(left.positions,right.positions,box=ts.dimensions).min()) for left,right in pairs]
        if rows and abs(rows[-1][0]-time)<1e-8:
            if not np.allclose(rows[-1][1:],values,atol=1e-8,rtol=0):raise ValidationError('Conflicting duplicate trajectory frame')
            continue
        rows.append([time,*values])
    data=np.asarray(rows)
    # GROMACS XVG times use seven significant digits; allow text rounding, not a frame shift.
    time_tolerance_ns=max(1e-8,5.1*10**(np.floor(np.log10(max(1.,float(np.max(np.abs(tv))))))-7))
    if len(data)!=len(tv) or not np.allclose(data[:,0],tv,atol=time_tolerance_ns,rtol=0) or not np.isfinite(data).all():
        atomic_json(out.with_name(out.name+'_failure.json'),{'status':'CONTACT_RECHECK_FAILED_CLOSED','trajectory_rows':len(data),'stored_rows':len(tv),'trajectory_times':data[:,0].tolist(),'stored_times':tv.tolist(),'nonfinite_locations':np.argwhere(~np.isfinite(data)).tolist()})
        raise ValidationError('Trajectory and stored distance timestamps differ; never truncate')
    errorv=float(np.abs(data[:,1]-storedv).max());errork=float(np.abs(data[:,3]-storedk).max())
    curves=[]
    for cutoff in [3.5,4.,4.5,5.,5.5,6.]:
        curves.append({'cutoff_A':cutoff,'legacy_named_groups_all_atom_fraction':float(np.mean((data[:,1]<cutoff)&(data[:,3]<cutoff))),'true_VHL_POI_all_atom_fraction':float(np.mean((data[:,2]<cutoff)&(data[:,3]<cutoff))),'true_VHL_POI_heavy_atom_fraction':float(np.mean((data[:,4]<cutoff)&(data[:,5]<cutoff)))})
    if any(sha256(p)!=hashes[str(p)] for p in paths):raise ValidationError('Original input changed during re-analysis')
    out.mkdir(parents=True);np.savetxt(out/'contacts_A.tsv',data,delimiter='\t',header='time_ns\tlegacy_label_VHL_all_A\ttrue_VHL_all_A\tPOI_all_A\ttrue_VHL_heavy_A\tPOI_heavy_A',comments='')
    report={'status':'RETROSPECTIVE_CONTACT_REANALYSIS_NOT_NEW_MD','poi':a.poi,'pose_id':a.pose_id,'frames':len(data),'time_window_ns':[float(data[0,0]),float(data[-1,0])],'time_alignment_tolerance_ns':time_tolerance_ns,'time_max_difference_ns':float(np.abs(data[:,0]-tv).max()),'identified_proteins':{k:{n:v for n,v in x.items() if n!='atoms'} for k,x in identified.items()},'legacy_VHL_group_moltypes':sorted(set(map(str,legacy.moltypes))),'group_atoms':{'ligand_all':len(ligand),'ligand_heavy':len(heavy(ligand)),'true_VHL_all':len(true_vhl),'true_VHL_heavy':len(heavy(true_vhl)),'POI_all':len(poi),'POI_heavy':len(heavy(poi))},'periodic_distance_reproduction_max_error_A':{'legacy_VHL':errorv,'POI':errork},'stored_GROMACS_distances_reproduced_within_0_002_A':errorv<.002 and errork<.002,'cutoff_curves':curves,'source_input_sha256':hashes,'sequence_cache_sha256':sha256(cache),'source_sha256':sha256(Path(__file__)),'original_outputs_modified':False,'new_MD_or_docking_executed':False,'interpretation':'This quantifies atom/role-definition sensitivity for one preselected existing trajectory. It neither replaces saved results nor estimates degradation probability.'}
    atomic_json(out/'summary.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['source_input_sha256','identified_proteins']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['project','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--poi',required=True);p.add_argument('--pose-id',required=True);p.add_argument('--time-end-ns',type=float,default=20.);main(p.parse_args())
