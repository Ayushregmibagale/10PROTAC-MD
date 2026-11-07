#!/usr/bin/env python3
"""Read-only cutoff/selection audit of existing SJF8240 random-screen outputs.

Never executes original scripts, docking, GROMACS or MD. Keep the historical
5 Å same-frame rule primary. Other cutoffs are sensitivity checks without labels.
"""
import argparse,csv,io,json,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,sha256,atomic_json
from gromacs_energy import named_energy_series


def distance_series(path,legend,allow_identical_duplicates=False):
    text=path.read_text()
    if not re.search(r'@\s+xaxis\s+label\s+"Time \(ns\)"',text) or not re.search(r'@\s+yaxis\s+label\s+"Distance \(nm\)"',text):
        raise ValidationError('Explicit ns/nm XVG axes required: '+str(path))
    if allow_identical_duplicates:
        kept=[];previous=None
        for line in text.splitlines():
            stripped=line.strip()
            if not stripped or stripped.startswith(('#','@')):kept.append(line);continue
            values=np.asarray(list(map(float,stripped.split())))
            if previous is not None and values[0]==previous[0]:
                if not np.array_equal(values,previous):raise ValidationError('Conflicting duplicate timestamp; no automatic repair')
                continue
            kept.append(line);previous=values
        text='\n'.join(kept)
    time,columns=named_energy_series(text,[legend])
    values=columns[legend]*10.0
    if np.any(values<0):raise ValidationError('Negative distance')
    return time,values


def joint_statistics(tv,dv,tk,dk,cutoff):
    if not np.isfinite(cutoff) or cutoff<=0:raise ValidationError('Positive finite cutoff required')
    arrays=[np.asarray(x,dtype=float) for x in [tv,dv,tk,dk]]
    if any(x.ndim!=1 or not len(x) or not np.isfinite(x).all() for x in arrays) or len({len(x) for x in arrays})!=1:
        raise ValidationError('Finite equally sized one-dimensional series required; never truncate mismatched data')
    tv,dv,tk,dk=arrays
    if not np.allclose(tv,tk,atol=1e-8,rtol=0) or np.any(np.diff(tv)<=0):raise ValidationError('Distance series timestamps are not aligned and increasing')
    if np.any(dv<0) or np.any(dk<0):raise ValidationError('Negative distance')
    v=dv<cutoff;k=dk<cutoff;joint=v&k
    pv=float(v.mean());pk=float(k.mean());pj=float(joint.mean())
    return {'frames':len(tv),'cutoff_A':float(cutoff),'vhl_fraction':pv,'poi_fraction':pk,'same_frame_bridge_fraction':pj,'product_of_marginals':pv*pk,'joint_minus_product':pj-pv*pk,'exact_threshold_vhl_frames':int(np.sum(dv==cutoff)),'exact_threshold_poi_frames':int(np.sum(dk==cutoff))}


def index_groups(path):
    groups={};name=None
    for line in path.read_text().splitlines():
        line=line.split(';',1)[0].strip()
        if line.startswith('['):
            name=line.strip('[] ')
            if name in groups:raise ValidationError('Duplicate index group')
            groups[name]=[]
        elif line and name:groups[name].extend(map(int,line.split()))
    chosen={name:groups.get(name,[]) for name in ['Other','VHL','Kinase']}
    if any(not ids or len(set(ids))!=len(ids) or min(ids)<1 for ids in chosen.values()):raise ValidationError('Invalid or missing contact index group')
    names=list(chosen)
    if any(set(chosen[a])&set(chosen[b]) for i,a in enumerate(names) for b in names[i+1:]):raise ValidationError('Overlapping ligand/protein contact groups')
    return chosen


def selected_gro_atoms(path,groups):
    wanted=set().union(*map(set,groups.values()));selected={};last=max(wanted)
    with path.open() as handle:
        handle.readline();count=int(handle.readline())
        if last>count:raise ValidationError('Contact index exceeds GRO atom count')
        for index in range(1,last+1):
            line=handle.readline()
            if len(line)<44:raise ValidationError('Truncated GRO input')
            if index in wanted:selected[index]=(line[5:10].strip(),line[10:15].strip())
    out={}
    for name,ids in groups.items():
        atoms=[selected[i] for i in ids]
        hydrogens=sum(re.sub(r'^\d+','',atom).startswith('H') for _,atom in atoms)
        out[name]={'atoms':len(ids),'hydrogen_name_count':hydrogens,'hydrogen_detection':'GRO atom-name heuristic; exact elemental/mass identity requires topology review','residue_names':sorted({res for res,_ in atoms})}
    return out


def audit(project,output,cutoffs,pois=None,allow_identical_duplicates=False,time_end_ns=None):
    project=project.resolve();output=output.resolve()
    if output.exists():raise ValidationError('Output exists; use a fresh audit directory')
    if output.is_relative_to(project) or output.is_relative_to(Path('/mnt/d/PROTAC degradation Prediction')):raise ValidationError('Audit output must be outside original projects')
    if 5.0 not in cutoffs or len(set(cutoffs))!=len(cutoffs):raise ValidationError('Explicit unique cutoffs must include the unchanged historical primary 5 Å rule')
    manifest=project/'random_screen/manifest.tsv';manifest_hash=sha256(manifest)
    with manifest.open() as handle:declared=[r for r in csv.DictReader(handle,delimiter='\t') if r['stage']=='ANALYSIS_DONE']
    all_declared_count=len(declared)
    if pois is not None:
        absent=set(pois)-{r['poi'] for r in declared}
        if absent:raise ValidationError('Requested targets absent from completed manifest: '+str(sorted(absent)))
        declared=[r for r in declared if r['poi'] in pois]
    if not declared or len({(r['poi'],r['pose_id']) for r in declared})!=len(declared):raise ValidationError('Nonempty unique completed-pose manifest required')
    rows=[];poses=[];failures=[]
    for row in declared:
        poi=row['poi'];pose=row['pose_id']
        if any(Path(value).name!=value or value in {'.','..'} for value in [poi,pose]):raise ValidationError('Unsafe manifest identifier')
        folder=project/'analysis/random_screen'/poi/pose;gro=project/'random_screen'/poi/pose/'topol/complex_ions.gro'
        paths=[folder/'dist_protac_vhl.xvg',folder/'dist_protac_kinase.xvg',folder/'metrics.csv',folder/'index.ndx',gro]
        try:
            hashes={str(p):sha256(p) for p in paths}
            rawv=np.loadtxt(paths[0],comments=('#','@'),ndmin=2);rawk=np.loadtxt(paths[1],comments=('#','@'),ndmin=2)
            if rawv.shape!=rawk.shape or rawv.shape[1]!=2 or not np.isfinite(rawv).all() or not np.isfinite(rawk).all() or not np.allclose(rawv[:,0],rawk[:,0],atol=1e-8,rtol=0):raise ValidationError('Raw distance series do not reconcile; never truncate')
            raw_primary=float(np.mean((rawv[:,1]*10<5.)&(rawk[:,1]*10<5.)))
            text=paths[2].read_text();csvrows=list(csv.DictReader(line for line in text.splitlines() if not line.startswith('#')))
            values=np.array([[float(r[k]) for k in ['time_ns','d_protac_vhl_A','d_protac_kinase_A']] for r in csvrows])
            if values.shape!=(len(rawv),3) or not np.allclose(values[:,0],rawv[:,0],atol=5.01e-5,rtol=0) or not np.allclose(values[:,1:],np.column_stack([rawv[:,1],rawk[:,1]])*10,atol=5.01e-5,rtol=0):raise ValidationError('Saved CSV does not reconcile to paired XVG data and units')
            match=re.search(r'# Bridging .*?:\s*([\d.]+)%',text)
            if not match:raise ValidationError('Historical bridging header missing')
            header_pct=float(match[1]);delta=abs(header_pct-100*raw_primary)
            if delta>0.050001:raise ValidationError('Historical reported bridging percentage does not reconcile')
            tv,dv=distance_series(paths[0],'Other-VHL',allow_identical_duplicates);tk,dk=distance_series(paths[1],'Other-Kinase',allow_identical_duplicates)
            canonical_frames=len(tv)
            joint_statistics(tv,dv,tk,dk,5.)
            if time_end_ns is not None:
                keep=tv<=time_end_ns+1e-8;tv=tv[keep];dv=dv[keep];tk=tk[keep];dk=dk[keep]
            primary=joint_statistics(tv,dv,tk,dk,5.)
            groups=index_groups(paths[3]);selection=selected_gro_atoms(gro,groups)
            step=np.diff(tv)
            poses.append({'poi':poi,'pose_id':pose,'source_pdb':row['source_pdb'],'raw_frames':len(rawv),'canonical_full_frames':canonical_frames,'removed_identical_duplicate_frames':len(rawv)-canonical_frames,'raw_time_last_ns':float(rawv[-1,0]),'frames':len(tv),'time_first_ns':float(tv[0]),'time_last_ns':float(tv[-1]),'median_sample_interval_ns':float(np.median(step)) if len(step) else None,'uniform_sample_spacing':bool(np.allclose(step,step[0],atol=1e-8,rtol=0)) if len(step) else True,'legacy_header_bridging_pct':header_pct,'recomputed_bridging_pct':100*primary['same_frame_bridge_fraction'],'index_group_composition':selection,'input_sha256':hashes})
            for cutoff in cutoffs:rows.append(dict(poi=poi,pose_id=pose,**joint_statistics(tv,dv,tk,dk,cutoff)))
            if any(sha256(p)!=hashes[str(p)] for p in paths):raise ValidationError('Original input changed during audit')
        except (OSError,ValueError,KeyError) as exc:
            failures.append({'poi':poi,'pose_id':pose,'reason':str(exc)})
    if sha256(manifest)!=manifest_hash:raise ValidationError('Original manifest changed during audit')
    aggregate=[]
    for poi in sorted({r['poi'] for r in rows}):
        for cutoff in cutoffs:
            group=[r for r in rows if r['poi']==poi and r['cutoff_A']==cutoff]
            if not group:continue
            weights=np.asarray([r['frames'] for r in group])
            aggregate.append({'poi':poi,'cutoff_A':cutoff,'poses':len(group),'total_frames':int(weights.sum()),'equal_pose_mean_bridge_fraction':float(np.mean([r['same_frame_bridge_fraction'] for r in group])),'pooled_frame_bridge_fraction':float(np.average([r['same_frame_bridge_fraction'] for r in group],weights=weights))})
    output.mkdir(parents=True)
    for name,table in [('per_pose_cutoff.tsv',rows),('per_target_cutoff.tsv',aggregate)]:
        if table:
            with (output/name).open('w') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(table[0]),delimiter='\t');writer.writeheader();writer.writerows(table)
    atomic_json(output/'pose_provenance.json',poses)
    report={'status':'COMPLETE_READ_ONLY_EXISTING_OUTPUT_AUDIT' if not failures else 'INCOMPLETE_AUDIT_REVIEW_REQUIRED','original_project':str(project),'manifest_sha256':manifest_hash,'all_manifest_completed_poses':all_declared_count,'selected_targets':pois,'identical_duplicate_policy':'drop exact repeats in derivative only' if allow_identical_duplicates else 'reject','time_end_ns':time_end_ns,'declared_completed_poses':len(declared),'audited_poses':len(poses),'failures':failures,'primary_rule':'(min ligand-to-VHL distance <5 Å) AND (min ligand-to-POI distance <5 Å) at aligned timestamps','cutoff_sensitivity_A':cutoffs,'original_outputs_modified':False,'new_PRosettaC_or_MD_executed':False,'labels_used_for_cutoff_selection':False,'contact_scope':'Actual stored Other/VHL/Kinase index groups; composition recorded, not silently redefined as heavy atoms','limitations':['Reuses stored distances; does not independently recompute periodic distances or prove correct index-to-TPR correspondence.','GRO hydrogen counts use atom names; topology is needed for exact element/mass identity.','Pose-sampled frame occupancy is not equilibrium binding probability, degradation probability or a free energy.','Correlated frames and bridge-selected starting poses preclude treating every frame/pose as an independent biological replicate.'],'source_sha256':sha256(Path(__file__)),'outputs_sha256':{p.name:sha256(p) for p in output.iterdir() if p.is_file()}}
    atomic_json(output/'summary.json',report);print(json.dumps({k:v for k,v in report.items() if k in ['status','declared_completed_poses','audited_poses','failures']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--project',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--cutoffs',type=float,nargs='+',default=[3.5,4.,4.5,5.,5.5,6.]);parser.add_argument('--poi',nargs='+');parser.add_argument('--allow-identical-duplicates',action='store_true');parser.add_argument('--time-end-ns',type=float);args=parser.parse_args();audit(args.project,args.output,args.cutoffs,args.poi,args.allow_identical_duplicates,args.time_end_ns)
