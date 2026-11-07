"""Time/window/image contracts for retained direct E2 descriptors."""
import numpy as np
from protac_qc import ValidationError
from e2_geometry import validate_time_alignment


def summarize_frames(frames,metric_times,dv,dk,start_ns,end_ns,duplicates='reject',image_policy='reject'):
    if not np.isfinite([start_ns,end_ns]).all() or start_ns<0 or end_ns<=start_ns:
        raise ValidationError('Explicit finite increasing E2 time window required')
    if duplicates not in {'reject','identical'} or image_policy not in {'reject','audit-only'}:
        raise ValidationError('Explicit duplicate/image policies required')
    required=['time_ns','min_E2_direct_A']+[f'reach_{c}A_direct' for c in [50,80,100]]+[n+'_gap_A' for n in ['ELOC_VHL_CA','VHL_POI_CA','ELOC_POI_CA']]
    if frames.dtype.names is None or set(required)-set(frames.dtype.names):raise ValidationError('Incomplete E2 frame packet')
    times=np.atleast_1d(frames['time_ns']);error=validate_time_alignment(metric_times,times,5.1e-5)
    dv=np.asarray(dv,dtype=float);dk=np.asarray(dk,dtype=float)
    if dv.shape!=times.shape or dk.shape!=times.shape or not np.isfinite(np.r_[dv,dk]).all() or np.any(dv<0) or np.any(dk<0):raise ValidationError('Invalid or unmatched contact arrays')
    for name in required:
        if not np.isfinite(frames[name]).all():raise ValidationError('Nonfinite E2 descriptor: '+name)
    if np.any(frames['min_E2_direct_A']<0) or any(np.any((frames[f'reach_{c}A_direct']<0)|(frames[f'reach_{c}A_direct']>1)) for c in [50,80,100]):raise ValidationError('Invalid E2 distance or reach fraction')
    keep=np.ones(len(times),dtype=bool);duplicates_dropped=0
    for i in range(1,len(times)):
        if times[i]!=times[i-1]:continue
        if duplicates=='reject':raise ValidationError('Repeated timestamp requires explicit identical-duplicate policy')
        if dv[i]!=dv[i-1] or dk[i]!=dk[i-1] or any(frames[name][i]!=frames[name][i-1] for name in required):raise ValidationError('Conflicting duplicate frame or contact values')
        keep[i]=False;duplicates_dropped+=1
    window=(times>=start_ns-5.1e-5)&(times<=end_ns+5.1e-5)&keep
    if not window.any() or abs(float(times[window][0])-start_ns)>5.1e-5 or abs(float(times[window][-1])-end_ns)>5.1e-5:raise ValidationError('Trajectory does not cover both requested window endpoints')
    gaps={name:float(np.max(frames[name+'_gap_A'][window])) for name in ['ELOC_VHL_CA','VHL_POI_CA','ELOC_POI_CA']}
    flagged=any(value>.001 for value in gaps.values())
    if flagged and image_policy=='reject':raise ValidationError('Protein image sensitivity detected; coherent assembly review required before E2 release')
    bridge=(dv<5)&(dk<5);values={'mean_min_E2_A':float(frames['min_E2_direct_A'][window].mean())}
    counts={'bridge':int(np.sum(window&bridge)),'nonbridge':int(np.sum(window&~bridge))}
    for cutoff in [50,80,100]:
        for name,state in [('bridge',bridge),('nonbridge',~bridge)]:
            selected=window&state;values[f'reach_{cutoff}A_{name}']=float(frames[f'reach_{cutoff}A_direct'][selected].mean()) if selected.any() else None
        b=values[f'reach_{cutoff}A_bridge'];n=values[f'reach_{cutoff}A_nonbridge'];values[f'delta_reach_{cutoff}A_paired_within_pose']=b-n if b is not None and n is not None else None
    return {'status':'IMAGE_REVIEW_REQUIRED_AUDIT_ONLY' if flagged else 'DESCRIPTIVE_E2_SUMMARY_NOT_TRAINING_RELEASE','selected_frames':int(window.sum()),'window_ns':[start_ns,end_ns],'timestamp_max_difference_ns':error,'duplicate_rows_dropped_full_input':duplicates_dropped,'duplicate_policy':duplicates,'image_policy':image_policy,'maximum_CA_image_gap_A':gaps,'state_counts':counts,'legacy_ELOB_POI_bridge_fraction':float(bridge[window].mean()),'state_definition':'Historical same-frame ligand-to-Elongin-B and ligand-to-POI all-atom distances, strict <5 A; not true VHL contact','values':values,'frame_weighted_descriptive_summary':True,'automatic_training_release':False,'catalytic_probability_claim':False}
