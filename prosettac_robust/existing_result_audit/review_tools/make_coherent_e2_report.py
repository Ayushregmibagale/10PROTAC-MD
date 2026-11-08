"""Export descriptive E2 imaging tables and vector figures; never fit a model."""
from pathlib import Path
import collections,csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,atomic_json,ValidationError
source=ROOT/'existing_result_audit/ligand_connected_e2_cohort_v1';output=ROOT/'existing_result_audit/coherent_e2_publication_v1'
if output.exists():raise ValidationError('Fresh publication artifact output required')
paths=[source/'summary.json',source/'per_case.json',ROOT/'existing_result_audit/retained_stage_evidence_v1/per_stage.json'];hashes={str(p):sha256(p) for p in paths};summary=json.loads(paths[0].read_text());cases=json.loads(paths[1].read_text());stages=json.loads(paths[2].read_text())
if summary['cases_completed']!=141 or summary['failures'] or any(c['frames']!=201 for c in cases):raise ValidationError('Incomplete coherent analysis')
exceptions={(r['poi'],r['pose_id']) for r in stages if r['stage']=='em' and r['log'].get('available') and not r['log']['positive_EM_convergence_statement']}
poses=[];disagreements=[]
for case in cases:
 poi=case['poi'];pose=case['pose_id'];path=source/'poses'/poi/pose/'frames.tsv';contacts=ROOT/'existing_result_audit/contact_core_v1/poses'/poi/pose/'contacts_A.tsv';hashes.update({str(p):sha256(p) for p in [path,contacts]});frame=np.genfromtxt(path,delimiter='\t',names=True);contact=np.genfromtxt(contacts,delimiter='\t',names=True)
 if not np.allclose(frame['time_ns'],contact['time_ns'],atol=5.1e-5,rtol=0):raise ValidationError('Unaligned figure/bridge frames')
 mismatch=frame['ligand_contact_image_agrees']==0;bridge=(contact['true_VHL_heavy_A']<5)&(contact['POI_heavy_A']<5)
 if mismatch.any():disagreements.append({'poi':poi,'pose_id':pose,'frames':int(mismatch.sum()),'true_VHL_heavy_bridge_frames_among_disagreements':int((mismatch&bridge).sum())})
 row={'poi':poi,'pose_id':pose,'EM_force_criterion_exception':int((poi,pose) in exceptions),'frames':case['frames'],'ligand_contact_agreement_frames':case['ligand_contact_image_agree_frames'],'ligand_reference_agreement_frames':case['ligand_reference_image_agree_frames']}
 for variant in ['old','contact_image','reference_image','ligand_image']:row[variant+'_mean_min_E2_A']=case['comparisons'][variant]['mean_min_E2_A']
 row['ligand_minus_old_mean_min_E2_A']=row['ligand_image_mean_min_E2_A']-row['old_mean_min_E2_A'];poses.append(row)
aggregates=[];conditional=[];sensitivity=[]
for poi in sorted({r['poi'] for r in poses}):
 selected=[r for r in poses if r['poi']==poi];row={'poi':poi,'poses':len(selected)}
 for variant in ['old','contact_image','reference_image','ligand_image']:row[variant+'_mean_min_E2_A']=float(np.mean([r[variant+'_mean_min_E2_A'] for r in selected]))
 row['ligand_minus_old_mean_min_E2_A']=row['ligand_image_mean_min_E2_A']-row['old_mean_min_E2_A'];aggregates.append(row)
 clean=[r for r in selected if not r['EM_force_criterion_exception']]
 if len(clean)!=len(selected):sensitivity.append({'poi':poi,'all_poses':len(selected),'without_EM_exception_poses':len(clean),'old_all_mean_A':row['old_mean_min_E2_A'],'old_without_exception_mean_A':float(np.mean([r['old_mean_min_E2_A'] for r in clean])),'ligand_all_mean_A':row['ligand_image_mean_min_E2_A'],'ligand_without_exception_mean_A':float(np.mean([r['ligand_image_mean_min_E2_A'] for r in clean])),'main_analysis_drops_exception':False})
 for variant in ['old','ligand_image']:
  for state in ['legacy_ELOB_POI','true_VHL_POI_heavy']:
   for cutoff in [50,80,100]:
    values=[r['comparisons'][variant]['conditional'][state][f'delta_{cutoff}A_paired'] for r in cases if r['poi']==poi];available=[x for x in values if x is not None]
    conditional.append({'poi':poi,'image_convention':variant,'state_mask':state,'threshold_A':cutoff,'total_poses':len(values),'paired_state_poses':len(available),'mean_within_pose_bridge_minus_nonbridge':float(np.mean(available)) if available else None,'automatic_GNN_training_release':False})
output.mkdir(parents=True)
def table(name,rows):
 with (output/name).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
for name,rows in [('per_pose_imaging.tsv',poses),('per_target_imaging.tsv',aggregates),('paired_state_sensitivity.tsv',conditional),('EM_exception_sensitivity.tsv',sensitivity),('image_disagreement_states.tsv',disagreements)]:table(name,rows)
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig,(left,right)=plt.subplots(1,2,figsize=(11,4.8),gridspec_kw={'width_ratios':[1.15,1]});palette=plt.get_cmap('tab10')
for i,poi in enumerate(sorted({r['poi'] for r in poses})):
 selected=[r for r in poses if r['poi']==poi];left.scatter([r['old_mean_min_E2_A'] for r in selected],[r['ligand_image_mean_min_E2_A'] for r in selected],s=20,color=palette(i),alpha=.8,label=poi.replace('_full',''),edgecolors='none')
maximum=max(max(r['old_mean_min_E2_A'],r['ligand_image_mean_min_E2_A']) for r in poses)+4
left.plot([0,maximum],[0,maximum],color='.5',lw=1,ls='--');left.set(xlim=(0,maximum),ylim=(0,maximum),xlabel='Original fragment images: mean minimum distance (Å)',ylabel='Ligand-connected images: mean minimum distance (Å)',title='A  Individual retained poses (n=141)');left.legend(frameon=False,ncol=2,fontsize=8,loc='lower right')
y=np.arange(len(aggregates));old=[r['old_mean_min_E2_A'] for r in aggregates];new=[r['ligand_image_mean_min_E2_A'] for r in aggregates]
for i,(a,b) in enumerate(zip(old,new)):right.plot([a,b],[i,i],color='.65',lw=1.5)
right.scatter(old,y,color='.5',s=28,label='Original');right.scatter(new,y,color='#1767a6',s=28,label='Ligand-connected (candidate)');right.set_yticks(y,[r['poi'].replace('_full','') for r in aggregates]);right.invert_yaxis();right.set(xlabel='Equal-pose mean minimum modeled E2–lysine distance (Å)',title='B  Target summaries');right.legend(frameon=False,fontsize=8,loc='lower left')
for ax in [left,right]:ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.12)
fig.suptitle('Periodic image convention changes E2 geometric descriptors',fontsize=12);fig.text(.5,.015,'Descriptive reanalysis; no independent-replica confidence intervals or new GNN predictions.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.045,1,.94])
for ext in ['pdf','svg','png']:fig.savefig(output/('e2_imaging_sensitivity.'+ext),dpi=300,bbox_inches='tight')
plt.close(fig)
if any(sha256(Path(p))!=h for p,h in hashes.items()):raise ValidationError('Figure source changed')
report={'status':'DESCRIPTIVE_PUBLICATION_ARTIFACTS_NOT_NEW_MODEL_PERFORMANCE','cases':141,'targets':len(aggregates),'frames_each':201,'poses_with_mean_E2_change_over_0_001_A':sum(abs(r['ligand_minus_old_mean_min_E2_A'])>.001 for r in poses),'image_disagreements':disagreements,'independent_replica_CI_claimed':False,'old_rows_or_predictions_changed':False,'new_MD_or_GNN_training':False,'automatic_training_release':False,'source_sha256':hashes,'script_sha256':sha256(Path(__file__)),'outputs_sha256':{p.name:sha256(p) for p in output.iterdir() if p.is_file()},'limitations':['These compare geometric image conventions on the same historical simulated structures.','New descriptor semantics do not inherit the old GNN performance metrics.','Paired-state summaries omit undefined contrasts explicitly and report their denominators; this is not negative-label imputation.','The EM-exception omission is a disclosed sensitivity table; main outputs retain every pose.']};atomic_json(output/'summary.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['source_sha256','outputs_sha256','limitations']},indent=2))
