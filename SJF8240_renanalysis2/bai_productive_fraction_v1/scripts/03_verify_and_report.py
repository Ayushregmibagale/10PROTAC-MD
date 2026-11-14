"""Independent Boolean replay,model-count arithmetic and qualified paper comparison."""
from pathlib import Path
import json,hashlib,itertools,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/sjf-bai-mpl')
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];SOURCE=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2/crl2_ubiquitination_v1')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rule(distance,exposure,clash):return distance<=16 and exposure>=.1 and clash<=.8
def pair_auc(labels,values):
 positives=[v for v,y in zip(values,labels) if y==1];negatives=[v for v,y in zip(values,labels) if y==0];return float(np.mean([float(a>b)+.5*float(a==b) for a in positives for b in negatives]))
def main():
 assert rule(16,.1,.8) and not rule(16.001,.1,.8) and not rule(16,.099,.8) and not rule(16,.1,.801)
 scoring=json.loads((ROOT/'validation/scoring_receipt.json').read_text())
 for name,h in scoring['outputs_sha256'].items():assert sha(ROOT/'results'/name)==h
 cases=json.loads((SOURCE/'completion.json').read_text())['records'];poses=pd.read_csv(ROOT/'results/pose_scores.tsv',sep='\t');replayed=[]
 for case in cases:
  loaded=np.load(SOURCE/'results/poses'/case['poi']/case['pose']/'packet.npz');p={key:loaded[key] for key in loaded.files};loaded.close();names=p['spec_names'].tolist();indices=[names.index('closed_observed__VHL'),names.index('open_RBX1_transferred_module__VHL')];productive=[]
  for frame in range(201):
   yes=False
   for si in indices:
    if p['target_scaffold_max_overlap_A'][si,frame]>.8:continue
    for lys in range(len(p['lys_resids'])):
     if rule(float(p['donor_distance_A'][si,frame,1,lys]),float(p['nz_exposed_fraction'][frame,lys]),float(p['target_scaffold_max_overlap_A'][si,frame])):yes=True;break
    if yes:break
   productive.append(yes)
  expected=poses[(poses.poi==case['poi'])&(poses.pose==case['pose'])&(poses['tail']=='tail_fit75')&(poses.distance_A==16)&(poses.exposure_fraction==.1)&(poses.scaffold_scope=='closed_plus_modeled_open')].iloc[0];score=sum(productive)/201;assert abs(score-expected.productive_frame_fraction)<1e-12;assert bool(productive[0])==bool(expected.frame0_productive);replayed.append({'poi':case['poi'],'pose':case['pose'],'productive_frames':int(sum(productive)),'all_frames':201})
 targets=pd.read_csv(ROOT/'results/primary_target_scores.tsv',sep='\t');cluster=pd.read_csv(ROOT/'results/cluster_target_scores.tsv',sep='\t');members=pd.read_csv(ROOT/'results/cluster_membership.tsv',sep='\t');flags=pd.read_csv(ROOT/'results/frame0_pose_flags.tsv',sep='\t');c=cluster[cluster.scaffold_scope=='closed_plus_modeled_open'];auc=pair_auc(targets.label,targets.productive_frame_fraction);static_auc=pair_auc(c.label,c.member_weighted_productive_fraction)
 for row in c.itertuples():
  entries=members[members.poi==row.poi];f=flags[(flags.poi==row.poi)&(flags.scaffold_scope=='closed_plus_modeled_open')].set_index('pose');total=0;equal=[]
  for _,g in entries.groupby('cluster'):
   representative=g.representative_pose.iloc[0];productive=bool(f.loc[representative,'frame0_productive']);total+=len(g)*productive;equal.append(productive)
  assert abs(total/len(entries)-row.member_weighted_productive_fraction)<1e-12;assert abs(sum(equal)/len(equal)-row.cluster_representative_productive_fraction)<1e-12
 verification={'status':'PASS','Boolean_frame_replays':replayed,'cases':141,'frames':28341,'threshold_boundaries':'inclusive16A,exposure>=0.1,clash<=0.8','independent_MD_AUROC':auc,'independent_cluster_member_weighted_AUROC':static_auc,'independent_method':'Explicit scalar Python loops for every pose/frame; pairwise AUROC counting; independent cluster-member sum arithmetic. Geometry reused from independently replayed CRL2 packets.','source_script_sha256':sha(__file__),'not_full_structural_paper_reproduction':True};(ROOT/'validation/independent_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
 text=f'''# Bai-style productive-fraction test — 14 November 2026

The paper's central productive/unproductive rule has now been implemented as a separate VHL adaptation. It does not improve the current eight-target ranking: MD productive-frame fraction AUROC={auc:.3f},versus0.750 for the earlier modified CRL2 score. The closer clustered static-starting-pose fraction gives AUROC={static_auc:.3f}. Neither result is a validated biological classifier or a literal replication of the CRBN paper.

## Rule and denominator

A configuration is productive if at least one target lysine is exposed and≤16Å from ubiquitin's C-terminus,and the target does not clash with the chosen ligase scaffold. This implementation defines the atoms as lysineNZ and modeled ubiquitinGly76 carbonylC,exposure as NZ accessible-sphere fraction≥0.1,and a clash as maximum heavy-atom van-der-Waals overlap>0.8Å. These quantitative atom/exposure/clash definitions are explicit local choices; the available paper methods do not specify them sufficiently to claim exact reproduction.

The extra CE→NZ orientation requirement and the per-frame PROTAC bridging requirement have both been removed. A pose/frame qualifies if either declared VHL scaffold permits it. The main terminal model is fit75;fit74 is a separately reported sensitivity,not a second chance used to increase the score. Genuine clashes remain unproductive and remain in the denominator. No labels,cutoffs or scaffold weights were selected to improve the result.

Productive fraction=number of productive configurations/number of evaluated configurations. A count of lysines,contact atoms,raw docking poses or a docking energy is not this score.

## Results

| Target | User label | MD productive fraction,% | Clustered-starting-pose member-weighted fraction,% |
|---|---|---:|---:|
'''
 for r in targets.itertuples():
  static=c[c.poi==r.poi].iloc[0].member_weighted_productive_fraction;text+=f'| {r.poi} | '+('Degrader' if r.label else 'Non-degrader')+f' | {100*r.productive_frame_fraction:.3f} | {100*static:.3f} |\n'
 primarymetric=pd.read_csv(ROOT/'results/metrics.tsv',sep='\t');m=primarymetric[(primarymetric['tail']=='tail_fit75')&(primarymetric.exposure_fraction==.1)&(primarymetric.distance_A==16)&(primarymetric.scaffold_scope=='closed_plus_modeled_open')&(primarymetric.score=='productive_frame_fraction')].iloc[0]
 text+=f'''\nAt the unvalidated descriptive rule score>0,the MD score has sensitivity{m.sensitivity_at_nonzero:.3f},specificity{m.specificity_at_nonzero:.3f},balanced accuracy{m.balanced_accuracy_at_nonzero:.3f}. Its exact two-sided target-label permutation p={m.exact_two_sided_p:.6f}. There are eight independent labelled targets,not28,341 independent observations. AXL remains nonzero and MAPK14_full remains zero. This is a retrospective geometric ranking,not LOPO-trained GNN performance.

An explicit3Å anchor-fit QC sensitivity produces a different MD AUROC;it is saved in the full table and is not promoted as the new primary result. Closed-only versus closed-plus-modeled-open scaffolds,both terminal placements,exposure cutoffs0.05/0.1/0.2 and16/20Å distances are all retained in the tables.

## Closer correspondence to the paper

For the pose-count analysis,one production-frame0 configuration per retained MD pose was aligned on trueVHL and clustered at3Å target Cα RMSD. Complete-linkage clustering was used,requiring matched sequence/coordinate correspondence. Each cluster's geometric medoid was evaluated;its productive/unproductive classification was propagated to all retained members. Both equal-cluster and member-weighted fractions are reported. This avoids passing off an MD frame percentage as a clustered-pose percentage.

These are the141 retained starting configurations of completed MD,not the paper's complete Rosetta ternary-model population. The available set does not include a comparable complete clustered Rosetta ensemble for all eight proteins. The paper uses MMTSB clustering and a top interface-energy representative;our retained-pose diagnostic uses complete-linkage and a coordinate medoid. That difference is recorded rather than hidden. No new docking was run.

## What still differs from a literal Bai reproduction

- The published system is CRBN/CRL4A;this adaptation is VHL/CRL2 because SJF8240 recruits VHL.
- The local closed template is a trapped mutant with incomplete sidechains. The open catalytic-module placement is modeled.
- DonorGly76 is absent from the deposited catalytic structure;its endpoint remains modeled. The unpassed native-reference calibration from the earlier analysis is not repaired by simply changing the scoring rule.
- Source chemistry and actual-coordinate stereochemistry of historical trajectories remain separate from supplier-SJF8240 equivalence.
- The full original Rosetta model populations,author-selected representative energies and fully specified exposure/clash/atom criteria are not reproduced here.

Therefore this is the paper-style productive-fraction hypothesis tested on our available data,not an exact reproduction or evidence that the published method fails. A faithful VHL analogue needs validated complete catalytic reference geometry and a suitable full ternary ensemble. Experimental ubiquitination data remain the appropriate validation target.

## Saved files

- config/protocol.json:criteria and explicit reproduction gaps.
- results/primary_target_scores.tsv:MD productive fractions for the eight targets.
- results/frame0_pose_flags.tsv and cluster_target_scores.tsv:pose-count comparison.
- results/cluster_membership.tsv:every cluster,representative and member weight.
- results/metrics.tsv and target_scores_all_variants.tsv:all predefined sensitivities.
- validation/independent_verification.json:all141 pose/28,341 frame Boolean checks,pairwise AUROC verification and member-count arithmetic.
- validation/frame0_coordinate_sources.json:source hashes for all frame0 geometries.
- scripts/:reproducible scoring,clustering and verification/report generation.

Original inputs and prior analyses are unchanged. No new MD or GNN training. SIK3/RPS6KB2 remain unknown and are excluded from labelled evaluation.

Source: Bai et al.(2022), https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/ . Exact modeled-scaffold provenance is in ../crl2_ubiquitination_v1/validation/scaffold_preparation.json.
'''
 (ROOT/'FINDINGS.md').write_text(text)
 from docx import Document
 doc=Document()
 for line in text.splitlines():
  if line.startswith('# '):doc.add_heading(line[2:],0)
  elif line.startswith('## '):doc.add_heading(line[3:],1)
  elif line:doc.add_paragraph(line)
 doc.save(ROOT/'FINDINGS.docx');print('PASS:141 poses,28,341 frames; MD AUROC',auc,';clustered representative member-weighted AUROC',static_auc)
if __name__=='__main__':main()
