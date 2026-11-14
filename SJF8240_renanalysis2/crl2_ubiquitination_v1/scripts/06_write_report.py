"""Write results with software success distinct from mechanistic validation."""
from pathlib import Path
import json,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/sjf-crl2-mpl')
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
 independent=json.loads((ROOT/'validation/independent_verification.json').read_text());assert independent['status']=='PASS';control=json.loads((ROOT/'validation/controls.json').read_text());targets=pd.read_csv(ROOT/'results/primary_target_scores.tsv',sep='\t');metrics=pd.read_csv(ROOT/'results/metrics.tsv',sep='\t');primary=metrics[(metrics.scaffold=='closed_observed')&(metrics.anchor=='VHL')];new=primary[primary.score=='productivity_joint_fraction'].iloc[0];text='''# Full-CRL2 VHL compatibility test — 14 November 2026

The primary retrospective target-ranking AUROC is0.750, compared with0.333 for true-VHL bridging fraction and0.583 for the prior nearest-E2 reach50 fraction. This is a new geometric compatibility analysis, not the earlier temporal GNN result that also had a0.750 matched-pair AUROC. There was no GNN training, no LOPO model fitting and no new MD.

All141 retained MD poses from the eight user-labelled targets were analyzed on the same201-frame grid0–20ns (28,341 frames). Six labels are degrader and two are non-degrader. SIK3/RPS6KB2 were not used as negatives or validation. Each target score is an equal mean over poses; frames,lysines and reference alternatives are not independent experimental examples.

## Prespecified primary score

The primary scaffold is8RX0, fitted through sequence-matched true-VHL Cα atoms. For each frame and each of two modeled ubiquitin terminal placements, at least one lysine must be within20Å of modeled donor Gly76 carbonylC, have NZ accessible-sphere fraction≥0.10, and have CE→NZ pointing into the donor-facing hemisphere (cosine≥0). Target–new-scaffold heavy-atom overlap must be≤0.8Å and anchor RMSD≤3Å. The primary joint score also requires simultaneous PROTAC/true-VHL and PROTAC/target heavy contacts<5Å. This is still an any-ligand contact proxy, not atom-mapped cognate-head binding validation.

The scalar is the fraction of sampled frame×endpoint combinations satisfying all these conditions, averaged equally over retained poses. It is not a thermodynamic probability, catalytic rate, cellular degradation probability or measured percentage degradation. Genuine target–scaffold clashes remain in the denominator as failures. Input QC failures would be reported separately rather than called biological negatives.

| Target | Label | Poses | Bridge fraction,% | Joint compatibility,% |
|---|---|---:|---:|---:|
'''
 for row in targets.itertuples():text+=f'| {row.poi} | '+('Degrader' if row.label else 'Non-degrader')+f' | {row.poses} | {100*row.bridge_fraction:.2f} | {100*row.productivity_joint_fraction:.3f} |\n'
 text+='\n## Ranking and operational classifications\n\n| Score | Target AUROC |\n|---|---:|\n'
 for score in ['bridge_fraction','reach50_fraction','distance_only_fraction','exposed_reachable_fraction','oriented_reachable_fraction','geometry_fraction','productivity_joint_fraction','conditional_productivity_fraction']:
  row=primary[primary.score==score].iloc[0];text+=f'| {score} | {row.AUROC:.3f} |\n'
 text+=f'''\nThe joint-score ranking wins9 of12 positive–negative target pairs when ties count0.5. Its exact two-sided AUROC permutation p-value is{new.exact_two_sided_AUROC_permutation_p:.6f} (all28 allocations of two negative labels). That is not convincing statistical evidence with eight targets. At the explicitly descriptive, unvalidated rule score>0, sensitivity is{new.sensitivity_descriptive:.3f}, specificity{new.specificity_descriptive:.3f} and balanced accuracy{new.balanced_accuracy_descriptive:.3f}. The score misses MAPK14_full and RIPK2_full and gives AXL a nonzero score. No decision threshold was tuned to improve classification.

Conditional-on-bridge scores are separately reported. The0.875 geometry-only AUROC is an ablation without the joint bridging condition, not the prespecified primary result. Its performance is not grounds to switch the primary endpoint after seeing the labels.

## Reference construction and limitations

8RX0 is a trapped UBE2R1-mutant CRL2 complex with incompletely resolved sidechains. Donor ubiquitin stops at Gly75; Gly76 is not deposited. Complete ubiquitin1UBQ was locally fitted using terminal backbone correspondences. Fits starting70,72 and73 generated inconsistent peptide distances and were rejected. Fits74 and75 preserved G75C–G76N within1.1–1.6Å and were retained as two explicitly modeled endpoints, not observed catalytic atoms. The initial failed attempt remains archived. The20Å primary cutoff was chosen as a coarse hypothesis informed by the VHL literature, not copied from the CRBN paper's16Å criterion or selected from these target results.

An open-scaffold sensitivity uses8RWZ with the closed catalytic module transferred through the shared RBX1 RING region40–104 (42 common Cα atoms,RMSD2.899Å). An initial all-RBX1 fit was poor (5.310Å); domain-based fitting was defined and frozen before target outcome scoring. This alternative is a constructed model, not another deposited active catalytic state. The source file and preparation receipts expose this distinction. Missing catalytic-module sidechains were not invented, so the steric assessment is an incomplete-atom lower bound. No energetic or reaction-pathway meaning is assigned to the orientation heuristic.

Scaffold/anchor sensitivity matters: at the same score settings,the open model gives AUROC0.875 with VHL anchoring but0.583 with ELOC anchoring. The closed model gives0.750 for either anchor. All120 predefined distance/exposure/orientation/scaffold/anchor combinations are saved; no maximum is promoted as validation.

## Native control and scientific release

Analytic controls for rigid transforms,inverse transforms,isolated/buried NZ exposure and a translated far negative pass. Independent replays using MDAnalysis unwrap,Bio.SVDSuperimposer and dense SciPy distances reproduce six raw frames from ABL1/DDR2,including NZ exposure and target–scaffold overlap. Target AUROCs were independently verified by counting positive–negative pairs. All141 input trajectories/topologies/PDBs match the prior source-bound audit and all28,341 frame contacts replay within0.002Å of the independent true-VHL contact audit.

However,the native MZ1 structural projection control does not recover a close,correctly oriented known positive under the primary20Å/hemisphere rule. Projected K456 distances are about22.5–25.0Å to the two modeled endpoints,while other reported sites may be close but have unfavorable frozen rotamers. This is a calibration/representation limitation,not a software PASS. Static5T35 rotamers and the missing donor tail are not a catalytic transition-state model. This test therefore does not release a validated degradation filter,does not justify promoting the new score to a GNN feature as proven biology,and does not approve a new MD shortlist.

Historical molecule identity/stereochemistry limitations remain: these are the actual archived trajectories,not supplier-equivalent SJF8240 validation. Only the retrospective structural hypothesis was tested. Binding-site provenance and correct-head placement remain separate gates.

## Interpretation and next step

This is a promising ranking signal worth investigating,not reliable separation of binders that degrade from binders that do not. The next technical improvement is to validate/reconstruct the complete E2~Ub catalytic endpoint and sidechain context against an appropriate native control,then repeat the frozen target comparison. Do not tune endpoint position,scaffold weights,cutoffs or thresholds to force ABL1/AXL below the six positives. Broader GNN or400-kinase screening should wait for those geometric controls or clearly remain exploratory.

## Sources

- Bai et al.(2022),CRL4A ubiquitination-ensemble concept: https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/
- Crowe et al.(2024),VHL catalytic machinery and ubiquitination-zone evidence: https://pmc.ncbi.nlm.nih.gov/articles/PMC11468923/
- Deposited coordinates: https://www.rcsb.org/structure/8RX0 ; https://www.rcsb.org/structure/8RWZ ; https://www.rcsb.org/structure/5T35 ; https://www.rcsb.org/structure/1UBQ

Source/configuration hashes,reference preparation,controls,per-pose packets,per-frame primary scores,target scores,all sensitivities and independent checks are saved separately. No original trajectory or prior result was overwritten.
'''
 (ROOT/'FINDINGS.md').write_text(text)
 from docx import Document
 doc=Document()
 for line in text.splitlines():
  if line.startswith('# '):doc.add_heading(line[2:],0)
  elif line.startswith('## '):doc.add_heading(line[3:],1)
  elif line:doc.add_paragraph(line)
 doc.save(ROOT/'FINDINGS.docx')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained');yy=np.arange(len(targets));colors=['#d08059' if label==0 else '#447fa3' for label in targets.label];axes[0].barh(yy,100*targets.bridge_fraction,color=colors);axes[1].barh(yy,100*targets.productivity_joint_fraction,color=colors)
 for ax in axes:ax.set_yticks(yy,targets.poi);ax.invert_yaxis()
 axes[0].set_xlabel('Mean bridged frames (%)');axes[0].set_title('True-VHL bridging: AUROC0.333');axes[1].set_xlabel('Joint geometric compatibility (%)');axes[1].set_title('Primary CRL2 score: AUROC0.750');fig.suptitle('Retrospective eight-target test; native control not validated');fig.savefig(ROOT/'results/target_comparison.png',dpi=180);fig.savefig(ROOT/'results/target_comparison.pdf');plt.close(fig)
 (ROOT/'validation/scientific_status.json').write_text(json.dumps({'extraction_complete':True,'software_numerical_verification':'PASS','native_primary_distance_orientation_control_pass':bool(control['known_positive_primary_distance_orientation_pass']),'biological_validation':False,'validated_degradation_classifier':False,'automatic_GNN_training_release':False,'primary_AUROC':float(new.AUROC),'exact_two_sided_AUROC_permutation_p':float(new.exact_two_sided_AUROC_permutation_p),'historical_simulated_chemistry':True},indent=2)+'\n');print('Saved complete report; primary AUROC0.750; scientific release remains false.')
if __name__=='__main__':main()
