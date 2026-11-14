"""Summarize completed native controls,target comparison,GNN readiness and kinome audit."""
from pathlib import Path
import json,os,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/sjf-native2-mpl')
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 verify=json.loads((ROOT/'validation/independent_verification.json').read_text());assert verify['status']=='PASS';gnn=json.loads((ROOT/'validation/gnn_runtime_check.json').read_text());assert gnn['status']=='PASS';native=pd.read_csv(ROOT/'results/native_controls.tsv',sep='\t');target=pd.read_csv(ROOT/'results/primary_target_scores.tsv',sep='\t');metrics=pd.read_csv(ROOT/'results/native_metrics.tsv',sep='\t');primary=metrics[(metrics.reference=='8R5H_observed_SY8')&(metrics.anchor=='VHL')&(metrics.score=='productive_fraction')].iloc[0];kin=json.loads((ROOT/'validation/kinome_receipt.json').read_text());assert len(target)==8 and kin['poses']==723;known_sites=native[native.lysine.isin([368,445,456])];known20=bool(known_sites.primary20A_geometry.any());r1=native[(native.template=='8RX0_G76_fit75')&(native.lysine==456)];r2=native[(native.template=='8R5H_observed_SY8')&(native.lysine==445)]
 status={'native_reference_geometric_control_20A':'PASS' if known20 else 'FAIL','native_backbone_correspondence_corrected':True,'exact_published18p4A_reproduction':False,'literal_Bai_CRBN_reproduction':False,'wildtype_ubiquitin_carbonyl_resolved':False,'primary_AUROC':float(primary.AUROC),'primary_exact_target_permutation_p':float(primary.exact_two_sided_AUROC_p),'software_independent_verification':'PASS','GNN_graphs_prepared':16,'GNN_finite_forward':'PASS','GNN_training_run':False,'GNN_optimization_steps':0,'new_MD':False,'kinome_annotation_poses':723,'kinome_geometry_compatible':kin['raw_geometry_productive'],'kinome_scientifically_eligible':kin['scientifically_eligible'],'validated_degrader_classifier':False,'original_inputs_modified':False,'automatic_production_release':False};(ROOT/'validation/scientific_status.json').write_text(json.dumps(status,indent=2)+'\n')
 text=f'''# Native catalytic-reference validation and follow-through — 14 November 2026

All four requested workstreams were assessed: native control/reference correction,reference-ensemble checks,the full eight-target MD comparison,and GNN-input/wider-kinome assessment. The native structural check improves,but the independently sourced primary descriptor does not separate the eight targets: AUROC{primary.AUROC:.3f}. No new biological classifier is released.

## 1. Correction to the earlier native-control interpretation

The earlier control superposed crystallographic5T35 BRD4 through VHL and carried over the crystal distal backbone/sidechain placement. That is not the target geometry present in the catalytic cryo-EM model. It made K456 appear23–25Å from the modeled donor endpoint. The corrected control keeps the deposited catalytic BRD4 backbone and transfers only missing sidechain heavy atoms through residue-local N/CA/C fits,with matching amino acids and RMSD≤0.25Å. Existing observed sidechain coordinates remain untouched. Crystal copiesA/E are retained as a sensitivity. No rotamers were optimized toward the donor or the eight labels.

After that correction,K456 in8RX0 is approximately{r1.NZ_distance_A.min():.2f}–{r1.NZ_distance_A.max():.2f}Å from the prior modeled Gly76 endpoint,with exposed NZ. This is approximate consistency with the published18.4Å mechanistic illustration,not an exact numerical reproduction: our missing-tail and sidechain models remain distinct. The original native-control failure must therefore be qualified as a representation/correspondence problem,not evidence that the native VHL system is unproductive.

The old hemisphere orientation rule still rejects these frozen K456 rotamers,so it is not reinstated as a mandatory reaction criterion. A static lysine rotamer is not an experimentally established transition-state attack geometry.

## 2. Independently sourced catalytic-reference test

We downloaded8R5H from the Li et al.(2024) VHL–MZ1–BRD4–UBE2R2 study and its SY8 chemical-component definition. It contains a chemically trapped three-way junction. PDB connectivity explicitly links E2Cys93SG to SY8C3 and UbGly75C to SY8N. CCD bonds give the path N–C4–C3. Thus C3 is a deposited,protein-tethered reaction-center surrogate;it is not native ubiquitinGly76 carbonylC. SY8's oxygen is attached elsewhere in the artificial linker. We do not silently relabel this as a resolved physiological thioester.

The C93SG–SY8C3 distance is approximately1.705Å. With the native observed scaffold,target and missing-sidechain procedure,the known K445 site is{r2.NZ_distance_A.min():.2f}Å from that surrogate,exposed,and within the coarse20Å geometry criterion. A500Å-translated target fails reach. The trapped K368C site is retained asCYS and never counted as a wild-type lysine.

Three declared references were constructed before eight-target rescoring: deposited8R5H/SY8,8RX0 with the prior modeled Gly76 fit75,and8RX0 with SY8 transferred through the common E2 active-loop backbone85–100. The active-loop fit is0.558Å. No endpoints were fitted to target degradation outcomes. Two anchors,VHL and ELOC,are compared. Missing atoms/components remain explicit:8R5H has no resolvedNEDD8 and8RX0 has multiple backbone-only catalytic components. These are not complete atomic reaction simulations.

## 3. Full eight-target MD result

All141 retained cases and201 saved frames each were processed:28,341 frames. The primary reference/anchor was frozen as8R5H/SY8 with VHL anchoring. The primary rule is: at least one target lysineNZ within20Å of C3,accessible-sphere fraction≥0.10,and maximum target/new-scaffold heavy-atom overlap≤0.8Å. No orientation or bridging condition is imposed on that primary productive fraction. Joint-bridge and fit3Å variants are separate sensitivities. Target-induced clashes remain in the denominator;frames and lysines are not independent labelled observations.

| Target | User label | Primary productive fraction,% | Joint bridged/productive fraction,% |
|---|---|---:|---:|
'''
 for r in target.itertuples():text+=f'| {r.poi} | '+('Degrader' if r.label else 'Non-degrader')+f' | {100*r.productive_fraction:.3f} | {100*r.joint_bridge_productive_fraction:.3f} |\n'
 text+=f'''\nPrimary target-ranking AUROC={primary.AUROC:.3f};exact two-sided permutation p={primary.exact_two_sided_AUROC_p:.3f} across all28 allocations of two negative labels. All eight targets have nonzero scores. The unvalidated score>0 classifier therefore has sensitivity1,specificity0,balanced accuracy0.5. ABL1 is geometrically compatible despite being a user-confirmed non-degrader. No threshold was tuned to remove it.

| Reference | Anchor | AUROC of the same productive-fraction rule |
|---|---|---:|
'''
 for r in metrics[metrics.score=='productive_fraction'].itertuples():text+=f'| {r.reference} | {r.anchor} | {r.AUROC:.3f} |\n'
 text+='''
The earlier0.750 result remains a different exploratory analysis with different endpoint/scaffold choices and additional rules;it was not validated by the new controls. Some alternative reference/anchor combinations score higher,but selecting those after examining eight labels would overstate performance. The stable conclusion is reference dependence and inadequate current evidence for a useful degradation classifier. This does not prove that full-machinery physics cannot help,or that a GNN cannot use finer patterns.

The exposure is the NZ accessible-sphere fraction from180 Fibonacci points and1.4Å probe,not an experimental SASA measurement. Overlaps use declared heavy-atom radii and omit missing sidechain atoms;they are not MolProbity scores. Distances to the SY8 surrogate are not catalytic probabilities.

## 4. GNN input assessment

Prepared16 matched four-window graphs:eight with five provisional native-reference lysine-node descriptors and eight zero-padded controls. Existing29 node features and every edge are byte-identical between controls and augmented inputs;both have34 features. Added columns are mean NZ–SY8C3 distance/100Å,reachable20Å fraction,exposure fraction,exposed/reachable/scaffold-compatible fraction,and joint bridge/productive fraction. These are explicitly engineered descriptors,not new full-scaffold nodes or observed ubiquitination labels. Each new descriptor uses equal pose weighting over the full201-frame data,while the original base features remain unchanged.

The existing TemporalEncoder was software-tested with only the first GAT input width changed29→34,using identical34-input architecture forboth variants. All16 graphs pass finite untrained CPU forward checks;there were zero optimizer steps. An initial unguardedPyG import hit an incompatible optional native extension;the existing process-local native-PyG guard resolves it without modifying installed packages. No trained predictions or replacement weights were generated.

This completes input/runtime readiness assessment,not GNN validation. A new training campaign is not automatically justified by a descriptor whose independent primary ranking is chance-level and whose historical chemistry remains unresolved. If a controlled exploratory GNN ablation is subsequently pursued,it must hold out entire proteins,retain the same architecture/splits/seeds/epochs forboth variants,and report all folds and CPU replay. Eight labelled targets cannot establish generalization regardless of graph/frame count.

## 5. Wider kinome assessment

'''
 text+=f'''Annotated all{kin['poses']} original QC-shortlist poses across{kin['accessions']} accessions with the deposited native-reference surrogate. {kin['raw_geometry_productive']} pass the raw exposed/reachable/no-added-scaffold-clash test. They are hypotheses about geometry only,not predicted binding or degradation. No unknown target was assigned a negative experimental label.

There are{kin['scientifically_eligible']} scientifically eligible new-MD candidates in this unrepaired archive. All723 inherit the historical molecule-connectivity mismatch and shared internal linker-overlap problem;original head-pocket placement validation is unresolved. New compatibility annotation cannot cure those defects. The existing ten-pose local repairs were a separate bounded historical-model diagnostic and are not substituted into the723-pose population. The top-ten table is labelled DIAGNOSTIC_ONLY and must not become a production shortlist.

## Verification and saved evidence

- Native PDB/CCD connectivity,local sidechain-transfer receipts,and every native lysine distance/exposure/clash are saved.
- All141 MD source trajectories/topologies/reference PDBs were rehashed against the earlier source-bound audit. Contact and frame-grid replay covers all28,341 frames.
- Independent MDAnalysis-unwrapping/Bio.SVDSuperimposer/dense-distance replay checks six raw frames from ABL1/DDR2 for distances,overlaps and exposure. Four target AUROCs are independently verified by pair counting.
-16 matched GNN graphs retain original base tensors and edges;finite untrained forward checks pass.
-723 archived kinome source PDBs were hash-checked and sequence-role checked by the existing parser/context contracts.
- Original structures,trajectories,labels,models and previous outputs remain unchanged. No new docking,MD production or GNN training was run.

## Decision

The native-control correspondence issue is corrected and qualitative20Å reference access is recovered. The full native-reference rescoring does not separate the eight labels. The GNN data path is prepared but untrained,and the wider archive remains unsuitable for supplier-SJF8240 selection. The justified next scientific work is an independent assay/chemical-provenance review and a validated candidate ensemble,not further tuning of geometric thresholds on these same eight labels. The chemistry-correct and historical-model branches remain separate. No general claim about binder-versus-degrader prediction is supported.

## Sources

- Crowe et al.(2024),VHL catalytic architecture: https://pmc.ncbi.nlm.nih.gov/articles/PMC11468923/
- Li et al.(2024),VHL/UBE2R2 substrate priming: https://pmc.ncbi.nlm.nih.gov/articles/PMC10997478/
- Bai et al.(2022),productive-fraction concept: https://pmc.ncbi.nlm.nih.gov/articles/PMC9019245/
- Deposited8R5H andCCD: https://www.rcsb.org/structure/8R5H ; https://files.rcsb.org/ligands/view/SY8.cif
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
 fig,ax=plt.subplots(figsize=(9,5),layout='constrained');ax.barh(target.poi,100*target.productive_fraction,color=['#c77c54' if y==0 else '#497da2' for y in target.label]);ax.invert_yaxis();ax.set_xlabel('Native-reference geometric productive fraction (%)');ax.set_title('Independent reference: eight-target AUROC0.500;not a degradation classifier');fig.savefig(ROOT/'results/native_target_comparison.png',dpi=180);fig.savefig(ROOT/'results/native_target_comparison.pdf');plt.close(fig)
 print('Saved all-workstream findings;native control geometry recovered;primary AUROC0.500;no scientific production release.')
if __name__=='__main__':main()
