from pathlib import Path
import json,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/sjf-kinome-mpl')
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
 d=pd.read_csv(ROOT/'pose_checks.tsv',sep='\t');t=pd.read_csv(ROOT/'previous_top10_comparison.tsv',sep='\t');s=pd.read_csv(ROOT/'cutoff_sensitivity.tsv',sep='\t');verification=json.loads((ROOT/'independent_verification.json').read_text());assert verification['status']=='PASS'
 rows=[]
 for r in t.itertuples():rows.append({'previous_rank':r.previous_rank,'gene':r.gene,'uid':r.uid,'correct_head_contacts_4A':'yes' if r.heads_contact_4A else 'no','balanced_contacts_4A':r.balanced_head_contacts_4A,'balanced_contacts_4p5A':r.balanced_head_contacts_4p5A,'balanced_contacts_5A':r.balanced_head_contacts_5A,'maximum_intercomponent_overlap_A':round(r.maximum_intercomponent_overlap_A,3),'internal_linker_overlap_A':round(r.ligand_internal_nonbonded_max_overlap_A,3),'assessment':'Do not advance this pose; clashes and untrusted pocket references','source_pdb':r.path,'source_sha256':r.pdb_sha256})
 compact=pd.DataFrame(rows);compact.to_csv(ROOT/'TOP10_ASSESSMENT.tsv',sep='\t',index=False)
 for tag in ['4A','4p5A','5A']:
  x=d[d['heads_contact_'+tag]&(d.maximum_intercomponent_overlap_A<=.8)].sort_values(['balanced_head_contacts_'+tag,'maximum_intercomponent_overlap_A'],ascending=[False,True]);x.to_csv(ROOT/f'diagnostic_external_only_{tag}.tsv',sep='\t',index=False)
 text='''# Stronger kinome-pose assessment — 13 November 2026

None of the previous top-ten distinct-kinase poses is suitable to advance unchanged on these checks. The previous contact-count ranking should be treated as superseded for selection, not as evidence that these kinases bind or degrade.

## What was checked

All 723 unique poses that passed the original operational screen were checked, across 198 retained source pools. There were no per-pose processing errors. This is a refinement of the previous shortlist, not a fresh search of all rejected, missing or discarded inputs. The old top ten remains preserved for direct comparison.

The saved full-PROTAC Rosetta topology was parsed with RDKit. Both head subgraphs have a unique atom set in the full 79-heavy-atom molecule; graph symmetries do not change the selected atoms. Contacts were measured separately for the correct VHL head and kinase head at 4, 4.5 and 5 Å, requiring three atoms from each head. These measure correct-head proximity to the protein, not independently established pocket engagement.

Intercomponent overlaps use declared element-specific van der Waals radii: C 1.70, N 1.55, O 1.52, F 1.47 and S 1.80 Å, with other supported elements recorded in protocol.json. Overlap is radius1 + radius2 minus distance. Sensitivity thresholds of 0.4, 0.6 and 0.8 Å are operational heavy-atom screens; they are not hydrogen-aware MolProbity clashscores, and no hydrogen-bond exceptions or energy minimization were applied. Large carbon/oxygen overlaps remain explicit concerns. Both reconstructed-full-E3 and observed-fragment protein–protein clashes are recorded separately.

## Why the conclusion changed

1. All ten earlier leaders have maximum intercomponent overlaps above 0.8 Å. The original fixed 1.5 Å minimum-distance cutoff admitted these clashes and favored crowded poses through contact counts.
2. All 723 poses contain an internal PROTAC nonbonded overlap above 0.8 Å under the saved topology. O11–C16, five bonds apart, are approximately 2.03 Å apart in the previous top ten, corresponding to approximately 1.19 Å of overlap. This is not a directly bonded or 1–3-neighbor exception. Covalent bond lengths themselves pass the broad 1.0–2.1 Å integrity check. This finding concerns stored static coordinates; existing MD trajectories were not revalidated here.
3. Both PT0 and PT1 reference-head geometries fail the topology-based bond-length check in all 198 reviewed pools. For example, TXK reference head bonded pairs extend to 8.48 Å for PT0 and 35.22 Å for PT1. The provenance may involve atom-name/coordinate handling, but the root cause has not been established. Consequently the old reference-derived pocket regions and any RMSD against those references cannot establish correct pocket engagement. The new script refuses to claim otherwise. Reconstructing full VHL does not repair the head references.

'''
 text+='## Comparison of the original top ten\n\n| Kinase | Correct-head contacts at 4 Å (at least 3 each) | Largest intercomponent overlap, Å | Decision |\n|---|---|---:|---|\n'
 for r in compact.itertuples():text+=f'| {r.gene} | {r.correct_head_contacts_4A} | {r.maximum_intercomponent_overlap_A:.3f} | Do not advance this pose |\n'
 text+='\n## Cutoff sensitivity\n\n| Contact cutoff, Å | Correct-head proximity only, poses | External overlap ≤0.6 Å, poses | External overlap ≤0.8 Å, poses | All integrity checks, poses |\n|---|---:|---:|---:|---:|\n'
 for cutoff,tag in [(4,'4A'),(4.5,'4p5A'),(5,'5A')]:
  h=d['heads_contact_'+tag];text+=f'| {cutoff} | {int(h.sum())} | {int((h&(d.maximum_intercomponent_overlap_A<=.6)).sum())} | {int((h&(d.maximum_intercomponent_overlap_A<=.8)).sum())} | 0 |\n'
 text+='''
The empty qualified top-ten files are deliberate: the criteria were not relaxed to manufacture ten winners. PAK4 and IGF1R pass correct-head proximity at 4 Å with the lenient external-overlap threshold, but fail the shared internal geometry check and lack trustworthy reference-pocket validation. They are diagnostic repair candidates only, not approved binders or MD priorities.

## Recommended next action

Repair or regenerate the PT0/PT1 reference coordinates with chemically verified atom correspondence, verify native VHL-ligand and kinase-warhead interactions against trustworthy reference structures, and repair/minimize the shared PROTAC geometry before reranking. Preserve the current results as an audit trail. Repeat the filters after repair and use ensemble/MD stability only on geometrically credible poses. This requires a separate structure-repair stage; no docking, molecule repair, new MD or GNN fit was performed in this comparison. Failure of a pose is not failure of the kinase to bind, and no experimental degradation labels are inferred for prospective targets.

## Saved evidence and reproducibility

- TOP10_ASSESSMENT.tsv: compact decisions and exact source paths/hashes.
- previous_top10_comparison.tsv: detailed old-versus-new pose metrics.
- pose_checks.tsv: all 723 pose checks.
- cutoff_sensitivity.tsv: all nine cutoff/overlap combinations.
- head_mapping_receipts.json: chemical atom maps, reference failures and input hashes.
- protocol.json and completion.json: parameters, accounting and output hashes.
- independent_verification.json: independent raw-PDB/NumPy replay of the ten kinase overlap measurements, shared linker distances, threshold arithmetic and output checksums; PASS.
- scripts/check_engagement.py, scripts/verify_results.py and scripts/write_report.py: reproducible analysis, independent checks and report generation.

The overlap methodology is informed by structural-validation principles, but this custom heavy-atom audit is not MolProbity: https://pmc.ncbi.nlm.nih.gov/articles/PMC1933162/ . Full hydrogen placement and chemically informed contact validation are a subsequent option, not a result claimed here.
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
 fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained');yy=np.arange(len(t));axes[0].barh(yy,t.balanced_head_contacts_4A,color='#367da2');axes[0].axvline(3,color='black',ls='--',lw=1);axes[0].set_yticks(yy,t.gene);axes[0].invert_yaxis();axes[0].set_xlabel('Smaller correct-head contact count at 4 Å');axes[0].set_title('Contact extent')
 axes[1].barh(yy,t.maximum_intercomponent_overlap_A,color='#b95850');axes[1].axvline(.8,color='black',ls='--',lw=1,label='Lenient overlap limit');axes[1].set_yticks(yy,t.gene);axes[1].invert_yaxis();axes[1].set_xlabel('Largest intercomponent vdW overlap (Å)');axes[1].set_title('All original top-ten poses fail');axes[1].legend(loc='lower right');fig.suptitle('SJF8240 kinome: previous top-ten pose reassessment');fig.savefig(ROOT/'top10_comparison.png',dpi=180);fig.savefig(ROOT/'top10_comparison.pdf');plt.close(fig)
 print('Saved findings, Word document, comparison tables and figures.')
if __name__=='__main__':main()
