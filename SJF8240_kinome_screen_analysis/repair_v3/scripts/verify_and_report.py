"""Replay repaired SDF/PDB geometry and write conservative before/after assessment."""
from pathlib import Path
import os,sys,json
sys.dont_write_bytecode=True;os.environ['OPENBLAS_NUM_THREADS']='1';os.environ.setdefault('MPLCONFIGDIR','/tmp/sjf-kinome-mpl')
import numpy as np,pandas as pd
from rdkit import Chem
BASE=Path(os.environ.get('KINOME_ANALYSIS_ROOT',str(Path(__file__).resolve().parents[2])));ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'engagement_v2/scripts'));import check_engagement as c
sys.path.insert(0,str(Path(__file__).resolve().parent));import repair_legacy_poses as repair;import rebuild_references as rebuild

def main():
 src=ROOT/'stronger_relaxation';receipts=json.loads((src/'legacy_repair_receipts.json').read_text());old=pd.read_csv(BASE/'engagement_v2/previous_top10_comparison.tsv',sep='\t');pools=json.loads((BASE/'config/pools.json').read_text());backend=c.load_backend(BASE);cfg=json.loads((BASE/'config/protocol.json').read_text());checks=[];assessments=[]
 for rr in receipts:
  gene=rr['gene'];out=src/'legacy_repairs'/gene;source=Path(rr['path']);c.require(c.sha(source)==rr['pdb_sha256'],'Source data changed');mol=Chem.MolFromMolFile(str(out/'repaired_historical_ligand.sdf'),removeHs=True);c.require(mol is not None,'SDF parsing failed');x=mol.GetConformer().GetPositions();names=json.loads((out/'atom_name_map.json').read_text());st=Chem.Mol(mol);Chem.AssignAtomChiralTagsFromStructure(st,replaceExistingTags=True);c.require(Chem.MolToInchiKey(st)==rr['archived_coordinate_inchikey'],'Chemical graph or3D stereo changed');c.require(len(x)==79,'Wrong ligand atom count')
  pool=next(p for p in pools if p['uid']==rr['uid'] and p['pool']==rr['pool']);context=backend.model_context(pool,cfg);d,g,poi,vhl,vfull,protein,pels=repair.make_complex(source,context,backend,cfg);elements=[a.GetSymbol().upper() for a in mol.GetAtoms()];radii=np.array([c.RADII[e] for e in elements]);rp=np.array([c.RADII[e] for e in pels]);dist=np.sqrt(np.sum((x[:,None,:]-protein[None,:,:])**2,axis=2));external=float((radii[:,None]+rp[None,:]-dist).max());graph=Chem.GetDistanceMatrix(mol);mask=np.triu(graph>2,1);internal_dist=np.sqrt(np.sum((x[:,None,:]-x[None,:,:])**2,axis=2));internal=float((radii[:,None]+radii[None,:]-internal_dist)[mask].max());bondlength=np.array([np.linalg.norm(x[b.GetBeginAtomIdx()]-x[b.GetEndAtomIdx()]) for b in mol.GetBonds()]);c.require(abs(external-rr['external_max_overlap_A'])<.001,'External numeric replay failed');c.require(abs(internal-rr['internal_max_overlap_A'])<.001,'Internal numeric replay failed');n_e3=len(vfull)+len(g['elob']['xyz'])+len(g['eloc']['xyz']);e3=protein[:n_e3];ep=protein[n_e3:];er=rp[:n_e3];pr=rp[n_e3:];interdist=np.sqrt(np.sum((e3[:,None,:]-ep[None,:,:])**2,axis=2));protein_overlap=float((er[:,None]+pr[None,:]-interdist).max());previous=old[old.gene==gene].iloc[0];c.require(abs(protein_overlap-previous.E3_POI_max_overlap_A)<1e-9,'Protein arrangement changed')
  # Write a complete, explicitly role-mapped heavy-atom derivative; no raw files changed.
  lines=[];serial=0
  for group,coords,chain in [(g['elob'],g['elob']['xyz'],'A'),(g['eloc'],g['eloc']['xyz'],'B'),(vhl,vfull,'C'),(poi,poi['xyz'],'D')]:
   residue_names={r:n for r,n,_ in group['ca']}
   for name,resid,element,p in zip(group['names'],group['resids'],c.protein_elements(group),coords):
    serial+=1;xx,yy,zz=p;lines.append(f'ATOM  {serial:5d} {name:>4s} {residue_names[int(resid)]:>3s} {chain}{int(resid):4d}    {xx:8.3f}{yy:8.3f}{zz:8.3f}  1.00  0.00          {element:>2s}  ')
   lines.append('TER')
  for i,p in enumerate(x):
   serial+=1;xx,yy,zz=p;lines.append(f'HETATM{serial:5d} {names[str(i)]:>4s} PTC X   1    {xx:8.3f}{yy:8.3f}{zz:8.3f}  1.00  0.00          {elements[i]:>2s}  ')
  lines.append('END');complexpath=out/'repaired_full_E3_complex.pdb';complexpath.write_text('\n'.join(lines)+'\n');parsed=backend.parse(complexpath);c.require(len(parsed['ligand'])==79 and len(parsed['groups'])==4,'Written full complex invalid');roleca={group['key'][1]:len(group['ca']) for group in parsed['groups']};c.require(roleca['C']==140,'True VHL incomplete')
  # Verify all protein coordinates match frozen originals after PDB rounding.
  desired=[g['elob']['xyz'],g['eloc']['xyz'],vfull,poi['xyz']]
  for chain,expected in zip('ABCD',desired):c.require(np.allclose(next(z['xyz'] for z in parsed['groups'] if z['key'][1]==chain),expected,atol=.000501,rtol=0),'Fixed protein moved')
  ligand_ok=bool(rr['convergence_ok'] and rr['stereo_preserved'] and bondlength.min()>=1. and bondlength.max()<=2.1 and external<=.6 and internal<=.6);protein_ok=protein_overlap<=.6;row={'gene':gene,'uid':rr['uid'],'before_external_overlap_A':previous.maximum_intercomponent_overlap_A,'after_ligand_external_overlap_A':external,'before_internal_overlap_A':previous.ligand_internal_nonbonded_max_overlap_A,'after_internal_overlap_A':internal,'fixed_E3_POI_overlap_A':protein_overlap,'converged':rr['convergence_ok'],'same_historical_graph_and_stereo':True,'bond_min_A':float(bondlength.min()),'bond_max_A':float(bondlength.max()),'maximum_heavy_atom_displacement_A':rr['max_displacement_A'],'ligand_geometry_pass_06':ligand_ok,'fixed_protein_geometry_pass_06':protein_ok,'supplier_SJF8240':False,'pocket_engagement_validated':False,'corrected_MD':False,'repaired_complex':str(complexpath),'repaired_complex_sha256':c.sha(complexpath)}
  for tag in ['4A','4p5A','5A']:row['balanced_head_contacts_'+tag]=rr['balanced_correct_head_contacts_'+tag];row['historical_geometry_candidate_'+tag]=bool(ligand_ok and protein_ok and rr['balanced_correct_head_contacts_'+tag]>=3)
  row['assessment']='Historical geometry candidate only; chemical identity and pocket validation unresolved' if row['historical_geometry_candidate_4A'] else 'Do not advance; residual geometry or correct-head contact failure'
  assessments.append(row);checks.append({'gene':gene,'independent_external_overlap_A':external,'independent_internal_overlap_A':internal,'independent_fixed_protein_overlap_A':protein_overlap,'source_pdb_sha256':c.sha(source),'repaired_sdf_sha256':c.sha(out/'repaired_historical_ligand.sdf'),'fixed_proteins_verified':True,'gradient_check_max_absolute_error':max(z['absolute_error'] for z in rr['gradient_verification'])})
 a=pd.DataFrame(assessments);a.to_csv(ROOT/'REPAIR_COMPARISON.tsv',sep='\t',index=False);rebuildrows=pd.read_csv(ROOT/'reference_rebuild_summary.tsv',sep='\t');reference_evidence=json.loads((ROOT/'reference_rebuild_receipts.json').read_text());c.require(len(rebuildrows)==198 and rebuildrows.vhl_reference_written.all(),'Incomplete VHL repair output')
 for r in reference_evidence:
  p=ROOT/'supplier_reference_candidates'/r['pool']/'vhl_head.sdf';m=Chem.MolFromMolFile(str(p),removeHs=True);rebuild.validate_head(m);c.require(m.GetNumAtoms()==30,'Wrong VHL head');c.require(c.sha(p)==r['files_sha256']['vhl_head.sdf'],'Reference drift')
 refca=[a for a in rebuild.pdb_atoms(rebuild.REF/'protein.pdb') if a['name']=='CA'];head=Chem.MolFromMolFile(str(rebuild.REF/'poi_head.sdf')).GetConformer().GetPositions();_,_,fit=rebuild.placement.align_domain(refca,refca,head);c.require(fit['pocket_rmsd_A']<1e-8,'Alignment positive control failed')
 verification={'status':'PASS','checked_legacy_repairs':len(checks),'VHL_heads_checked':198,'kinase_alignment_self_control_pocket_RMSD_A':fit['pocket_rmsd_A'],'independent_checks':checks,'limitations':'Numeric replay and chemical identity/coordinate checks, not binding/degradation validation. Reuses the reviewed protein-role reconstruction.','sources_sha256':{'verification_script':c.sha(__file__),'repair_script':c.sha(Path(__file__).parent/'refine_legacy_poses.py'),'reference_script':c.sha(Path(__file__).parent/'rebuild_references.py')}};(ROOT/'independent_repair_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
 selected=a[a.historical_geometry_candidate_4A].gene.tolist();identity=json.loads((ROOT/'chemical_identity_audit.json').read_text());text=f'''# Kinome geometry repairs — 14 November 2026

Repairs and verification are complete for the bounded ten-pose diagnostic and VHL references. This does not complete a validated SJF8240 kinome reranking.

## Completed repairs

- Preserved every original docking PDB, trajectory, label and analysis output.
- Performed two recorded local repair stages for the original ten leading poses: MMFF94s plus explicit positional restraints and overlap penalties, fixed proteins, deterministic seed42 and at most two starting perturbations per stage. No docking or MD production was run. The second stage increased steric penalties without changing the evaluation gates. Its finite-difference gradient checks passed. Optimization energies are repair objectives, not binding energies or cross-target ranking scores.
- Regenerated 198 VHL reference heads by chemical atom maps from deposited4W9H/3JF bound-parent coordinates and fitted them through sequence-verified ELOC. All198 pass bond-geometry checks. They are separate derivatives, not overwrites of the old PT0 files.
- Attempted sequence and CE structural transfer of a chemically correct foretinib head from3LQ8 into all198 shortlisted archived kinase models. All198 failed the existing domain/pocket fit requirements. No threshold was relaxed to invent a successful placement. The same alignment code passes an identical-structure positive control with near-zero pocket RMSD.
- Saved chemically correct native MET/foretinib and VHL parent references as a separate control bundle. These are independently bound parents, not an observed ternary complex.

## Before/after result

The repair assessment requires maximum intercomponent and nonbonded ligand overlap≤0.6Å; the internal-ligand check is stricter than the earlier0.8Å diagnostic limit. It also requires sensible bonded distances1.0–2.1Å, preserved archived graph/3D stereo, optimizer convergence and at least three atoms from each chemical head within4Å. Fixed protein–protein clashes are included separately, so a ligand-only repair does not incorrectly rescue an overlapping protein arrangement.

{len(selected)} historical poses pass these operational geometry/contact checks: {', '.join(selected)}. They are useful only for studying the archived molecular model. They are not validated binders or supplier-matched SJF8240 candidates, and their pocket engagement remains unverified.

| Kinase | Ligand overlap after repair, Å | Internal overlap after repair, Å | Fixed protein overlap, Å | Smaller head contact count at4Å | Geometry/contact result |
|---|---:|---:|---:|---:|---|
'''
 for r in a.itertuples():text+=f'| {r.gene} | {r.after_ligand_external_overlap_A:.3f} | {r.after_internal_overlap_A:.3f} | {r.fixed_E3_POI_overlap_A:.3f} | {r.balanced_head_contacts_4A} | '+('Pass for archived-model diagnostics only' if r.historical_geometry_candidate_4A else 'Fail')+' |\n'
 text+='''
EPHA1 did not meet the optimizer convergence criterion in the stronger repair stage, despite its reduced overlaps. TXK remains strained: one covalent bond is approximately2.107Å, above the unchanged2.1Å bound, and it lacks enough4Å contacts. MAP3K14 retains a protein–protein overlap above0.6Å. RPS6KB1, PDIK1L and PIK3CG retain substantial fixed protein–protein overlaps. MAP3K13 loses correct-head proximity. These failures are structural outcomes, not non-binding or non-degradation labels. Alternative cutoffs4.5Å and5Å are recorded in REPAIR_COMPARISON.tsv; they do not fix molecular provenance or pocket-validation failures.

## Molecular identity — do not relabel existing results

'''
 text+=f"The saved historical ligand graph has InChIKey `{identity['historical_graph_inchikey']}` without stereo; the historical intended SMILES is `{identity['historical_intended_inchikey']}`. The actual archived3D coordinates in these ten poses imply `MCJFROJTDORRBB-OVEQSRFISA-N`. Repairs preserve that actual3D graph/stereo. The supplier-verified SJF8240 is `{identity['supplier_inchikey']}`. These have different chemical connectivity, not merely different conformations. A geometry minimization cannot convert one into the other.\n\n"
 text+='''The historical reference-preparation code assumes that raw PDB coordinate ordering matches an independently parsed reference SMILES ordering, then overwrites part of an independently embedded molecule. That is an unsafe atom-correspondence procedure and provides a concrete mechanism for the head geometry failures. The new reference route uses CCD atom IDs and graph checks. The exact full historical execution chain has not been reconstructed, so the code finding is not claimed as the sole proven cause.

## What remains before an actual SJF8240 top ten

The VHL reference geometry and bounded local repair scripts are fixed and checked. A defensible supplier-matched top ten still requires suitable target-pocket structures/conformations and chemically correct ternary construction. Existing archived kinase models cannot be certified through the tested native-head transfer. That work is a new candidate-generation branch, not a repair of completed MD labels or a simple rerun of the old contact ranking. Do not launch new MD from these historical repaired poses as though they were supplier SJF8240.

## Files

- REPAIR_COMPARISON.tsv: all before/after results and three cutoff checks.
- stronger_relaxation/legacy_repairs/<gene>/: repaired ligand SDF, observed-fragment PDB, complete reconstructed-E3 PDB, atom maps and first-party receipts.
- legacy_repairs/ and probe_archive/: preserved weaker-stage repairs and initial probe.
- supplier_reference_candidates/:198 repaired VHL heads and explicit kinase-placement failure receipts.
- supplier_native_MET_control/: correct chemical/native-parent control inputs; production_ready remains false.
- chemical_identity_audit.json, reference_rebuild_summary.tsv, independent_repair_verification.json: provenance, failures and checks.
- scripts/: complete source for identity auditing, reference rebuilding, local repair, stronger refinement and verification/reporting.

Numerical replay, file integrity, actual coordinate stereochemistry, fixed protein coordinates, all198 rebuilt VHL head geometries and a reference-alignment positive control pass. This is software/geometry validation, not prospective biological validation.
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
 yy=np.arange(len(a));fig,ax=plt.subplots(figsize=(10,5),layout='constrained');ax.barh(yy-.17,a.before_external_overlap_A,height=.3,label='Before: largest ligand or protein overlap',color='#ba6460');ax.barh(yy+.17,np.maximum(a.after_ligand_external_overlap_A,a.fixed_E3_POI_overlap_A),height=.3,label='After: includes fixed protein clashes',color='#437da0');ax.axvline(.6,color='black',ls='--',label='Primary overlap threshold');ax.set_yticks(yy,a.gene);ax.invert_yaxis();ax.set_xlabel('Maximum intercomponent heavy-atom overlap (Å)');ax.set_title('Bounded local repairs of the historical molecular model');ax.legend(loc='lower right');fig.savefig(ROOT/'repair_comparison.png',dpi=180);fig.savefig(ROOT/'repair_comparison.pdf');plt.close(fig)
 sourcefiles=[p for p in (ROOT/'scripts').glob('*.py')];(ROOT/'SOURCE_MANIFEST.json').write_text(json.dumps({str(p.relative_to(ROOT)):c.sha(p) for p in sourcefiles},indent=2)+'\n');print('PASS: 10 independent repaired-geometry replays and198 VHL heads. Historical4A candidates:',', '.join(selected))
if __name__=='__main__':main()
