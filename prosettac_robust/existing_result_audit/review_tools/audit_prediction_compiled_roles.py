"""Sequence-identify prediction-target compiled proteins without legacy indices."""
from pathlib import Path
import json,sys,importlib.util
import MDAnalysis as mda
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,atomic_json,ValidationError
spec=importlib.util.spec_from_file_location('role_identity',ROOT/'existing_result_audit/review_tools/audit_existing_protein_roles.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
base=ROOT/'existing_result_audit/scientific_signoff_v1';roster=base/'prediction_pose_roster.json';seq=base/'all19_sequences.json';mapping=base/'protein_groups_19.json';groups=json.loads(mapping.read_text());entries={r['primaryAccession']:r['sequence']['value'] for r in json.loads(seq.read_text())['results']};out=ROOT/'existing_result_audit/prediction_compiled_roles_v3'
if out.exists():raise ValidationError('Fresh output required')
rows=[]
for case in json.loads(roster.read_text()):
 p=Path('/home/ayushregmibagale/protac-md/random_screen')/case['poi']/case['pose_id']/'prod/prod.tpr';h=sha256(p);u=mda.Universe(str(p));found=[]
 for kind in sorted(set(map(str,u.atoms.moltypes))):
  if not kind.startswith('Protein'):continue
  ca=u.select_atoms('moltype '+kind+' and name CA');sequence=''.join(seq1(a.resname,undef_code='X') for a in ca)
  found.append(dict(mod.identify(sequence,{k:entries[k] for k in ['P40337','Q15370','Q15369',groups[case['poi']]]}),moltype=kind,CA_atoms=len(ca)))
 if sha256(p)!=h:raise ValidationError('Source drift')
 if {r['uniprot'] for r in found}!={'P40337','Q15370','Q15369',groups[case['poi']]}:raise ValidationError('Missing/ambiguous protein role')
 rows.append(dict(case,proteins=found,tpr=str(p),tpr_sha256=h))
out.mkdir(parents=True);atomic_json(out/'per_pose_roles.json',rows)
report={'status':'COMPILED_PROTEIN_IDENTIFICATION_COMPLETE_LEGACY_INDICES_NOT_AVAILABLE','poses_identified':len(rows),'targets':sorted({r['poi'] for r in rows}),'all_expected_proteins_present':True,'legacy_index_roles_verified':False,'new_MD_executed':False,'original_files_modified':False,'input_sha256':{str(p):sha256(p) for p in [roster,seq,mapping,ROOT/'existing_result_audit/review_tools/audit_existing_protein_roles.py']},'source_sha256':sha256(Path(__file__)),'records_sha256':sha256(out/'per_pose_roles.json')};atomic_json(out/'summary.json',report);print(json.dumps(report,indent=2))
