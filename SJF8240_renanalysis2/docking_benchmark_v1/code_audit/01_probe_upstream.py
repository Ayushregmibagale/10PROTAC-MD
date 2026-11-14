"""Read-only upstream regression probes, isolated outputs and immutable source hashes."""
from pathlib import Path
import json,subprocess,sys,shutil,hashlib,importlib
ROOT=Path(__file__).resolve().parents[1];UP=ROOT/'upstream_PROTAC_ternary';OUT=ROOT/'code_audit/probes';OUT.mkdir(exist_ok=True);sys.path.insert(0,str(UP))
from PROTAC_ternary import PDBContainer,read_atom_file
from rdkit import Chem

def main():
 example=UP/'Example1_ternary_model_prediction';results=[]
 for name,extra in [('example_default',[]),('negative_cutoff',['-c','-1']),('nan_cutoff',['-c','nan']),('zero_alignment_iterations',['-ai','0'])]:
  out=OUT/name;out.mkdir(exist_ok=True);cmd=[sys.executable,str(UP/'ternary_model_prediction.py'),'-d',str(example/'docking_decoy.pdb'),'-l',str(example/'linker_conformer.pdb'),'-da',str(example/'decoy_atom_list.txt'),'-la',str(example/'linker_atom_list.txt'),'-wd',str(example/'decoy_atom_list_delete.txt'),'-ld',str(example/'linker_atom_list_delete.txt'),'-t','default',*extra];run=subprocess.run(cmd,cwd=out,text=True,capture_output=True,timeout=60);(out/'stdout.txt').write_text(run.stdout);(out/'stderr.txt').write_text(run.stderr);results.append({'probe':name,'returncode':run.returncode,'ternary_outputs':len(list(out.glob('ternary*.pdb'))),'rmsd_text':(out/'rmsd.txt').read_text() if (out/'rmsd.txt').exists() else None,'command':cmd})
 # Missing deletion labels are ignored by upstream rather than validated.
 model=PDBContainer(str(example/'docking_decoy.pdb'),id_format=1);n=model.structure.GetNumAtoms();model.delete_atoms(['DOES_NOT_EXIST C1']);results.append({'probe':'nonexistent_delete_atom','silently_ignored':model.structure.GetNumAtoms()==n})
 malformed=OUT/'malformed_atoms.txt';malformed.write_text('XXXX C1\n')
 try:read_atom_file(str(malformed));value='accepted'
 except Exception as error:value=type(error).__name__+':'+str(error)
 results.append({'probe':'malformed_atom_format','outcome':value})
 # Duplicate unqualified labels overwrite the dictionary entry; demonstrate with two isolated repeated records.
 pdb=OUT/'duplicate_label.pdb';pdb.write_text('HETATM    1  C1  UNL A   1       0.000   0.000   0.000  1.00  0.00           C  \nHETATM    2  C1  UNL B   2      20.000   0.000   0.000  1.00  0.00           C  \nEND\n');m=PDBContainer(str(pdb),id_format=1);results.append({'probe':'duplicate_unqualified_atom_label','atoms':m.structure.GetNumAtoms(),'dictionary_entries':len(m.atom_ids),'selected_index':m.atom_ids.get('UNL C1')})
 (OUT/'RESULTS.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
if __name__=='__main__':main()
