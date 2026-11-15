"""Check static Rosetta head typing and preservation of native AKT phospho-Thr.

This loads inputs only: no docking, optimization, MD or degradation prediction.
Run in the project's licensed Python 3.10 docking runtime.
"""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from rdkit import Chem

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine',type=Path,required=True);p.add_argument('--head-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(a.engine));from official_parameters import generate
    report={'status':'RUNNING','scope':'static parameter/load compatibility only','MD':False,'started_unix':time.time(),
       'source_script_sha256':sha(Path(__file__)),
       'source_hashes':{f:sha(a.head_root/f) for f in ('vhl_head.sdf','poi_modeled_head.sdf','AKT1_parent_protein.pdb')}}
    try:
        params=[]
        for role,name in enumerate(('vhl_head.sdf','poi_modeled_head.sdf')):
            pdb,param=generate(str(a.head_root/name),'PT'+str(role),str(a.output/('PT'+str(role))))
            params.append(param)
        protein=a.head_root/'AKT1_parent_protein.pdb'
        lines=[x for x in protein.read_text().splitlines() if x.startswith(('ATOM  ','HETATM'))]
        if not any(x[17:20]=='TPO' and int(x[22:26])==308 for x in lines):raise ValueError('Source lacks expected phospho-Thr308')
        combo=a.output/'AKT1_PT1_input.pdb'
        combo.write_text('\n'.join(lines)+'\n'+(a.output/'PT1_0001.pdb').read_text()+'\nEND\n')
        import pyrosetta
        # Output paths are deliberately space-free for the Rosetta option parser.
        if any(' ' in x for x in params):raise ValueError('Use a space-free output path for Rosetta extra_res_fa')
        flags='-mute all -constant_seed -jran 42 -ignore_unrecognized_res false -load_PDB_components false -extra_res_fa '+params[0]+' -extra_res_fa '+params[1]
        report['rosetta_flags']=flags
        pyrosetta.init(flags)
        pose=pyrosetta.pose_from_file(str(combo));info=pose.pdb_info()
        phospho=[i for i in range(1,pose.total_residue()+1) if info.chain(i)=='D' and info.number(i)==308]
        if len(phospho)!=1 or not pose.residue(phospho[0]).has('P'):raise ValueError('Rosetta lost or changed phosphate at residue308')
        ligands=[i for i in range(1,pose.total_residue()+1) if pose.residue(i).name3()=='PT1']
        if len(ligands)!=1:raise ValueError('AKT head missing/ambiguous after load')
        head=Chem.SDMolSupplier(str(a.head_root/'poi_modeled_head.sdf'),removeHs=True)[0]
        provenance=json.loads((a.output/'PT1_parameter_provenance.json').read_text());mapped=provenance['all_atom_index_to_name']
        coordinates=np.array([list(pose.residue(ligands[0]).xyz(mapped[str(i)])) for i in range(head.GetNumAtoms())])
        displacement=float(np.linalg.norm(coordinates-head.GetConformer().GetPositions(),axis=1).max())
        if displacement>.003:raise ValueError('Head coordinates changed during parameter/load round trip')
        input_ca={(x[21],int(x[22:26])):np.array([float(x[30:38]),float(x[38:46]),float(x[46:54])])
                  for x in lines if x[12:16].strip()=='CA'}
        errors=[]
        for i in range(1,pose.total_residue()+1):
            key=(info.chain(i),info.number(i))
            if key in input_ca:
                if not pose.residue(i).has('CA'):raise ValueError('Source backbone atom lost')
                errors.append(float(np.linalg.norm(np.array(list(pose.residue(i).xyz('CA')))-input_ca[key])))
        if len(errors)!=len(input_ca) or max(errors,default=0)>.003:raise ValueError('Source protein CA coordinates/count changed')
        pose.dump_pdb(str(a.output/'AKT1_PT1_loaded.pdb'))
        report.update(status='STATIC_INPUT_LOAD_PASSED_NOT_SCIENTIFIC_RELEASE',protein_CA_count=len(errors),
          max_protein_CA_displacement_A=max(errors),head_heavy_atoms=head.GetNumAtoms(),max_head_displacement_A=displacement,
          phospho_residue_type=pose.residue(phospho[0]).name(),phosphate_atom_preserved=True,
          parameter_hashes={Path(x).name:sha(Path(x)) for x in params},
          limitations=['Input loading is not a binding-mode validation','No kinase-domain transfer or ternary docking performed',
                      'Protonation/assay-state review remains pending; no dominant solution-state population claimed'])
    except BaseException as exc:
        report.update(status='FAILED',error=repr(exc));raise
    finally:
        report['finished_unix']=time.time();(a.output/'ROSETTA_INPUT_PROBE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
