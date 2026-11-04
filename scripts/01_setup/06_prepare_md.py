#!/usr/bin/env python3
"""Build a mapped AMBER99SB-ILDN/GAFF2/TIP3P system and consistent MD stages.

Based on protac-md's actual staged-release workflow. No -missing or blanket
-maxwarn; ligand placement uses chemical atom maps, not sequential PDB order.
Normal mode prepares the system. Production execution is a separate command.
"""
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,pdb_atoms,atomic_json,sha256
spec=importlib.util.spec_from_file_location('ligand_parameters',Path(__file__).with_name('04_parameterize_ligands.py'))
ligparams=importlib.util.module_from_spec(spec);spec.loader.exec_module(ligparams)

def gmx(command,cwd,log,stdin=None):
    env=os.environ.copy();env['GMX_MAXBACKUP']='-1'
    result=subprocess.run(command,cwd=cwd,input=stdin,text=True,capture_output=True,env=env)
    (Path(cwd)/log).write_text(result.stdout+result.stderr)
    if result.returncode:raise ValidationError('GROMACS failed: '+' '.join(command)+'; see '+str(Path(cwd)/log))
    return result

def write_mdps(directory,seed,production_ns=20):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    nonbonded='''cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
coulombtype = PME
vdwtype = Cut-off
DispCorr = EnerPres
pme-order = 4
fourierspacing = 0.12
pbc = xyz
'''
    em='integrator = steep\nnsteps = 50000\nemtol = 500\nemstep = 0.01\nconstraints = none\n'+nonbonded
    (directory/'em.mdp').write_text(em)
    (directory/'ions.mdp').write_text(em.replace('coulombtype = PME','coulombtype = Cut-off'))
    dynamics='''integrator = md
dt = 0.002
constraints = h-bonds
constraint-algorithm = lincs
lincs-iter = 2
lincs-order = 4
tcoupl = V-rescale
tc-grps = System
tau-t = 1.0
ref-t = 300
comm-mode = Linear
comm-grps = System
nstcomm = 100
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = 5000
compressed-x-grps = SOLUTE
nstenergy = 1000
nstlog = 1000
'''+nonbonded
    nvt=f'define = -DPOSRES -DPOSRES_LIG\nnsteps = 100000\ncontinuation = no\n'+dynamics+f'pcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = {seed}\nld-seed = {seed}\n'
    (directory/'nvt.mdp').write_text(nvt)
    pressure='pcoupl = C-rescale\npcoupltype = isotropic\ntau-p = 5.0\nref-p = 1.0\ncompressibility = 4.5e-5\nrefcoord-scaling = com\n'
    for name,defines in [('npt_restrained','-DPOSRES -DPOSRES_LIG'),('npt','-DPOSRES')]:
        (directory/(name+'.mdp')).write_text(f'define = {defines}\nnsteps = 250000\ncontinuation = yes\n'+dynamics+pressure+f'gen-vel = no\nld-seed = {seed}\n')
    (directory/'prod.mdp').write_text(f'nsteps = {round(production_ns*500000)}\ncontinuation = yes\n'+dynamics+
        pressure.replace('pcoupl = C-rescale','pcoupl = Parrinello-Rahman').replace('refcoord-scaling = com\n','')+f'gen-vel = no\nld-seed = {seed}\n')

def ligand_coordinates(pdb,selector,smiles,mapping):
    atoms=pdb_atoms(pdb)
    chosen=[a for a in atoms if a['record']==selector['record'] and a['chain'] in selector['chains']
            and a['element'] not in {'H','D'} and ('resname' not in selector or a['resname']==selector['resname'])]
    pose=Chem.MolFromPDBBlock('\n'.join(a['line'] for a in chosen)+'\nEND\n',sanitize=False,removeHs=True)
    full=Chem.MolFromSmiles(smiles)
    pose=AllChem.AssignBondOrdersFromTemplate(full,pose)
    Chem.AssignStereochemistryFrom3D(pose,replaceExistingTags=True)
    if Chem.MolToInchiKey(pose)!=Chem.MolToInchiKey(full):raise ValidationError('Starting ligand differs from parameterized stereochemistry')
    matches=pose.GetSubstructMatches(full,useChirality=True)
    if not matches:raise ValidationError('No stereochemical atom map from pose to ligand topology')
    match=matches[0];placed=Chem.Mol(full);placed.RemoveAllConformers();conf=Chem.Conformer(full.GetNumAtoms())
    for i,j in enumerate(match):conf.SetAtomPosition(i,pose.GetConformer().GetAtomPosition(j))
    placed.AddConformer(conf)
    hyd=Chem.AddHs(placed,addCoords=True)
    ff=AllChem.UFFGetMoleculeForceField(hyd)
    for i in range(full.GetNumAtoms()):ff.AddFixedPoint(i)
    if ff.Minimize(maxIts=2000):raise ValidationError('Mapped ligand hydrogen minimization did not converge')
    if hyd.GetNumAtoms()!=len(mapping['atom_index_to_name']):raise ValidationError('Hydrogen/topology atom-count mismatch')
    positions=hyd.GetConformer().GetPositions()
    if np.max(np.linalg.norm(positions[:full.GetNumAtoms()]-placed.GetConformer().GetPositions(),axis=1))>1e-6:
        raise ValidationError('Heavy atoms moved during hydrogen placement')
    coords={name:positions[int(index)]/10 for index,name in mapping['atom_index_to_name'].items()}
    return coords,{'pose_atom_mapping':list(match),'symmetry_match_count':len(matches),
                   'heavy_atoms_fixed':True,'hydrogens':'reconstructed and minimized in the bound frame'}

def build(root,compound,pose,protein,selection,work,gmx_path,seed,production_ns,parameter_state='neutral_reference'):
    work=Path(work).resolve()
    if work.exists() and any(work.iterdir()):raise ValidationError('MD preparation directory already used; preserve it and use a fresh directory')
    work.mkdir(parents=True,exist_ok=True)
    if Path(parameter_state).name!=parameter_state or parameter_state in {'.','..'}:
        raise ValidationError('Unsafe parameter-state identifier')
    params=root/'ligand_params'/compound/parameter_state
    report=json.loads((params/'parameter_validation.json').read_text())
    if report['status']!='PARAMETERS_VERIFIED':raise ValidationError('Ligand parameters are not verified')
    for filename,digest in report['files'].items():
        if sha256(params/filename)!=digest:raise ValidationError('Ligand parameter artifact changed: '+filename)
    repair=json.loads(protein.with_suffix('.repair.json').read_text())
    if repair['status']!='PROTEIN_MODEL_COMPLETED' or sha256(protein)!=repair['output_sha256']:
        raise ValidationError('Complete, verified protein model required before MD')
    mapping=json.loads((params/'atom_mapping.json').read_text())
    if report.get('atom_mapping_sha256')!=sha256(params/'atom_mapping.json'):
        raise ValidationError('Ligand atom mapping is not bound to its parameter validation receipt')
    coordinates,ligand_map=ligand_coordinates(pose,selection['protac'],mapping['smiles'],mapping)
    shutil.copy2(protein,work/'protein_input.pdb')
    ff_source=gmx_path.resolve().parent.parent/'share/gromacs/top/amber99sb-ildn.ff'
    if not (ff_source/'forcefield.itp').exists():
        raise ValidationError('Cannot locate the explicit GROMACS AMBER99SB-ILDN data directory')
    shutil.copytree(ff_source,work/'amber99sb-ildn.ff')
    gmx([str(gmx_path),'pdb2gmx','-f','protein_input.pdb','-o','protein_processed.pdb','-p','topol.top',
         '-ff','amber99sb-ildn','-water','tip3p','-ignh','-merge','no'],work,'pdb2gmx.log')
    gmx([str(gmx_path),'editconf','-f','protein_processed.pdb','-o','protein.gro'],work,'protein_gro.log')
    protein_atoms=pdb_atoms(work/'protein_processed.pdb')
    natoms=len(protein_atoms)
    gro=(work/'protein.gro').read_text().splitlines()
    if int(gro[1])!=natoms:raise ValidationError('Processed protein PDB/GRO atom-order contract violated')
    itp=params/report['topology'];sections=ligparams.itp_sections(itp)
    atoms=sections['atoms'];name=sections['moleculetype'][0][0]
    ligand_itp=itp.read_text()+'\n#ifdef POSRES_LIG\n#include "posre_lig.itp"\n#endif\n'
    (work/'ligand.itp').write_text(ligand_itp)
    (work/'posre_lig.itp').write_text('[ position_restraints ]\n; ai funct fcx fcy fcz\n'+
        ''.join(f'{int(a[0])} 1 1000 1000 1000\n' for a in atoms if mapping['elements'][a[4]]!='H'))
    top=(work/'topol.top').read_text();ff='#include "amber99sb-ildn.ff/forcefield.itp"'
    local_ff='#include "./amber99sb-ildn.ff/forcefield.itp"'
    if top.count(ff)+top.count(local_ff)!=1:raise ValidationError('Cannot identify one local protein force-field include')
    top=top.replace(local_ff,ff)
    if top.count(ff)!=1:raise ValidationError('Cannot place GAFF2 atomtypes before all molecule definitions')
    top=top.replace(ff,ff+'\n#include "ligand.itp"',1).rstrip()+f'\n{name} 1\n'
    (work/'topol.top').write_text(top)
    ligand_lines=[]
    for i,a in enumerate(atoms,1):
        xyz=coordinates[a[4]]
        ligand_lines.append(f'{1:5d}{"PTC":<5s}{a[4]:>5s}{(natoms+i)%100000:5d}'+''.join(f'{v:8.3f}' for v in xyz))
    merged=['Mapped protein + '+compound,str(natoms+len(atoms)),*gro[2:2+natoms],*ligand_lines,gro[-1]]
    (work/'merged.gro').write_text('\n'.join(merged)+'\n')
    indices={role:[i+1 for i,a in enumerate(protein_atoms) if a['chain'] in selection[role]['chains'] and a['element'] not in {'H','D'}]
             for role in ['vhl','poi']}
    indices['protac']=[natoms+i+1 for i,a in enumerate(atoms) if mapping['elements'][a[4]]!='H']
    if not all(indices.values()):raise ValidationError('An MD atom-selection group is empty')
    indices['solute']=list(range(1,natoms+len(atoms)+1))
    with (work/'index.ndx').open('w') as f:
        for role,group in indices.items():
            f.write(f'[ {role.upper()} ]\n')
            for i in range(0,len(group),15):f.write(' '.join(map(str,group[i:i+15]))+'\n')
    write_mdps(work/'mdp',seed,production_ns)
    gmx([str(gmx_path),'editconf','-f','merged.gro','-o','boxed.gro','-bt','dodecahedron','-d','1.5'],work,'box.log')
    gmx([str(gmx_path),'solvate','-cp','boxed.gro','-cs','spc216.gro','-o','solvated.gro','-p','topol.top'],work,'solvate.log')
    gmx([str(gmx_path),'grompp','-f','mdp/ions.mdp','-c','solvated.gro','-p','topol.top','-o','ions.tpr'],work,'grompp_ions.log')
    # Group names, not brittle numeric group IDs. Solvent replacement preserves solute indices.
    gmx([str(gmx_path),'genion','-s','ions.tpr','-o','system.gro','-p','topol.top','-pname','NA','-nname','CL',
         '-neutral','-conc','0.15','-seed',str(seed)],work,'genion.log','SOL\n')
    system_count=int((work/'system.gro').read_text().splitlines()[1])
    with (work/'index.ndx').open('a') as f:
        f.write('[ System ]\n')
        for start in range(1,system_count+1,15):
            f.write(' '.join(map(str,range(start,min(start+15,system_count+1))))+'\n')
    gmx([str(gmx_path),'grompp','-f','mdp/em.mdp','-c','system.gro','-p','topol.top','-n','index.ndx','-o','em.tpr'],work,'grompp_em.log')
    result={'status':'TOPOLOGY_PREPARED','compound':compound,'pose_sha256':sha256(pose),'protein_repair':repair,
        'parameter_state':parameter_state,'ligand_microstate':report['microstate'],
        'ligand_formal_charge':report['formal_charge'],
        'ligand_parameter_validation_sha256':sha256(params/'parameter_validation.json'),
        'ligand_atom_mapping_sha256':sha256(params/'atom_mapping.json'),
        'seed':seed,'temperature_K':300,'salt_molar':0.15,'box_padding_nm':1.5,
        'force_fields':['AMBER99SB-ILDN','GAFF2','TIP3P'],'ligand_mapping':ligand_map,
        'atom_indices_one_based':indices,'solute_atoms':natoms+len(atoms),
        'production_ns':production_ns,'complete_equilibration_required':True,
        'files':{str(f.relative_to(work)):sha256(f) for f in work.rglob('*') if f.is_file() and not f.name.endswith('.log')}}
    atomic_json(work/'md_preparation.json',result)
    print(json.dumps({k:result[k] for k in ['status','compound','solute_atoms','seed']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--compound',required=True)
    p.add_argument('--pose',type=Path,required=True);p.add_argument('--protein',type=Path,required=True)
    p.add_argument('--selection',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--gmx',type=Path,required=True);p.add_argument('--seed',type=int,required=True)
    p.add_argument('--production-ns',type=float,default=20)
    p.add_argument('--parameter-state',default='neutral_reference')
    a=p.parse_args();build(a.root,a.compound,a.pose,a.protein,json.loads(a.selection.read_text()),a.work,a.gmx,a.seed,a.production_ns,a.parameter_state)
