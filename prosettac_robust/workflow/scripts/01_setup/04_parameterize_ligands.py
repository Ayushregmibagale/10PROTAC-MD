#!/usr/bin/env python3
"""Generate separate AM1-BCC/GAFF2 parameter sets with explicit atom-name maps.

    No charge-method fallback, no template borrowed from another PROTAC. Uses the
    exact supplier neutral-reference graph; solution microstates are a distinct
    scientific design choice and are recorded rather than implicitly changed.
"""
import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from rdkit import Chem
import rdkit
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,atomic_json,sha256

def atom_type(atom):
    symbol=atom.GetSymbol()
    if symbol=='C':return 'C.ar' if atom.GetIsAromatic() else 'C.2' if atom.GetHybridization()==Chem.HybridizationType.SP2 else 'C.3'
    if symbol=='N':
        if atom.GetIsAromatic():return 'N.ar'
        if any(n.GetSymbol()=='C' and any(b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(n).GetSymbol()=='O'
                    for b in n.GetBonds()) for n in atom.GetNeighbors()):return 'N.am'
        return 'N.4' if atom.GetFormalCharge()>0 else 'N.2' if atom.GetHybridization()==Chem.HybridizationType.SP2 else 'N.3'
    if symbol=='O':return 'O.2' if any(b.GetBondType()==Chem.BondType.DOUBLE for b in atom.GetBonds()) else 'O.3'
    if symbol=='S':return 'S.2' if atom.GetIsAromatic() else 'S.3'
    if symbol in {'H','F','Cl','Br','I'}:return symbol
    raise ValidationError('Unsupported ligand element: '+symbol)

def prepare_input(root,compound,destination,seed,state_id=None,assessment=None):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    supplier=json.loads((root/'config/molecules.json').read_text())['compounds'][compound]
    smiles=(root/'molecules'/compound/'full.smi').read_text().strip()
    microstate='supplier_neutral_reference';assessment_sha256=None
    if state_id:
        if assessment is None:raise ValidationError('Explicit assessment file required for a predicted microstate')
        report=json.loads(Path(assessment).read_text())
        state=report['compounds'][compound]['states'][state_id]
        smiles=state['smiles'];microstate='predicted_'+state_id;assessment_sha256=sha256(assessment)
    full=Chem.MolFromSmiles(smiles)
    expected=Chem.MolToInchiKey(full)
    parent=rdMolStandardize.Uncharger().uncharge(Chem.Mol(full)) if state_id else full
    if Chem.MolToInchiKey(parent)!=supplier['inchikey']:raise ValidationError('Supplier parent graph/stereochemistry mismatch')
    embedding_attempts=[]
    for trial in range(8):
        mol=Chem.AddHs(full);p=AllChem.ETKDGv3();p.randomSeed=seed+173*trial;p.numThreads=1;p.enforceChirality=True
        p.useRandomCoords=trial>0
        accepted=AllChem.EmbedMolecule(mol,p)>=0
        embedding_attempts.append({'seed':p.randomSeed,'random_coordinates':p.useRandomCoords,'accepted':accepted})
        if accepted:break
    else:raise ValidationError('Parameterization conformer embedding failed for all recorded seeds')
    if not AllChem.UFFHasAllMoleculeParams(mol):raise ValidationError('Missing preoptimization parameters')
    if AllChem.UFFOptimizeMolecule(mol,maxIters=5000)!=0:raise ValidationError('Parameter conformer did not converge')
    check=Chem.RemoveHs(Chem.Mol(mol));Chem.AssignStereochemistryFrom3D(check,replaceExistingTags=True)
    if Chem.MolToInchiKey(check)!=expected:raise ValidationError('Parameterization stereochemistry changed')
    names={i:a.GetSymbol()+str(i+1) for i,a in enumerate(mol.GetAtoms())}
    if any(len(n)>4 for n in names.values()):raise ValidationError('Atom name exceeds PDB capacity')
    lines=['@<TRIPOS>MOLECULE',compound,f'{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0','SMALL','USER_CHARGES','','@<TRIPOS>ATOM']
    for i,atom in enumerate(mol.GetAtoms()):
        x,y,z=mol.GetConformer().GetAtomPosition(i)
        lines.append(f'{i+1:7d} {names[i]:<4} {x:12.6f} {y:12.6f} {z:12.6f} {atom_type(atom):<5} 1 PTC {0.0:10.6f}')
    lines.append('@<TRIPOS>BOND')
    for i,bond in enumerate(mol.GetBonds(),1):
        order='ar' if bond.GetIsAromatic() else str(int(bond.GetBondTypeAsDouble()))
        lines.append(f'{i:6d} {bond.GetBeginAtomIdx()+1:5d} {bond.GetEndAtomIdx()+1:5d} {order}')
    lines+=['@<TRIPOS>SUBSTRUCTURE','1 PTC 1 TEMP 0 **** **** 0 ROOT']
    source=destination/'input.mol2';source.write_text('\n'.join(lines)+'\n')
    w=Chem.SDWriter(str(destination/'input.sdf'));w.write(mol);w.close()
    meta={'supplier_inchikey':supplier['inchikey'],'parameter_molecule_inchikey':expected,
        'preparation_versions':{'python':sys.version,'rdkit':rdkit.__version__},
        'smiles':smiles,'seed':seed,'embedding_attempts':embedding_attempts,'protonation_assessment_sha256':assessment_sha256,
        'formal_charge':sum(a.GetFormalCharge() for a in full.GetAtoms()),
        'charge_method':'AM1-BCC','force_field':'GAFF2','microstate':microstate,
        'atom_index_to_name':names,'heavy_atom_count':full.GetNumAtoms(),
        'elements':{names[i]:a.GetSymbol() for i,a in enumerate(mol.GetAtoms())},
        'bond_names':[[names[b.GetBeginAtomIdx()],names[b.GetEndAtomIdx()]] for b in mol.GetBonds()],
        'input_sha256':sha256(source),'status':'INPUT_PREPARED'}
    atomic_json(destination/'atom_mapping.json',meta)
    return source,meta

def itp_sections(path):
    sections={};current=None
    for line in Path(path).read_text().splitlines():
        line=line.split(';',1)[0].strip()
        if not line or line.startswith('#'):continue
        if line.startswith('['):current=line.strip('[] ').lower();sections.setdefault(current,[]);continue
        if current:sections[current].append(line.split())
    return sections

def validate_parameters(destination,compound,meta):
    files=list(Path(destination).glob('*.acpype/'+compound+'_GMX.itp'))
    if len(files)!=1:raise ValidationError('Missing/ambiguous generated ligand topology')
    itp=files[0];sections=itp_sections(itp);atoms=sections.get('atoms',[])
    expected=set(meta['atom_index_to_name'].values())
    found=[a[4] for a in atoms]
    if set(found)!=expected:
        # AmberTools uppercases two-character element atom names (Cl -> CL).
        # Permit only this explicit, bijective alias; connectivity is checked below.
        aliases={name:name.upper() for name in expected-set(found)
                 if meta['elements'][name] in {'Cl','Br'} and name.upper() in set(found)-expected}
        if len(found)==len(expected) and {aliases.get(n,n) for n in expected}==set(found):
            original=json.loads(json.dumps(meta))
            meta=dict(meta,tool_atom_name_aliases=aliases,
                original_input_atom_index_to_name=meta['atom_index_to_name'],
                atom_index_to_name={k:aliases.get(v,v) for k,v in meta['atom_index_to_name'].items()},
                elements={aliases.get(k,k):v for k,v in meta['elements'].items()},
                bond_names=[[aliases.get(n,n) for n in bond] for bond in meta['bond_names']])
            atomic_json(Path(destination)/'atom_mapping_before_name_normalization.json',original)
            atomic_json(Path(destination)/'atom_mapping.json',meta)
            expected=set(meta['atom_index_to_name'].values())
    if len(found)!=len(expected) or set(found)!=expected:raise ValidationError('Parameter atom names/count differ from mapped molecule')
    charges=[float(a[6]) for a in atoms];masses=[float(a[7]) for a in atoms]
    if not all(math.isfinite(c) for c in charges+masses) or min(masses)<=0:
        raise ValidationError('Invalid charge/mass in generated topology')
    if abs(sum(charges)-meta['formal_charge'])>0.001:raise ValidationError('Partial charges do not sum to formal molecular charge')
    numbered={int(a[0]):a[4] for a in atoms}
    bonds={frozenset([numbered[int(b[0])],numbered[int(b[1])]]) for b in sections.get('bonds',[])}
    if bonds!={frozenset(b) for b in meta['bond_names']}:raise ValidationError('Topology connectivity differs from source atom mapping')
    lines=(itp.parent/(compound+'_bcc_gaff2.mol2')).read_text().splitlines()
    start=lines.index('@<TRIPOS>ATOM')+1
    stop=lines.index('@<TRIPOS>BOND')
    aliases=meta.get('tool_atom_name_aliases',{})
    coords={aliases.get(x[1],x[1]):list(map(float,x[2:5])) for line in lines[start:stop] if (x:=line.split())}
    if set(coords)!=expected:raise ValidationError('Final parameter coordinates do not match mapped atoms')
    output=Chem.AddHs(Chem.MolFromSmiles(meta['smiles']));conf=Chem.Conformer(output.GetNumAtoms())
    for index,name in meta['atom_index_to_name'].items():conf.SetAtomPosition(int(index),coords[name])
    output.AddConformer(conf);output=Chem.RemoveHs(output)
    Chem.AssignStereochemistryFrom3D(output,replaceExistingTags=True)
    if Chem.MolToInchiKey(output)!=meta.get('parameter_molecule_inchikey',meta['supplier_inchikey']):
        raise ValidationError('Quantum-optimized output stereochemistry differs from supplier')
    for file in itp.parent.glob('*.frcmod'):
        if 'ATTN' in file.read_text() or 'need revision' in file.read_text().lower():
            raise ValidationError('Unresolved GAFF2 parameter requires inspection: '+str(file))
    report={'status':'PARAMETERS_VERIFIED','compound':compound,'force_field':'GAFF2','charges':'AM1-BCC',
            'formal_charge':meta['formal_charge'],'partial_charge_sum':sum(charges),'atoms':len(atoms),
            'microstate':meta['microstate'],'files':{str(p.relative_to(destination)):sha256(p)
                    for p in itp.parent.iterdir() if p.is_file()},'topology':str(itp.relative_to(destination))}
    report['output_stereochemistry_verified']=True
    report['atom_mapping_sha256']=sha256(Path(destination)/'atom_mapping.json')
    report['parameter_molecule_inchikey']=meta.get('parameter_molecule_inchikey',meta['supplier_inchikey'])
    atomic_json(Path(destination)/'parameter_validation.json',report)
    return report

def run(root,compound,amber_prefix,timeout,seed,state_id=None,assessment=None):
    if state_id and (Path(state_id).name!=state_id or not state_id.replace('_','').isalnum()):
        raise ValidationError('Unsafe microstate identifier')
    dest=root/'ligand_params'/compound/('microstate_'+state_id if state_id else 'neutral_reference')
    if (dest/'execution.json').exists():raise ValidationError('Parameter attempt already used; choose a new output root')
    source,meta=prepare_input(root,compound,dest,seed,state_id,assessment)
    executable=amber_prefix/'bin/acpype'
    env=os.environ.copy();env['PATH']=str(amber_prefix/'bin')+os.pathsep+env.get('PATH','')
    env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1'
    cmd=[str(executable),'-i',str(source.resolve()),'-b',compound,'-a','gaff2','-c','bcc',
         '-n',str(meta['formal_charge']),'-m','1','-o','gmx','-s',str(timeout)]
    atomic_json(dest/'execution.json',{'status':'RUNNING','command':cmd,'acpype_sha256':sha256(executable)})
    with (dest/'parameterization.log').open('w') as f:
        result=subprocess.run(cmd,cwd=dest,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=timeout+60)
    if result.returncode:
        atomic_json(dest/'execution.json',{'status':'FAILED','command':cmd,'exit_code':result.returncode})
        raise ValidationError('AM1-BCC/GAFF2 generation failed; inspect '+str(dest/'parameterization.log'))
    try:report=validate_parameters(dest,compound,meta)
    except Exception as exc:
        atomic_json(dest/'execution.json',{'status':'GENERATED_VALIDATION_FAILED','command':cmd,'exit_code':0,'reason':str(exc)})
        raise
    atomic_json(dest/'execution.json',{'status':'COMPLETE','command':cmd,'exit_code':0})
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--compound',required=True)
    p.add_argument('--amber-prefix',type=Path,required=True)
    p.add_argument('--timeout-seconds',type=int,default=3600)
    p.add_argument('--seed',type=int,default=104729)
    p.add_argument('--state-id');p.add_argument('--assessment',type=Path)
    a=p.parse_args();run(a.root,a.compound,a.amber_prefix,a.timeout_seconds,a.seed,a.state_id,a.assessment)
