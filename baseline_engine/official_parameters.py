"""Official Rosetta ligand parameters within the existing adapter's API.

Preserves every input conformer's coordinate frame and explicit heavy-index map.
This repairs chemical parameter generation; it does not sample docking poses.
"""
import hashlib,json,subprocess,sys
import os
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from rdkit import Chem

def generate(sdf_file,name,pdb_prefix,overwrite=True,conformers=False,nbr=-1,v_atoms_sdf=''):
    prefix=Path(pdb_prefix).resolve();params=prefix.with_suffix('.params')
    if params.exists() and not overwrite:raise RuntimeError('Parameter artifact already exists; use a fresh attempt')
    molecules=list(Chem.SDMolSupplier(str(sdf_file),removeHs=True))
    if not molecules or any(m is None or m.GetNumConformers()!=1 for m in molecules):raise RuntimeError('Invalid conformer input')
    expected=Chem.MolToInchiKey(molecules[0])
    if any(Chem.MolToInchiKey(m)!=expected for m in molecules):raise RuntimeError('Different graphs/stereochemistry in conformer ensemble')
    hyd=[Chem.AddHs(m,addCoords=True) for m in molecules]
    frame=None
    if conformers and len(hyd)>1:
        if not 0<=nbr<molecules[0].GetNumAtoms():raise RuntimeError('Multiple ligand conformers require an explicit fixed orientation center')
        frame=[nbr]+[n.GetIdx() for n in molecules[0].GetAtomWithIdx(nbr).GetNeighbors() if n.GetAtomicNum()>1]
        if len(frame)<3:raise RuntimeError('Rotamer orientation center needs two fixed heavy neighbors')
        reference=hyd[0].GetConformer().GetPositions()[frame]
        displacement=max(float(np.max(np.linalg.norm(h.GetConformer().GetPositions()[frame]-reference,axis=1))) for h in hyd)
        if displacement>0.003:raise RuntimeError('Ligand rotamer orientation frame is not invariant across conformers')
    for h in hyd:
        check=Chem.RemoveHs(Chem.Mol(h));Chem.AssignStereochemistryFrom3D(check,replaceExistingTags=True)
        if Chem.MolToInchiKey(check)!=expected:raise RuntimeError('Conformer 3D stereochemistry differs from input graph')
    source=prefix.with_name(prefix.name+'_explicit_H.sdf')
    writer=Chem.SDWriter(str(source));writer.SetKekulize(False)
    for h in hyd:writer.write(h)
    writer.close()
    script=Path(__file__).resolve().parent/'vendor/molfile_to_params.py'
    command=[sys.executable,str(script),str(source),'-n',name,'-p',str(prefix),'--chain','X']
    if overwrite:command.append('--clobber')
    if conformers:command.append('--conformers-in-one-file')
    if nbr!=-1:
        if not 0<=nbr<molecules[0].GetNumAtoms():raise RuntimeError('Neighbour index outside full ligand')
        command+=['--nbr_atom',str(nbr+1)]
    result=subprocess.run(command,text=True,capture_output=True)
    prefix.with_name(prefix.name+'_parameterization.log').write_text(result.stdout+result.stderr)
    if result.returncode:raise RuntimeError('Official Rosetta parameter generation failed; see log')
    pdb=Path(str(prefix)+('.pdb' if conformers else '_0001.pdb'))
    records=[l for l in pdb.read_text().splitlines() if l.startswith('HETATM')]
    if len(records)!=hyd[0].GetNumAtoms():raise RuntimeError('Parameter atom count mismatch')
    xyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in records])
    target=hyd[0].GetConformer().GetPositions();offset=xyz.mean(0)-target.mean(0)
    dist=np.linalg.norm((target+offset)[:,None,:]-xyz[None,:,:],axis=2)
    for i,a in enumerate(hyd[0].GetAtoms()):
        for j,line in enumerate(records):
            if a.GetSymbol().upper()!=line[76:78].strip().upper():dist[i,j]=1e9
    rows,cols=linear_sum_assignment(dist)
    if dist[rows,cols].max()>0.003:raise RuntimeError('Cannot map generated parameter names to all input atoms')
    mapping={int(i):records[int(j)][12:16].strip() for i,j in zip(rows,cols)}
    pdb_fields={int(i):records[int(j)][12:16] for i,j in zip(rows,cols)}
    def block(mol):
        lines=[]
        for i,a in enumerate(mol.GetAtoms()):
            x,y,z=mol.GetConformer().GetAtomPosition(i)
            lines.append(f'HETATM{i+1:5d} {pdb_fields[i]} {name:3s} X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a.GetSymbol():>2s}  ')
        return '\n'.join(lines)+'\nTER\n'
    # Official generation centers additional conformers; restore each exact input
    # frame. The explicit invariant NBR frame is also required during Rosetta's
    # conformer orientation; restoring Cartesian files alone is insufficient.
    pdb.write_text(block(hyd[0]))
    if conformers:
        rotamers=Path(str(prefix)+'_conformers.pdb')
        rotamers.write_text(''.join(block(m) for m in hyd[1:]))
        text=params.read_text();text='\n'.join(l for l in text.splitlines() if not l.startswith('PDB_ROTAMERS'))+'\n'
        if len(hyd)>1:text+='PDB_ROTAMERS '+str(rotamers)+'\n'
        params.write_text(text)
    digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    provenance={'method':'official Rosetta molfile_to_params; original baseline conformers',
        'source_sha256':digest(sdf_file),'generator_sha256':digest(script),'params_sha256':digest(params),
        'pdb_sha256':digest(pdb),'input_conformers':len(hyd),'heavy_atoms':molecules[0].GetNumAtoms(),
        'all_atom_index_to_name':mapping,'input_coordinate_frames_preserved':True,
        'neighbor_atom_zero_based':nbr,'rotamer_frame_atom_indices':frame,
        'rotamer_frame_max_displacement_A':displacement if frame is not None else None,
        'virtual_atoms_requested':bool(v_atoms_sdf),'virtual_atoms_used':False,
        'virtual_atom_note':'Compatibility reference only; fixed-head construction uses real atoms, and no physical vanadium atoms are added',
        'command':command}
    prefix.with_name(prefix.name+'_parameter_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    requested=Path(pdb_prefix)
    return str(requested)+('.pdb' if conformers else '_0001.pdb'),str(requested)+'.params'

def docked_heads(docked_path,head_paths,root=None):
    pdb=Path(str(docked_path).replace('.sdf','.pdb'))
    root=Path(root) if root is not None else Path(os.environ['PROTAC_ATTEMPT_ROOT'])
    records={}
    for line in pdb.read_text().splitlines():
        if line.startswith('HETATM') and line[17:20].strip() in {'PT0','PT1'}:
            key=(line[17:20].strip(),line[12:16].strip())
            if key in records:raise RuntimeError('Duplicate docked head atom name')
            records[key]=list(map(float,[line[30:38],line[38:46],line[46:54]]))
    result=[]
    for role,path in enumerate(head_paths):
        residue='PT'+str(role)
        meta=json.loads((root/(residue+'_parameter_provenance.json')).read_text())
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=meta['source_sha256']:
            raise RuntimeError('Head template differs from parameter provenance')
        molecule=Chem.SDMolSupplier(str(path),removeHs=True)[0]
        if molecule is None or molecule.GetNumAtoms()!=meta['heavy_atoms']:raise RuntimeError('Head graph/count mismatch')
        conf=Chem.Conformer(molecule.GetNumAtoms())
        for i in range(molecule.GetNumAtoms()):
            name=meta['all_atom_index_to_name'][str(i)]
            if (residue,name) not in records:raise RuntimeError('Docking lost a mapped head atom')
            conf.SetAtomPosition(i,records[(residue,name)])
        molecule.RemoveAllConformers();molecule.AddConformer(conf);result.append(molecule)
    return Chem.CombineMols(*result)
