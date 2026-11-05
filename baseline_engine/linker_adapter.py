"""Bound-head-preserving linker construction with explicit sampling rejections."""
import hashlib,json,os
from pathlib import Path
from rdkit import Chem
from rdkit.Geometry import Point3D
from constrained_linker import constrained_linker,match_heads,fixed_head_neighbor,DEFAULT_OUTER_ATTEMPTS
from official_parameters import docked_heads

def generate(head_paths,docked,smiles_path,output,anchor_a,virtual_path,n=100,homo_protac=False):
    from pyrosetta_adapter import publication_seed
    heads=[Chem.SDMolSupplier(str(p),removeHs=True)[0] for p in head_paths]
    if any(h is None for h in heads):raise RuntimeError('Invalid bound head graph')
    reference=docked_heads(docked,head_paths)
    positions=reference.GetConformer().GetPositions();cut=heads[0].GetNumAtoms()
    targets=[positions[:cut],positions[cut:]]
    full=Chem.MolFromSmiles(Path(smiles_path).read_text().split()[0])
    if full is None:raise RuntimeError('Invalid full PROTAC graph')
    mappings=match_heads(full,heads)
    neighbour=fixed_head_neighbor(full,mappings[0][0])
    # The original caller creates/removes a virtual-reference SDF. Keep that
    # interface artifact valid; no dummy atoms enter physical ligand parameters.
    rw=Chem.RWMol()
    for _ in range(3):rw.AddAtom(Chem.Atom(0))
    virtual=rw.GetMol();conf=Chem.Conformer(3);center=positions.mean(0)
    for i,delta in enumerate([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]):conf.SetAtomPosition(i,Point3D(*map(float,center+delta)))
    virtual.AddConformer(conf);writer=Chem.SDWriter(str(virtual_path));writer.write(virtual);writer.close()
    requested=min(n,int(os.environ.get('PROTAC_CONFORMERS',str(n))))
    if os.environ.get('PROTAC_DIAGNOSTIC')=='1':requested=min(requested,int(os.environ.get('PROTAC_DIAGNOSTIC_CONFORMERS','2')))
    trials=max(1,int(os.environ.get('PROTAC_LINKER_EMBED_TRIALS',str(DEFAULT_OUTER_ATTEMPTS))))
    diagnostics=[];accepted=0
    writer=Chem.SDWriter(str(output));writer.SetKekulize(False)
    for i in range(requested):
        seed=publication_seed('fixed_head_linker_'+str(i))
        molecule,result=constrained_linker(full,heads,targets,seed,attempts=trials)
        diagnostics.append(result)
        if molecule is None:
            # All requested conformers have identical fixed-head geometry. A
            # bounded failed feasibility search is a rejection, not a fallback.
            if accepted==0:break
            continue
        writer.write(molecule);accepted+=1
    writer.close()
    report={'status':'FEASIBLE_CONFORMERS_GENERATED' if accepted else 'GEOMETRY_SAMPLING_REJECTED',
        'requested_conformers':requested,'generated_conformers':accepted,'embedding_trials_per_search':trials,
        'conformer_searches_attempted':len(diagnostics),
        'conformer_searches_not_attempted':requested-len(diagnostics),
        'early_stop_after_failed_initial_feasibility_search':accepted==0,
        'bound_head_atoms_fixed':True,'fallback_to_free_conformer':False,
        'rotamer_neighbor_atom_zero_based':neighbour,
        'rotamer_orientation_frame':'Ring center and bonded heavy neighbors wholly inside fixed VHL head',
        'full_molecule_inchikey':Chem.MolToInchiKey(full),'diagnostics':diagnostics,
        'interpretation':'A failed bounded search does not establish biological nondegradation',
        'source_sha256':{str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in [*head_paths,smiles_path]}}
    Path(str(output)+'.sampling.json').write_text(json.dumps(report,indent=2)+'\n')
    return neighbour,virtual_path
