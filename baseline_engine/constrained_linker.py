"""Exact two-head linker constraints for the existing pRosettaC adapter.

Restores bound-head geometry during linker generation. No global/local protein
sampling or docking algorithm is replaced.
"""
import sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/lib'))
try:
    from protac_qc import ValidationError
except ImportError:
    class ValidationError(ValueError):pass

def stereo_key(mol):
    mol=Chem.RemoveHs(Chem.Mol(mol))
    Chem.AssignStereochemistryFrom3D(mol,replaceExistingTags=True)
    return Chem.MolToInchiKey(mol)

def match_heads(full,heads):
    """Full stereochemical substructures, disjoint atom sets; never MCS substitution."""
    matches=[full.GetSubstructMatches(h,useChirality=True,uniquify=True) for h in heads]
    pairs=[(tuple(a),tuple(b)) for a in matches[0] for b in matches[1] if not set(a)&set(b)]
    if not pairs:raise ValidationError('No disjoint stereochemical full-molecule head mapping')
    # Graph automorphisms are deterministic, with all alternatives retained for provenance.
    return sorted(set(pairs))

def fit(mobile,reference):
    mc,rc=mobile.mean(0),reference.mean(0)
    u,_,vt=np.linalg.svd((mobile-mc).T@(reference-rc))
    rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    return rotation,rc-mc@rotation

def constrained_linker(full,heads,targets,seed,attempts=20,rmsd_limit=0.35):
    """Return a valid full PROTAC or an explicit geometric sampling rejection.

    Both complete bound head fragments participate in the distance constraints.
    All head heavy atoms stay fixed during minimization. A failed embedding never
    falls back to a free conformer or an unrelated VHL template.
    """
    expected=Chem.MolToInchiKey(full)
    mappings=match_heads(full,heads)
    failures=[]
    for trial in range(attempts):
        a,b=mappings[trial%len(mappings)]
        indices=list(a+b);target=np.vstack(targets)
        mol=Chem.AddHs(Chem.Mol(full))
        params=AllChem.ETKDGv3()
        params.randomSeed=(int(seed)+trial)%2147483647
        params.enforceChirality=True
        params.numThreads=1
        params.useRandomCoords=True
        params.SetCoordMap({idx:Point3D(*map(float,p)) for idx,p in zip(indices,target)})
        if AllChem.EmbedMolecule(mol,params)<0:
            failures.append('embedding_failed');continue
        conf=mol.GetConformer();positions=conf.GetPositions()
        r,t=fit(positions[indices],target)
        positions=positions@r+t
        error=np.linalg.norm(positions[indices]-target,axis=1)
        if float(np.sqrt(np.mean(error**2)))>rmsd_limit or error.max()>1.0:
            failures.append('both_head_rmsd_exceeded');continue
        for i,p in enumerate(positions):conf.SetAtomPosition(i,p)
        for i,p in zip(indices,target):conf.SetAtomPosition(i,p)
        ff=AllChem.UFFGetMoleculeForceField(mol)
        if ff is None:raise ValidationError('Missing UFF parameters')
        for i in indices:ff.AddFixedPoint(i)
        if ff.Minimize(maxIts=2000)!=0:
            failures.append('linker_minimization_not_converged');continue
        mol=Chem.RemoveHs(mol)
        if stereo_key(mol)!=expected:
            failures.append('output_stereochemistry_changed');continue
        pos=mol.GetConformer().GetPositions()
        if np.max(np.linalg.norm(pos[indices]-target,axis=1))>1e-5:
            raise ValidationError('Constrained head atoms moved')
        # Reject stretched/compressed bonds before calling a protein score function.
        if any(not 0.9<np.linalg.norm(pos[x.GetBeginAtomIdx()]-pos[x.GetEndAtomIdx()])<2.1
               for x in mol.GetBonds()):
            failures.append('distorted_linker_bonds');continue
        return mol,{'accepted':True,'seed':params.randomSeed,'attempts':trial+1,
                    'head_mappings':[list(a),list(b)],'mapping_alternatives':len(mappings),
                    'initial_both_head_rmsd_A':float(np.sqrt(np.mean(error**2))),
                    'rejection_reasons':failures}
    return None,{'accepted':False,'seed':int(seed),'attempts':attempts,
                 'mapping_alternatives':len(mappings),'rejection_reasons':failures}

