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

DEFAULT_EMBED_ITERATIONS=20
DEFAULT_OUTER_ATTEMPTS=20

def fixed_head_neighbor(full,head_indices):
    """Choose a rotamer orientation center whose heavy neighbors are all fixed.

    Rosetta orients ligand conformers using NBR_ATOM and two bonded heavy atoms.
    A linker attachment atom can therefore rotate a conformer off both heads.
    A ring atom wholly inside the fixed head supplies an invariant frame.
    """
    fixed=set(head_indices)
    boundary=[i for i in fixed if any(n.GetIdx() not in fixed for n in full.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum()>1)]
    candidates=[]
    distances=Chem.GetDistanceMatrix(full)
    for i in sorted(fixed):
        atom=full.GetAtomWithIdx(i)
        neighbors={n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum()>1}
        if atom.IsInRing() and len(neighbors)>=2 and neighbors<=fixed:
            separation=min((distances[i,j] for j in boundary),default=0.)
            candidates.append((-separation,i))
    if not candidates:raise ValidationError('No invariant ring frame inside fixed head')
    return min(candidates)[1]

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

def distance_anchor_indices(full,heads,anchors):
    """Require exact disjoint stereochemical heads and unambiguous full anchors."""
    if len(heads)!=2 or len(anchors)!=2:
        raise ValidationError('Exactly two heads and anchors are required')
    if any(not 0<=i<h.GetNumAtoms() for h,i in zip(heads,anchors)):
        raise ValidationError('Anchor outside its exact heavy-atom head')
    pairs={(a[anchors[0]],b[anchors[1]]) for a,b in match_heads(full,heads)}
    if len(pairs)!=1:
        raise ValidationError('Ambiguous full-molecule anchor mapping')
    return next(iter(pairs))


def fit(mobile,reference):
    mc,rc=mobile.mean(0),reference.mean(0)
    u,_,vt=np.linalg.svd((mobile-mc).T@(reference-rc))
    rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    return rotation,rc-mc@rotation

def constrained_linker(full,heads,targets,seed,attempts=DEFAULT_OUTER_ATTEMPTS,rmsd_limit=0.35,embed_iterations=DEFAULT_EMBED_ITERATIONS):
    """Return a valid full PROTAC or an explicit geometric sampling rejection.

    Both complete bound head fragments participate in the distance constraints.
    All head heavy atoms stay fixed during minimization. A failed embedding never
    falls back to a free conformer or an unrelated VHL template.
    """
    if not 1<=attempts<=2000 or not 1<=embed_iterations<=2000:
        raise ValidationError('Explicit positive bounded embedding budgets are required')
    origin=np.vstack(targets).mean(axis=0)
    targets=[np.asarray(a,dtype=float)-origin for a in targets]
    expected=Chem.MolToInchiKey(full)
    mappings=match_heads(full,heads)
    failures=[];embedding_failures=[]
    for trial in range(attempts):
        a,b=mappings[trial%len(mappings)]
        indices=list(a+b);target=np.vstack(targets)
        mol=Chem.AddHs(Chem.Mol(full))
        params=AllChem.ETKDGv3()
        params.useExpTorsionAnglePrefs=True
        params.useBasicKnowledge=False
        params.randomSeed=(int(seed)+trial)%2147483647
        params.enforceChirality=True
        params.numThreads=1
        # RDKit's zero default means 10*all-atom-count INTERNAL retries per
        # call, multiplied by the outer seed budget. Make this hidden cost
        # explicit and deterministic; never reinterpret exhaustion as inactivity.
        params.maxIterations=int(embed_iterations)
        params.trackFailures=True
        params.useRandomCoords=True
        params.SetCoordMap({idx:Point3D(*map(float,p)) for idx,p in zip(indices,target)})
        embedded=AllChem.EmbedMolecule(mol,params)
        embedding_failures.append(list(params.GetFailureCounts()))
        if embedded<0:
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
        for i,point in enumerate(mol.GetConformer().GetPositions()+origin):mol.GetConformer().SetAtomPosition(i,point)
        return mol,{'embedding_method':'CENTERED_TORSION_DG_FIXED_HEADS_UFF','embedding_origin_A':origin.tolist(),'accepted':True,'seed':params.randomSeed,'attempts':trial+1,
                    'head_mappings':[list(a),list(b)],'mapping_alternatives':len(mappings),
                    'initial_both_head_rmsd_A':float(np.sqrt(np.mean(error**2))),
                    'embedding_max_iterations_per_attempt':embed_iterations,
                    'embedding_failure_counts_by_attempt':embedding_failures,
                    'rejection_reasons':failures}
    return None,{'embedding_method':'CENTERED_TORSION_DG_FIXED_HEADS_UFF','embedding_origin_A':origin.tolist(),'accepted':False,'seed':int(seed),'attempts':attempts,
                 'mapping_alternatives':len(mappings),'rejection_reasons':failures,
                 'embedding_max_iterations_per_attempt':embed_iterations,
                 'embedding_failure_counts_by_attempt':embedding_failures}

