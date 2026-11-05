#!/usr/bin/env python3
"""Model modified conjugated heads using fixed bound common scaffolds.

Coordinates outside the shared scaffold are explicitly modeled, never labelled
experimental. These analog poses remain blocked pending structural validation.
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError, pdb_atoms, min_distance, atomic_json, sha256

def fragment(mol,ids):
    rw=Chem.RWMol(mol)
    for i in reversed(range(mol.GetNumAtoms())):
        if i not in ids:rw.RemoveAtom(i)
    part=rw.GetMol();Chem.SanitizeMol(part)
    return part

def extract_head(full,name):
    if name=='SJF1521':query=Chem.MolFromSmarts('[c][O][CH2][CH2][O]');size=34
    elif name=='JNJ1013':query=Chem.MolFromSmarts('[C;R][O][CH2][C](=[O])[N]');size=29
    elif name=='SJF_delta':query=Chem.MolFromSmarts('[c][O][CH2][CH2][O]');size=39
    else:raise ValidationError('Unsupported modified head')
    matches=full.GetSubstructMatches(query)
    expected=2 if name=='SJF_delta' else 1
    if len(matches)!=expected:raise ValidationError('Ambiguous supplier attachment pattern')
    rw=Chem.RWMol(full)
    cuts=sorted({tuple(sorted((match[1],match[2]))) for match in matches})
    for a,b in cuts:
        if rw.GetBondBetweenAtoms(a,b).IsInRing():raise ValidationError('Cannot cut a ring')
        rw.RemoveBond(a,b)
    split=rw.GetMol();Chem.SanitizeMol(split)
    candidates=[ids for ids in Chem.GetMolFrags(split) if len(ids)==size]
    if len(candidates)!=1:raise ValidationError('Unexpected head partition')
    ids=candidates[0];head=fragment(full,ids)
    exit_ids=[a for cut in cuts for a in cut if a in ids]
    if len(exit_ids)!=1:raise ValidationError('Ambiguous head exit atom')
    return head,list(ids).index(exit_ids[0]),list(ids),cuts

def model(head,reference,protein,seeds):
    parameters=rdFMCS.MCSParameters()
    parameters.Timeout=20
    parameters.AtomTyper=rdFMCS.AtomCompare.CompareElements
    parameters.BondTyper=rdFMCS.BondCompare.CompareOrderExact
    parameters.AtomCompareParameters.RingMatchesRingOnly=True
    parameters.AtomCompareParameters.CompleteRingsOnly=True
    parameters.BondCompareParameters.RingMatchesRingOnly=True
    parameters.BondCompareParameters.CompleteRingsOnly=True
    mcs=rdFMCS.FindMCS([head,reference],parameters)
    if mcs.canceled or mcs.numAtoms<20:raise ValidationError('Insufficient exact common scaffold')
    query=Chem.MolFromSmarts(mcs.smartsString)
    ref_matches=reference.GetSubstructMatches(query)
    if not ref_matches:raise ValidationError('Cannot map bound scaffold')
    # Retain the complete shared scaffold and its deposited coordinates.
    ids=sorted(ref_matches[0]);core=fragment(reference,ids)
    target_matches=head.GetSubstructMatches(core)
    if not target_matches:raise ValidationError('Common scaffold cannot map to conjugated head')
    match=target_matches[0]
    protein_xyz=np.array([a['xyz'] for a in pdb_atoms(protein) if a['element'] not in {'H','D'}])
    expected_key=Chem.MolToInchiKey(head)
    candidates=[];failures=[]
    for seed in seeds:
        try:
            embedded=AllChem.ConstrainedEmbed(Chem.AddHs(head),core,randomseed=seed,useTethers=True)
            for i,j in enumerate(match):
                embedded.GetConformer().SetAtomPosition(j,core.GetConformer().GetAtomPosition(i))
            ff=AllChem.UFFGetMoleculeForceField(embedded)
            for j in match:ff.AddFixedPoint(j)
            converged=ff.Minimize(maxIts=2000)
            if converged:raise ValidationError('UFF did not converge')
            energy=float(ff.CalcEnergy())
            if not np.isfinite(energy):raise ValidationError('Nonfinite UFF energy')
            embedded=Chem.RemoveHs(embedded)
            Chem.AssignStereochemistryFrom3D(embedded,replaceExistingTags=True)
            if Chem.MolToInchiKey(embedded)!=expected_key:raise ValidationError('Modeled stereochemistry changed')
            distance=min_distance(embedded.GetConformer().GetPositions(),protein_xyz)
            if distance<1.5:raise ValidationError('Severe protein clash')
            core_error=float(np.max(np.linalg.norm(embedded.GetConformer().GetPositions()[list(match)]-
                                       core.GetConformer().GetPositions(),axis=1)))
            if core_error>1e-6:raise ValidationError('Bound scaffold moved during minimization')
            candidates.append((energy,seed,distance,embedded))
        except (ValueError,RuntimeError) as exc:failures.append({'seed':seed,'reason':str(exc)})
    if not candidates:raise ValidationError('No converged clash-free scaffold-preserving conformer; '+json.dumps(failures))
    energy,seed,distance,chosen=min(candidates,key=lambda x:(x[0],x[1]))
    meta={'method':'fixed_bound_scaffold_analog_model','experimental_coordinates_for_unmatched_atoms':False,
          'production_ready':False,'common_scaffold_atoms':len(match),'head_atoms':head.GetNumAtoms(),
          'core_reference_indices_zero_based':ids,'core_head_indices_zero_based':list(match),
          'reference_match_count':len(ref_matches),'head_match_count':len(target_matches),
          'symmetry_rule':'first RDKit exact graph match; all counts archived for review',
          'mcs_smarts':mcs.smartsString,'seeds_attempted':list(seeds),'seeds_valid':len(candidates),
          'selected_seed':seed,'UFF_energy_internal_only':energy,'min_protein_contact_A':distance,
          'failed_seeds':failures,'review_required':['modified_atoms_pose','symmetry_mapping','binding_mode','protein_environment_not_in_UFF_energy']}
    return chosen,meta

def prepare(root):
    refs=root/'reference_structures/validated'
    registry=json.loads((root/'config/molecules.json').read_text())['compounds']
    results={}
    for name,parent,role in [('SJF1521','SJF1521_parent','poi'),('JNJ1013','JNJ1013_parent','poi'),('SJF_delta','VH032_parent','vhl')]:
        full=Chem.MolFromSmiles(registry[name]['smiles'])
        head,anchor,ids,cuts=extract_head(full,name)
        reference=Chem.MolFromMolFile(str(refs/parent/'bound_ligand.sdf'),removeHs=False)
        destination=root/'modeled_reference_candidates'/name;destination.mkdir(parents=True,exist_ok=True)
        try:
            chosen,meta=model(head,reference,refs/parent/'proteins.pdb',range(1,13))
            meta.update(full_supplier_indices_zero_based=ids,full_supplier_cut_bonds_zero_based=cuts,
                        supplier_inchikey=registry[name]['inchikey'],anchor_one_based=anchor+1,
                        parent_sha256=sha256(refs/parent/'bound_ligand.sdf'))
            file=destination/(role+'_modeled_head.sdf')
            w=Chem.SDWriter(str(file));w.write(chosen);w.close()
            meta['head_sha256']=sha256(file)
            atomic_json(destination/(role+'_model_provenance.json'),meta)
            results[name]={'status':'MODELED_CANDIDATE_REQUIRES_REVIEW','role':role,
                           'common_atoms':meta['common_scaffold_atoms'],'valid_seeds':meta['seeds_valid']}
        except (ValueError,RuntimeError) as exc:
            results[name]={'status':'REJECTED','reason':str(exc)}
    atomic_json(root/'audits/modified_head_modeling.json',results)
    # Assemble complete candidate bundles while retaining the per-head provenance.
    # No production workflow is allowed to consume these implicitly.
    bundles={
        'SJF1521':('SJF1521_parent','VH032_parent','SJF8240'),
        'JNJ1013':('JNJ1013_parent','GSK215','GSK215'),
        'SJF_delta':('foretinib_parent','VH032_parent','SJF8240'),
    }
    for name,(poi_parent,vbc_parent,exact_head_parent) in bundles.items():
        if results[name]['status']=='REJECTED':continue
        import shutil
        destination=root/'modeled_reference_candidates'/name
        modeled_role=results[name]['role']
        other_role='vhl' if modeled_role=='poi' else 'poi'
        exact=json.loads((refs/exact_head_parent/'heads.json').read_text())[other_role]
        shutil.copy2(refs/exact_head_parent/(other_role+'_head.sdf'),destination/(other_role+'_head.sdf'))
        shutil.copy2(destination/(modeled_role+'_modeled_head.sdf'),destination/(modeled_role+'_head.sdf'))
        shutil.copy2(refs/poi_parent/'protein.pdb',destination/'protein.pdb')
        shutil.copy2(refs/vbc_parent/'proteins.pdb',destination/'proteins.pdb')
        modeled=json.loads((destination/(modeled_role+'_model_provenance.json')).read_text())
        heads={other_role:exact,modeled_role:{'anchor_one_based':modeled['anchor_one_based']},
               'production_ready':False,'method':'modeled_conjugated_head_candidate',
               'modeled_role':modeled_role,'model_provenance_sha256':sha256(destination/(modeled_role+'_model_provenance.json'))}
        for role in ['poi','vhl']:
            head=Chem.MolFromMolFile(str(destination/(role+'_head.sdf')))
            full=Chem.MolFromSmiles(registry[name]['smiles'])
            if not full.HasSubstructMatch(head,useChirality=True):raise ValidationError('Bundle head chemistry mismatch: '+name)
            heads[role]['file']=str((destination/(role+'_head.sdf')).relative_to(root))
            heads[role]['sha256']=sha256(destination/(role+'_head.sdf'))
        atomic_json(destination/'heads.json',heads)
    print(json.dumps(results,indent=2))
    if any(r['status']=='REJECTED' for r in results.values()):raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    prepare(p.parse_args().root)
