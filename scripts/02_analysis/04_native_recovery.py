#!/usr/bin/env python3
"""Measure protein-interface recovery after receptor alignment.

These are explicit project metrics, not DockQ/CAPRI classifications, free energies
or validation of a biological prediction. A native-frame control is not an
independent global redocking benchmark. Residue numbering must be compatible.
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,pdb_atoms,sha256,atomic_json
spec=importlib.util.spec_from_file_location('recovery_alignment',Path(__file__).resolve().parents[1]/'01_setup/01_place_warhead_batch.py')
alignment=importlib.util.module_from_spec(spec);spec.loader.exec_module(alignment)


def residue(atom):return (atom['chain'],atom['resnum'],atom['icode'],atom['resname'])


def ca_map(atoms,chains):
    result={}
    for atom in atoms:
        if atom['record']=='ATOM' and atom['name']=='CA' and atom['chain'] in chains:
            key=residue(atom)
            if key in result:raise ValidationError('Ambiguous Cα residue identity')
            result[key]=np.array(atom['xyz'])
    return result


def contacts(atoms,poi,vhl,cutoff):
    def selected(chains):
        return [a for a in atoms if a['record']=='ATOM' and a['chain'] in chains and a['element'] not in {'H','D'}]
    receptor,mobile=selected(poi),selected(vhl)
    if not receptor or not mobile:raise ValidationError('Missing interface protein atoms')
    tree=cKDTree(np.array([a['xyz'] for a in mobile]))
    neighbors=tree.query_ball_point(np.array([a['xyz'] for a in receptor]),cutoff)
    return {(residue(a),residue(mobile[j])) for a,indices in zip(receptor,neighbors) for j in indices}


def measure(reference,pose,selection,cutoff=5.0):
    if not np.isfinite(cutoff) or cutoff<=0:raise ValidationError('Contact cutoff must be finite and positive')
    ref,model=pdb_atoms(reference),pdb_atoms(pose)
    poi=set(selection['poi']['chains']);vhl=set(selection['vhl']['chains'])
    if poi&vhl:raise ValidationError('Receptor and mobile selections must be disjoint')
    paired={};coverage={}
    for role,chains,minimum in [('poi',poi,100),('vhl',vhl,50)]:
        r,m=ca_map(ref,chains),ca_map(model,chains);keys=sorted(r.keys()&m.keys())
        if len(keys)<minimum or len(keys)/max(len(r),1)<0.8:
            raise ValidationError('Insufficient compatible reference Cα coverage for '+role)
        paired[role]=(np.array([m[k] for k in keys]),np.array([r[k] for k in keys]))
        coverage[role]={'matched_CA':len(keys),'reference_CA':len(r),'fraction':len(keys)/len(r)}
    rotation,translation=alignment.kabsch(*paired['poi'])
    rmsd={role:float(np.sqrt(np.mean(np.sum((m@rotation.T+translation-r)**2,axis=1))))
          for role,(m,r) in paired.items()}
    native_contacts=contacts(ref,poi,vhl,cutoff);model_contacts=contacts(model,poi,vhl,cutoff)
    if not native_contacts:raise ValidationError('Reference has no protein-interface contacts at this cutoff')
    recovered=len(native_contacts&model_contacts)
    return {'status':'GEOMETRY_BENCHMARK_ONLY','reference_sha256':sha256(reference),'pose_sha256':sha256(pose),
        'selection':selection,'alignment':'POI Cα, compatible chain/residue number/insertion/residue name',
        'coverage':coverage,'poi_aligned_CA_rmsd_A':rmsd['poi'],'vhl_CA_rmsd_after_POI_fit_A':rmsd['vhl'],
        'contact_cutoff_A':cutoff,'contacts_definition':'Unique heavy-atom protein residue pairs within cutoff',
        'native_contact_pairs':len(native_contacts),'model_contact_pairs':len(model_contacts),
        'recovered_native_contact_pairs':recovered,'fraction_native_contacts':recovered/len(native_contacts),
        'independent_global_recovery_demonstrated':False,
        'limitations':['Protein orientation/interface metrics only; ligand geometry and chemistry need separate QC',
            'A native-frame or input-derived self comparison cannot establish independent predictive recovery',
            'No DockQ/CAPRI classification, affinity, degradation or thermodynamic probability inferred']}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['reference','pose','selection','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--contact-cutoff',type=float,default=5.0)
    a=p.parse_args();result=measure(a.reference,a.pose,json.loads(a.selection.read_text()),a.contact_cutoff)
    atomic_json(a.output,result);print(json.dumps(result,indent=2))
