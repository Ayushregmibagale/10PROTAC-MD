from pathlib import Path
import importlib.util,json,sys
import numpy as np
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'candidate_workflow/baseline_engine';sys.path.insert(0,str(BASE))
spec=importlib.util.spec_from_file_location('candidate_bound_head_DG',BASE/'constrained_linker.py');method=importlib.util.module_from_spec(spec);spec.loader.exec_module(method)
FIX=Path(__file__).parent/'fixtures'

def inputs(frame='10_1'):
    full=Chem.MolFromSmiles((FIX/'full.smi').read_text().split()[0]);heads=[Chem.SDMolSupplier(str(FIX/n),removeHs=True)[0] for n in ['vhl_head_H.sdf','poi_head_H.sdf']];p=np.asarray(json.loads((FIX/'head_coordinates.json').read_text())[frame]['coordinates_A']);return full,heads,[p[:30],p[30:]]

def assert_fixed(full,heads,targets,mol):
    mapping=method.match_heads(full,heads)[0];xyz=mol.GetConformer().GetPositions()
    for ids,expected in zip(mapping,targets):assert np.max(np.linalg.norm(xyz[list(ids)]-expected,axis=1))<1e-5
    assert method.stereo_key(mol)==Chem.MolToInchiKey(full)

def test_bound_head_geometry_preserved_in_original_protein_frame():
    full,heads,targets=inputs();mol,report=method.constrained_linker(full,heads,targets,155921,attempts=3)
    assert mol is not None and report['embedding_method']=='CENTERED_DG_FIXED_HEADS_UFF';assert_fixed(full,heads,targets,mol)

def test_large_translation_does_not_destroy_a_feasible_head_geometry():
    full,heads,targets=inputs();shift=np.array([1234.,-987.,456.]);targets=[p+shift for p in targets]
    mol,report=method.constrained_linker(full,heads,targets,155921,attempts=3)
    assert mol is not None;assert_fixed(full,heads,targets,mol)

def test_failed_geometry_remains_an_explicit_bounded_rejection():
    full,heads,targets=inputs('10_3');mol,report=method.constrained_linker(full,heads,targets,155921,attempts=3)
    assert mol is None and not report['accepted'] and report['attempts']==3
    assert report['embedding_method']=='CENTERED_DG_FIXED_HEADS_UFF'


def test_torsion_preserving_variant_keeps_exact_bound_head_geometry():
    path=ROOT/'candidate_torsion_workflow/baseline_engine/constrained_linker.py'
    spec=importlib.util.spec_from_file_location('candidate_torsion_bound_heads',path)
    torsion=importlib.util.module_from_spec(spec);spec.loader.exec_module(torsion)
    full,heads,targets=inputs()
    for shift in [np.zeros(3),np.array([1234.,-987.,456.])]:
        moved=[p+shift for p in targets]
        mol,report=torsion.constrained_linker(full,heads,moved,155921,attempts=3)
        assert mol is not None and report['embedding_method']=='CENTERED_TORSION_DG_FIXED_HEADS_UFF'
        assert_fixed(full,heads,moved,mol)
