import importlib.util,json,unittest
from pathlib import Path
import numpy as np
from rdkit import Chem
BASE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('ms143',BASE/'scripts/02_prepare_ms143_heads.py')
ms=importlib.util.module_from_spec(spec);spec.loader.exec_module(ms)

class MS143GeometryChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=BASE/'MS143/strict_review_v2'
        cls.meta=json.loads((cls.path/'PREPARATION.json').read_text())
        cls.parent=Chem.SDMolSupplier(str(cls.path/'AZD5363_bound_parent.sdf'),removeHs=True)[0]
        cls.head=Chem.SDMolSupplier(str(cls.path/'poi_modeled_head.sdf'),removeHs=True)[0]
        cls.protein=cls.path/'AKT1_parent_protein.pdb'
    def test_parent_and_selected_head_pass_declared_steric_screen(self):
        self.assertLessEqual(ms.overlap_review(self.parent,self.protein)['max_vdw_overlap_A'],.6)
        self.assertLessEqual(ms.overlap_review(self.head,self.protein)['max_vdw_overlap_A'],.6)
    def test_initial_coarse_candidate_is_rejected_by_steric_screen(self):
        bad=Chem.SDMolSupplier(str(BASE/'MS143/poi_modeled_head.sdf'),removeHs=True)[0]
        self.assertGreater(ms.overlap_review(bad,self.protein)['max_vdw_overlap_A'],1.)
    def test_bound_scaffold_and_head_3d_stereo_are_preserved(self):
        meta=self.meta['poi_model'];p=self.parent.GetConformer().GetPositions()[meta['core_reference_indices_zero_based']]
        q=self.head.GetConformer().GetPositions()[meta['core_head_indices_zero_based']]
        self.assertEqual(len(p),29);self.assertLess(float(np.linalg.norm(p-q,axis=1).max()),.001)
        expected=Chem.MolToInchiKey(self.head);check=Chem.Mol(self.head)
        Chem.AssignStereochemistryFrom3D(check,replaceExistingTags=True)
        self.assertEqual(Chem.MolToInchiKey(check),expected)
        self.assertFalse(self.meta['production_ready'])

if __name__=='__main__':unittest.main()
