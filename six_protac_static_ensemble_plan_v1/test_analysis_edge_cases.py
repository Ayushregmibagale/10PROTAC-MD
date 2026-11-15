"""Scientific reporting edge cases: zero valid poses and explicit chain roles."""
import importlib.util,json,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
from rdkit import Chem
BASE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('coverage',BASE/'scripts/04_audit_prepared_coverage.py')
coverage=importlib.util.module_from_spec(spec);spec.loader.exec_module(coverage)

class AnalysisEdgeCases(unittest.TestCase):
    def test_no_valid_poses_is_reported_without_invented_recovery(self):
        with tempfile.TemporaryDirectory(prefix='six-static-nohit-') as tmp:
            root=Path(tmp)
            for ref in ('5T35_MZ1','7PI4_GSK215'):
                run=root/'runs'/f'global_{ref}_seed42';(run/'Patchdock_Results').mkdir(parents=True)
                shutil.copy2(BASE/'reference_inputs'/ref/'selection.json',run/'selection.json')
                (run/'Patchdock_Results/score.sc').write_text('SEQUENCE: \nSCORE: total_score description\n')
                (run/'engine_status.json').write_text(json.dumps({'status':'COMPLETE_NO_VALID_POSES',
                    'stages':{'scoring':{'succeeded':0},'pose_qc':{'succeeded':0}},'started_unix':0,'finished_unix':1}))
            proc=subprocess.run([sys.executable,str(BASE/'scripts/01_assess_reference_recovery.py'),
                '--runs',str(root/'runs'),'--references',str(BASE/'reference_evaluation'),'--output',str(root/'output')],
                capture_output=True,text=True,timeout=30)
            self.assertEqual(proc.returncode,0,proc.stderr)
            report=json.loads((root/'output/REFERENCE_RECOVERY.json').read_text())
            for row in report['summary']:
                self.assertEqual(row['retained_QC_poses'],0);self.assertEqual(row['gate_pass_count'],0)
                self.assertIsNone(row['best_POI_RMSD_A']);self.assertFalse(row['biological_validation'])
            self.assertEqual(len((root/'output/REFERENCE_POSE_RECOVERY.tsv').read_text().splitlines()),1)
    def test_partner_chain_is_not_mistaken_for_vhl_binary_contacts(self):
        source=BASE/'reference_inputs/7PI4_GSK215';head=source/'vhl_head.sdf'
        mol=Chem.SDMolSupplier(str(head),removeHs=True)[0];x,y,z=mol.GetConformer().GetPositions()[0]
        ghost=f'ATOM  99999  CA  ALA D9999    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  '
        original='\n'.join(s for s in (source/'vbc.pdb').read_text().splitlines() if s.startswith('ATOM  '))
        with tempfile.TemporaryDirectory(prefix='six-static-roles-') as tmp:
            protein=Path(tmp)/'protein.pdb';protein.write_text(original+'\n'+ghost+'\nEND\n')
            restricted=coverage.sterics(head,protein,chains={'A','B','C'})
            allchains=coverage.sterics(head,protein)
            self.assertGreater(allchains['max_vdw_overlap_A'],3.)
            self.assertLess(restricted['max_vdw_overlap_A'],.8)

if __name__=='__main__':unittest.main()
