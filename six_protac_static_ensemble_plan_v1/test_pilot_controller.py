import hashlib,json,subprocess,sys,tempfile,unittest,os
from pathlib import Path
CONTROLLER=Path(__file__).resolve().parent/'scripts/03_run_ernie_reference_pilot.py'

class ControllerChecks(unittest.TestCase):
    def fixture(self,mode):
        tmp=tempfile.TemporaryDirectory(prefix='sixprotac-controller-');self.addCleanup(tmp.cleanup);r=Path(tmp.name)
        (r/'workflow/baseline_engine').mkdir(parents=True);(r/'inputs/ref').mkdir(parents=True)
        fake='''import json,sys,time,subprocess
from pathlib import Path
r=Path(sys.argv[sys.argv.index('--run-dir')+1])
'''
        if mode=='complete':fake+="(r/'engine_status.json').write_text(json.dumps({'status':'COMPLETE','stages':{}}))\n"
        if mode=='failure':fake+="raise SystemExit(7)\n"
        if mode=='timeout':fake+="p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)']); (r/'child.pid').write_text(str(p.pid)); time.sleep(60)\n"
        engine=r/'workflow/baseline_engine/run.py';engine.write_text(fake)
        cfg={'cpu_limit':1,'poll_seconds':.1,'max_group_RSS_bytes':1024**3,'max_output_bytes':1024**3,'min_free_bytes':0,
             'tasks':[{'name':'test','reference':'ref','seed':42,'protocol':'global','n_dock':1,'n_local':1,'conformers':1,'max_global_trials':10,'timeout_seconds':.5 if mode=='timeout' else 5}]}
        (r/'PILOT_CONFIG.json').write_text(json.dumps(cfg))
        (r/'BUNDLE_MANIFEST.json').write_text(json.dumps({'files':{'workflow/baseline_engine/run.py':hashlib.sha256(engine.read_bytes()).hexdigest()}}))
        return r
    def run_controller(self,r):return subprocess.run([sys.executable,str(CONTROLLER),'--root',str(r)],capture_output=True,text=True,timeout=20)
    def test_success_receipt(self):
        r=self.fixture('complete');p=self.run_controller(r);self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(json.loads((r/'receipts/test.json').read_text())['status'],'COMPLETE')
    def test_nonzero_is_failure(self):
        r=self.fixture('failure');p=self.run_controller(r);self.assertNotEqual(p.returncode,0)
        self.assertEqual(json.loads((r/'receipts/test.json').read_text())['exit_code'],7)
        self.assertEqual(json.loads((r/'STATUS.json').read_text())['status'],'FAILED_OR_STOPPED')
    def test_timeout_terminates_descendants(self):
        r=self.fixture('timeout');p=self.run_controller(r);self.assertNotEqual(p.returncode,0)
        self.assertEqual(json.loads((r/'receipts/test.json').read_text())['stop_reason'],'TIME_LIMIT')
        child=int((r/'runs/test/child.pid').read_text());stat=Path('/proc')/str(child)/'stat'
        if stat.exists():self.assertEqual(stat.read_text().rsplit(')',1)[1].split()[0],'Z')
    def test_source_drift_prevents_launch(self):
        r=self.fixture('complete');(r/'workflow/baseline_engine/run.py').write_text('raise SystemExit(0)')
        p=self.run_controller(r);self.assertNotEqual(p.returncode,0);self.assertFalse((r/'runs').exists())

if __name__=='__main__':unittest.main()
