"""Orchestration contracts only; these do not validate original-engine science."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from campaign import outputs, validate_tasks
from check_paths import sha


def test_changed_input_is_rejected(tmp_path):
    d = tmp_path/'inputs'; d.mkdir()
    p = d/'full.smi'; p.write_text('CC\n')
    tasks = {'tasks':[dict(id='drug__target',input='inputs',files={'full.smi':sha(p)})]}
    validate_tasks(tmp_path, tasks)
    p.write_text('CCC\n')
    with pytest.raises(ValueError, match='Input changed'):
        validate_tasks(tmp_path,tasks)


def test_missing_outputs_are_not_success(tmp_path):
    with pytest.raises(RuntimeError,match='without score.sc'):
        outputs(tmp_path)


def test_nonfinite_scores_are_rejected(tmp_path):
    (tmp_path/'Results').mkdir()
    (tmp_path/'Results/model.pdb').write_text('END\n')
    (tmp_path/'Patchdock_Results').mkdir()
    (tmp_path/'Patchdock_Results/score.sc').write_text('SCORE: nan model\n')
    with pytest.raises(RuntimeError,match='Nonfinite'):
        outputs(tmp_path)


def fake_original(tmp_path, command):
    source = tmp_path/'original'; source.mkdir()
    (source/'cluster/LOCAL').mkdir(parents=True)
    (source/'cluster/__init__.py').write_text('def getCluster(name): return None\n')
    (source/'cluster/LOCAL/__init__.py').write_text('')
    (source/'cluster/LOCAL/LOCAL.py').write_text('class LOCAL:\n def __init__(self,n_workers=None): self.n_workers=n_workers\n')
    (source/'main.py').write_text('import cluster\nfrom pathlib import Path\ndef main(name,args):\n Path("ORIGINAL_CALLED").write_text(name)\n cluster.getCluster("LOCAL").runBatchCommands(['+repr(command)+'],mem=1)\n')
    cfg=tmp_path/'config.json'
    cfg.write_text(json.dumps(dict(source=str(source),patchdock='/missing',rosetta='/missing',openbabel='/missing')))
    return subprocess.run([sys.executable,str(HERE/'run_original.py'),'--config',str(cfg),'--workers','2'],cwd=tmp_path,capture_output=True,text=True)


def test_original_entrypoint_called_without_replacement_engine(tmp_path):
    result=fake_original(tmp_path,'true')
    assert result.returncode==0,result.stderr
    assert (tmp_path/'ORIGINAL_CALLED').read_text()=='main.py'


def test_original_worker_failure_propagates(tmp_path):
    result=fake_original(tmp_path,'false')
    assert result.returncode!=0
    assert 'Original subprocess failed' in result.stderr
