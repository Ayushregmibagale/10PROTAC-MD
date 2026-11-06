#!/usr/bin/env python3
import os,sys,json,hashlib,importlib.util,argparse
from pathlib import Path
r=Path(__file__).resolve().parents[1];w=r/'workflow';sys.path[:0]=[str(w/'scripts/lib'),str(w/'baseline_engine')]
from protac_qc import read_params,sha256
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
folder=w/'prepared_inputs/SJF8240/MET_P08581';heads=[str(folder/n) for n in ['vhl_head.sdf','poi_head.sdf']];anchors=[int(x)-1 for x in read_params(folder/'Protac_params.txt')['Anchor atoms'].split()];results={}
for label,path in [('historical_patched_baseline',w/'validation_evidence/execution_sources/612dc20ac9706f93/baseline_engine/pyrosetta_adapter.py'),('isolated_hardened_sampler',w/'baseline_engine/pyrosetta_adapter.py')]:
 d=out/label;d.mkdir();os.chdir(d);os.environ.update(PROTAC_ATTEMPT_ROOT=str(d),PROTAC_SEED='155921',PROTAC_TASK_KEY='parent')
 s=importlib.util.spec_from_file_location(label,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
 bounds=m.sample_dist_free_conf(heads,anchors,str(folder/'full.smi'),n=200)
 results[label]={'source_sha256':sha256(path),'range_A':list(bounds),'histogram_sha256':sha256(d/'initial_distances.hist'),'conformer_sdf_sha256':sha256(d/'random_sampling.sdf')}
result={'status':'PAIRED_ACTUAL_DOCKING_ENVIRONMENT_DISTANCE_CONTROL','seed':155921,'requested_conformers':200,'results':results,'same_histogram':len({x['histogram_sha256'] for x in results.values()})==1,'same_range':results['historical_patched_baseline']['range_A']==results['isolated_hardened_sampler']['range_A'],'scope':'Initial distance sampling only; does not establish docking recovery, MD convergence or GNN predictive validity'}
(out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
