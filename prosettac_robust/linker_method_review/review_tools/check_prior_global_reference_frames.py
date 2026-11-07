"""Exploratory regression on all previously accepted GSK215 global frames.

Selected by historical QC acceptance, not by current candidate yield. This is
not an independent benchmark; both successful and non-native frames are kept.
"""
from pathlib import Path
import sys,json,importlib.util,hashlib
import numpy as np
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'workflow/baseline_engine';sys.path.insert(0,str(BASE))
from official_parameters import docked_heads
raw=Path('/tmp/10protac_finish/benchmarks/global_recovery_gsk215_seed155921_v27/engine_raw')
summary=json.loads((ROOT/'workflow/validation_evidence/global_recovery_gsk215_v27/recovery_summary.json').read_text())
headpaths=[raw/'vhl_head_H.sdf',raw/'poi_head_H.sdf']
heads=[Chem.SDMolSupplier(str(p),removeHs=True)[0] for p in headpaths]
full=Chem.MolFromSmiles((raw/'full.smi').read_text().split()[0]);n=heads[0].GetNumAtoms();rows=[]
for variant in ['candidate_workflow','candidate_torsion_workflow']:
 path=ROOT/'linker_method_review'/variant/'baseline_engine/constrained_linker.py'
 spec=importlib.util.spec_from_file_location(variant+'_regression',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 for item in summary['poses']:
  name=item['pose'];_,glob,local,_=name.split('_');p=raw/'Patchdock_Results'/('docked_'+glob+'_'+local+'.sdf')
  xyz=docked_heads(p,headpaths,root=raw).GetConformer().GetPositions()
  for seed in [392,104729,155921]:
   mol,report=module.constrained_linker(full,heads,[xyz[:n],xyz[n:]],seed,attempts=3)
   rows.append({'variant':variant,'historical_pose':name,'historical_vhl_RMSD_A':item['vhl_CA_rmsd_after_POI_fit_A'],'historical_native_contact_fraction':item['fraction_native_contacts'],'seed':seed,'input_sdf_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'accepted':mol is not None,'details':report})
result={'status':'HISTORICALLY_SELECTED_GLOBAL_FRAME_REGRESSION_NOT_INDEPENDENT_VALIDATION','selection':'Every QC-accepted pose in the retained GSK215 global seed155921 run, including non-native arrangements','rows':rows,'production_release':False}
(ROOT/'linker_method_review/evidence/prior_global_reference_frame_regression.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'cases':len(rows),'accepted':sum(x['accepted'] for x in rows),'by_variant':{v:sum(x['accepted'] for x in rows if x['variant']==v) for v in ['candidate_workflow','candidate_torsion_workflow']}},indent=2))
