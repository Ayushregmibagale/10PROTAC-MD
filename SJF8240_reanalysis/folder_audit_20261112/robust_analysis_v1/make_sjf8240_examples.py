"""Explicit migration adapter for existing SJF8240 graphs, not a ten-compound default."""
from pathlib import Path
import argparse,json
from contracts import *
def main(root,out):
 root=Path(root).resolve();out=Path(out);out.mkdir(parents=True,exist_ok=True);labelpath=root/'config/current_md_labels.json';labels=json.loads(labelpath.read_text())['labels'];evidence=root.parent/'prosettac_robust/existing_result_audit/protein_roles_core_complete/per_pose_roles.json'
 def binding(path):return {'path':str(path),'sha256':sha(path)}
 for model in ['original_md','lysine_temporal']:
  c={'schema_version':1,'compound_id':'SJF8240','units':{'coordinates':'angstrom','time':'ns'},'model':model,'evaluation':'leave_pair_out','training':{'epochs':300,'seed':42,'lr':.001,'weight_decay':.0001,'loss':'unweighted','threshold':.5},'roles':{r:{'uniprot':uid,'identity_evidence':binding(evidence)} for r,uid in [('vhl','P40337'),('elob','Q15370')]},'targets':[],'purpose':'Reusable protocol example; no new scientific training has been performed by this adapter.'}
  names=list(labels)+(['SIK3_full','RPS6KB2_full'] if model=='original_md' else [])
  for target in names:
   row={'target_id':target,'protein_group':target,'label':labels.get(target),'label_status':'confirmed_user' if target in labels else 'unknown','label_evidence':binding(labelpath) if target in labels else None,'graphs':{}}
   for role in c['roles']:
    path=root/'corrected_md_only_v1/inputs'/role/(target+'.npz') if model=='original_md' else root/'mechanistic_v3/inputs/temporal'/role/(target+'.npz')
    row['graphs'][role]={**binding(path),'compound_id':'SJF8240','target_id':target,'role':role}
   c['targets'].append(row)
  path=out/(model+'.json');atomic_json(path,c);validate_study(path)
 print('Validated two explicit SJF8240 example manifests; other compounds need their own labels, graphs and role evidence.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.project_root,a.output)
