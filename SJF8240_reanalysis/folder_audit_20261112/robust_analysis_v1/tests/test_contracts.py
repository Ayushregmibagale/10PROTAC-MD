from pathlib import Path
import sys,json,copy,subprocess,os
import numpy as np,pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from contracts import *
@pytest.fixture
def study(tmp_path):
 evidence=tmp_path/'evidence.json';labels={'p1':1,'p2':1,'n1':0,'n2':0};atomic_json(evidence,{'labels':labels});binding={'path':str(evidence),'sha256':sha(evidence)};targets=[]
 for target,label in labels.items():
  graph=tmp_path/(target+'.npz');np.savez(graph,x=np.array([[1,0,0],[0,1,1]],np.float32),edge_index=np.array([[0,1],[1,0]]),edge_attr=np.array([[.5,1,1],[.5,1,1]],np.float32),poi=np.array([target]),label=np.array([label]))
  targets.append({'target_id':target,'protein_group':target,'label':label,'label_status':'confirmed_user','label_evidence':binding,'graphs':{'vhl':{'path':str(graph),'sha256':sha(graph),'role':'vhl','target_id':target,'compound_id':'test'}}})
 c={'schema_version':1,'compound_id':'test','units':{'coordinates':'angstrom','time':'ns'},'model':'original_md','evaluation':'leave_pair_out','training':{'epochs':1,'seed':42,'lr':.001,'weight_decay':.0001,'loss':'balanced','threshold':.5},'roles':{'vhl':{'uniprot':'P40337','identity_evidence':binding}},'targets':targets};path=tmp_path/'study.json';atomic_json(path,c);return path,c

def test_valid_study(study):
 p,c=study;_,labels,paths=validate_study(p);assert len(labels)==4 and len(paths)==4
@pytest.mark.parametrize('change',['unknown_label','duplicate_alias','wrong_role','missing_graph','bad_label','wrong_compound','missing_evidence','wrong_units','missing_arm','bad_epoch'])
def test_reject_bad_study(study,change):
 p,c=study;r=c['targets'][0]
 if change=='unknown_label':r['label_status']='unknown'
 elif change=='duplicate_alias':r['protein_group']=c['targets'][1]['protein_group']
 elif change=='wrong_role':c['roles']['vhl']['uniprot']='Q15370'
 elif change=='missing_graph':r['graphs']['vhl']['path']='absent.npz'
 elif change=='bad_label':r['label']=False
 elif change=='wrong_compound':r['graphs']['vhl']['compound_id']='other'
 elif change=='missing_evidence':r['label_evidence']['sha256']='0'*64
 elif change=='wrong_units':c['units']['coordinates']='nm'
 elif change=='missing_arm':r['graphs']={}
 else:c['training']['epochs']=0
 atomic_json(p,c)
 with pytest.raises(ContractError):validate_study(p)

def test_hash_tamper(study):
 p,c=study;graph=Path(c['targets'][0]['graphs']['vhl']['path']);graph.write_bytes(graph.read_bytes()+b'x')
 with pytest.raises(ContractError,match='hash'):validate_study(p)

@pytest.mark.parametrize('kind',['float_index','oob','nan','bad_frequency','bad_role'])
def test_invalid_graph(tmp_path,kind):
 x=np.array([[1,0,0],[0,1,1]],float);ei=np.array([[0,1],[1,0]]);ea=np.array([[.5,1,1],[.5,1,1]])
 if kind=='float_index':ei=ei.astype(float)
 elif kind=='oob':ei[0,0]=2
 elif kind=='nan':x[0,2]=np.nan
 elif kind=='bad_frequency':ea[0,0]=1.1
 else:x[0,0]=0
 p=tmp_path/'bad.npz';np.savez(p,x=x,edge_index=ei,edge_attr=ea)
 with pytest.raises(ContractError):graph_arrays(p,3)

def test_empty_edges_are_valid(tmp_path):
 p=tmp_path/'empty.npz';np.savez(p,x=np.array([[1,0,0],[0,1,1]],float),edge_index=np.empty((2,0),int),edge_attr=np.empty((0,3)));assert graph_arrays(p,3)[0][1].shape==(2,0)

def test_pair_exclusion_and_train_only_weights():
 labels={'p1':1,'p2':1,'p3':1,'n1':0,'n2':0};fs=folds(labels,'leave_pair_out');assert len(fs)==6
 for f in fs:
  assert not set(f['held_out'])&set(f['training']);weights=balanced_weights(labels,f['training']);tot=[sum(weights[t] for t in f['training'] if labels[t]==v) for v in [0,1]];assert tot[0]==tot[1]

def test_pair_metric_ties_and_threshold():
 r=rank_metrics([{'positive_score':.5,'negative_score':.5},{'positive_score':.6,'negative_score':.4}]);assert r['matched_pair_AUROC']==.75 and r['specificity']==.5

def test_sequence_swap_is_rejected():
 refs={'P40337':'MATRVLESGPK','Q15370':'CWYFHDNQICA'};assert verify_sequence_identity(refs['P40337'],'P40337',refs)==1
 with pytest.raises(ContractError):verify_sequence_identity(refs['Q15370'],'P40337',refs)

def test_atomic_json_rejects_nan(tmp_path):
 p=tmp_path/'record.json';atomic_json(p,{'score':1})
 with pytest.raises(ValueError):atomic_json(p,{'score':float('nan')})
 assert json.loads(p.read_text())=={'score':1}

def test_fresh_training_resume_cpu_and_source_drift(study,tmp_path):
 p,c=study;output=tmp_path/'run';cmd=[sys.executable,str(ROOT/'run_study.py'),'--study',str(p),'--output',str(output),'--role','vhl'];env=os.environ.copy();env.pop('PYTHONPATH',None)
 for suffix in [[],[],['--verify']]:
  result=subprocess.run(cmd+suffix,env=env,text=True,capture_output=True);assert result.returncode==0,result.stderr
 assert json.loads((output/'cpu_verification.json').read_text())['folds']==4
 hashes={x.name:sha(x) for x in output.glob('*.pt')};assert len(hashes)==4
 c['training']['seed']=43;atomic_json(p,c);result=subprocess.run(cmd,env=env,text=True,capture_output=True);assert result.returncode!=0 and 'intent drift' in result.stderr
 assert hashes=={x.name:sha(x) for x in output.glob('*.pt')}

def test_temporal_window_correspondence(tmp_path):
 arrays={}
 for w in range(4):
  x=np.zeros((2,29),np.float32);x[0,0]=1;x[1,1]=1;x[1,2]=1;x[0,3]=1;x[1,11]=1;x[1,23]=1;x[1,24:]=[.2,.1,.8,.4,.1];arrays.update({f'x{w}':x,f'edge_index{w}':np.array([[0,1],[1,0]]),f'edge_attr{w}':np.array([[.5,1,1],[.5,1,1]])})
 p=tmp_path/'temporal.npz';np.savez(p,**arrays);assert len(graph_arrays(p,29,4))==4
 arrays['x3']=arrays['x3'][::-1];np.savez(p,**arrays)
 with pytest.raises(ContractError,match='identities'):graph_arrays(p,29,4)

def test_temporal_training_and_fresh_cpu(study,tmp_path):
 p,c=study;c['model']='lysine_temporal'
 for row in c['targets']:
  arrays={}
  for w in range(4):
   x=np.zeros((2,29),np.float32);x[0,0]=1;x[1,1]=1;x[0,3]=1;x[1,11]=1;x[1,23]=1;x[1,24:]=[.2,.1,.8,.4,.1];arrays.update({f'x{w}':x,f'edge_index{w}':np.array([[0,1],[1,0]]),f'edge_attr{w}':np.array([[.5,1,1],[.5,1,1]])})
  path=Path(row['graphs']['vhl']['path']);np.savez(path,**arrays);row['graphs']['vhl']['sha256']=sha(path)
 atomic_json(p,c);cmd=[sys.executable,str(ROOT/'run_study.py'),'--study',str(p),'--output',str(tmp_path/'temporal'),'--role','vhl'];env=os.environ.copy();env.pop('PYTHONPATH',None)
 for extra in [[],['--verify']]:
  result=subprocess.run(cmd+extra,env=env,capture_output=True,text=True);assert result.returncode==0,result.stderr

def test_unknown_is_prediction_only(study):
 p,c=study;row=copy.deepcopy(c['targets'][0]);row.update(target_id='unknown',protein_group='unknown',label=None,label_status='unknown',label_evidence=None);g=row['graphs']['vhl'];g['target_id']='unknown';path=p.parent/'unknown.npz'
 with np.load(g['path']) as z:arrays={k:v.copy() for k,v in z.items()}
 arrays['poi']=np.array(['unknown']);arrays['label']=np.array([-1]);np.savez(path,**arrays);g.update(path=str(path),sha256=sha(path));c['targets'].append(row);atomic_json(p,c);_,labels,_=validate_study(p);assert 'unknown' not in labels and all('unknown' not in f['training'] for f in folds(labels,'leave_pair_out'))

def test_anticorrelated_contacts_do_not_create_bridge_edge():
 # Ligand contacts only one protein in each frame: marginal product would be 0.25.
 joint=(joint_contact_observation([1],[10],5)+joint_contact_observation([10],[1],5))/2
 assert float(joint[0,0])==0

def test_periodic_distance_translation():
 from MDAnalysis.lib.distances import distance_array
 box=np.array([30,30,30,90,90,90]);a=np.array([[1.,2,3]]);b=np.array([[29.,2,3]])
 d=distance_array(a,b,box=box);shifted=distance_array(a+np.array([30,0,0]),b-np.array([0,30,0]),box=box);assert np.allclose(d,shifted,atol=1e-6) and abs(d[0,0]-2)<1e-5
