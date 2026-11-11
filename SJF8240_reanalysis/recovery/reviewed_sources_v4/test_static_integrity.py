from pathlib import Path
import importlib.util
import pytest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('candidate_builder',ROOT/'candidate_scripts_v4/01_build_static_graphs.py');builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)

def atom(name,chain,res=1,element='C'):
 return (f"ATOM  {res:5d} {name:^4} ALA {chain}{res:4d}    {1.:8.3f}{2.:8.3f}{3.:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2}  ").encode()

def fixture():
 data=b'\n'.join([atom('CA','A'),atom('CA','C'),atom('CA','D')]+[atom('C'+str(i),'X',i) for i in range(1,80)])+b'\n'
 cfg={'kin_chain':'D','kin_min':0};g=builder.parse(data,cfg)
 return data,cfg,{r:[x[0] for x in v] for r,v in g.items()}

def test_empty_is_excluded_explicitly():
 data,cfg,expected=fixture();g,rec,reason=builder.validated_pose(b'',cfg,expected)
 assert g is None and reason=='EMPTY_FILE'

def test_truncated_atomic_record_is_excluded():
 data,cfg,expected=fixture();g,rec,reason=builder.validated_pose(data+b'ATOM    20  CA  ALA ',cfg,expected)
 assert reason=='TRUNCATED_ATOM_RECORD'

def test_missing_nodes_are_excluded_for_both_role_graphs():
 data,cfg,expected=fixture();partial=b'\n'.join(data.splitlines()[1:])+b'\n'
 g,rec,reason=builder.validated_pose(partial,cfg,expected)
 assert reason=='INCOMPLETE_REFERENCE_NODE_SET'

def test_unexpected_node_identity_remains_fatal():
 data,cfg,expected=fixture();altered=data.replace(b' CA  ALA D   1',b' CA  ALA D   2')
 assert altered!=data
 with pytest.raises(ValueError,match='identity change'):builder.validated_pose(altered,cfg,expected)

def test_complete_pose_is_retained():
 data,cfg,expected=fixture();g,rec,reason=builder.validated_pose(data,cfg,expected)
 assert reason is None and len(g['ligand'])==79
