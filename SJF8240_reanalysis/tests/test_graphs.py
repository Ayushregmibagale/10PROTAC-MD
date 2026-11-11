from pathlib import Path
import sys,importlib.util
import numpy as np,pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts/lib'))
spec=importlib.util.spec_from_file_location('builder',ROOT/'scripts/01_build_static_graphs.py');builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
def group(chain,name,x):return [((chain,1,'',name,'C','ALA'),[x,0.,0.])]
def test_true_vhl_selection_changes_contact_without_changing_target_or_ligand():
 g={'elob':group('A','CA',50.),'vhl':group('C','CA',1.),'poi':group('D','CA',4.),'ligand':group('X','C1',0.)}
 old=builder.Accumulator(g,'elob');new=builder.Accumulator(g,'vhl');old.add(g);new.add(g)
 assert old.count[0,2]==0 and new.count[0,2]==1
 assert old.count[1,2]==new.count[1,2]==1
 assert old.keys[1:]==new.keys[1:]
def test_static_graph_rejects_changed_residue_correspondence():
 g={'elob':group('A','CA',50.),'vhl':group('C','CA',1.),'poi':group('D','CA',4.),'ligand':group('X','C1',0.)};acc=builder.Accumulator(g,'vhl')
 g['poi']=group('D','CB',4.)
 with pytest.raises(ValueError,match='correspondence'):acc.add(g)

def test_historical_coordinate_overflow_is_exact_and_explicit():
 line=b'HETATM    1  C1  PTC X   1     420.195-1083.518-174.581  1.00  0.00           C  '
 xyz,recovered=builder.coordinates(line)
 assert recovered and xyz==[420.195,-1083.518,-174.581]

def test_ambiguous_coordinate_record_is_rejected():
 line=b'HETATM    1  C1  PTC X   1     420.195-1083.518-174.581 garbage'
 with pytest.raises(ValueError,match='Ambiguous'):builder.coordinates(line)

def test_standard_coordinate_record_is_unchanged():
 line=b'HETATM    1  C1  PTC X   1      42.195 -83.518 -74.581  1.00  0.00           C  '
 xyz,recovered=builder.coordinates(line)
 assert not recovered and xyz==[42.195,-83.518,-74.581]


def test_valid_fixed_width_two_decimal_large_coordinate_preserved():
 line=b'HETATM    1 C1   PTC X   1     418.453-1083.32   5.221  1.00  0.00           C  '
 xyz,recovered=builder.coordinates(line)
 assert not recovered and xyz==[418.453,-1083.32,5.221]
