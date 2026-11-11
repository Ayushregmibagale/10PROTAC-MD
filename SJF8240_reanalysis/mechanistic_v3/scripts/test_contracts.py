from pathlib import Path
import sys,json,itertools
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'));sys.path.insert(0,str(ROOT/'scripts'))
from common import *

def test_pairs_exclude_both_and_equal_prior():
 p=config();pairs=list(itertools.product([x for x in p['cohort'] if p['labels'][x]],['ABL1','AXL']));assert len(pairs)==12
 for pos,neg in pairs:
  train=[x for x in p['cohort'] if x not in [pos,neg]];assert len(train)==6 and sum(p['labels'][x] for x in train)==5;assert not set(train)&{pos,neg}

def test_baseline_unchanged():
 source=ROOT.parent/'SJF8240_corrected_md_only_v1';assert (ROOT/'lib/model.py').read_bytes()==(source/'lib/model.py').read_bytes()
 for role in ['elob','vhl']:
  _,g=load_inputs(role)
  assert len(g)==10

def test_sasa_occluder_and_translation():
 import mdtraj as md
 top=md.Topology();c=top.add_chain();r=top.add_residue('LYS',c);top.add_atom('NZ',md.element.nitrogen,r);top.add_atom('C',md.element.carbon,r)
 xyz=np.array([[[0,0,0],[.25,0,0]]],dtype=np.float32);buried=md.shrake_rupley(md.Trajectory(xyz,top),atom_indices=[0])[0,0];alone=md.shrake_rupley(md.Trajectory(xyz[:,:1],top.subset([0])),atom_indices=[0])[0,0];translated=md.shrake_rupley(md.Trajectory(xyz+1,top),atom_indices=[0])[0,0]
 assert 0<buried<alone and abs(translated-buried)<1e-6

def test_temporal_graphs_and_predictions():
 if not (ROOT/'config/temporal_input_manifest.json').exists():
  import pytest;pytest.skip('Temporal graphs not yet built')
 from temporal_train import load_temporal
 from temporal_model import TemporalEncoder
 p,g=load_temporal('vhl','cpu');model=TemporalEncoder().eval()
 for target,d in g.items():
  assert d.x.shape[1]==29 and d.num_graphs==4
  assert torch.isfinite(d.x).all() and torch.isfinite(model(d)).all()
  assert all(bool((d.x[:,3:23].sum(1)==1).all()) for _ in [0])
  assert bool((d.x[d.x[:,24]>0,23]==1).all())
