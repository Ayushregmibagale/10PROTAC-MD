"""Receipt-bound solute trajectory loading for the corrected SJF8240 graph schema."""
import json
import numpy as np
from protac_qc import ValidationError,sha256

def load_mapped_production(preparation,topology,trajectory,receipt):
    meta=json.loads(preparation.read_text());r=json.loads(receipt.read_text())
    if r.get('status')!='COMPLETE' or r.get('intent',{}).get('stage')!='prod':
        raise ValidationError('Completed production MD receipt required for GNN features')
    if r['intent']['inputs'].get('md_preparation.json')!=sha256(preparation):raise ValidationError('MD atom mapping differs from executed system')
    if r['outputs'].get(trajectory.name)!=sha256(trajectory):raise ValidationError('Changed or unbound production trajectory')
    if r['outputs'].get(topology.name)!=sha256(topology):
        p=json.loads(topology.with_suffix('.provenance.json').read_text())
        if p.get('group')!='SOLUTE' or p.get('target_sha256')!=sha256(topology) or p.get('source_sha256') not in [h for n,h in r['outputs'].items() if n.endswith('.tpr')]:
            raise ValidationError('Reduced topology is not bound to production receipt')
    import MDAnalysis as mda
    u=mda.Universe(str(topology),str(trajectory));groups={};used=set()
    if len(u.atoms)<meta['solute_atoms']:raise ValidationError('Incomplete mapped solute topology')
    for role in ['vhl','poi','protac']:
        ids=np.asarray(meta['atom_indices_one_based'][role],dtype=int)-1
        if not len(ids) or ids.min()<0 or ids.max()>=len(u.atoms) or len(set(ids))!=len(ids) or used.intersection(ids):
            raise ValidationError('Invalid, duplicate or overlapping mapped heavy groups')
        used.update(ids);groups[role]=u.atoms[ids]
        if np.any(groups[role].masses<1.5):raise ValidationError('Hydrogen in mapped heavy group')
    return u,groups
