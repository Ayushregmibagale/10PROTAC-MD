"""Receipt-bound solute trajectory loading for the corrected SJF8240 graph schema."""
import json
import numpy as np
from protac_qc import ValidationError,sha256

def mapped_heavy_groups(universe,meta):
    """Validate executed one-based heavy-atom indices without silent coercion."""
    solute=meta.get('solute_atoms')
    if type(solute) is not int or solute<=0 or len(universe.atoms)<solute:
        raise ValidationError('Incomplete or invalid mapped solute topology')
    groups={};used=set()
    for role in ['vhl','poi','protac']:
        raw=np.asarray(meta.get('atom_indices_one_based',{}).get(role,[]))
        if raw.ndim!=1 or not len(raw) or not np.issubdtype(raw.dtype,np.integer):
            raise ValidationError('Mapped heavy groups require nonempty integer indices')
        if raw.min()<1 or raw.max()>solute or len(set(raw))!=len(raw) or used.intersection(raw):
            raise ValidationError('Invalid, duplicate or overlapping mapped heavy groups')
        used.update(raw);groups[role]=universe.atoms[raw-1]
        masses=np.asarray(groups[role].masses)
        if not np.isfinite(masses).all() or np.any(masses<1.5):
            raise ValidationError('Nonfinite mass or hydrogen in mapped heavy group')
    return groups

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
    from read_only_trajectory import ReadOnlyXTCReader
    if trajectory.suffix.lower()!='.xtc':raise ValidationError('Qualified read-only reader requires an XTC trajectory')
    u=mda.Universe(str(topology),str(trajectory),format=ReadOnlyXTCReader)
    groups=mapped_heavy_groups(u,meta)
    return u,groups
