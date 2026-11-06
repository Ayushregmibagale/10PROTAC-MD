#!/usr/bin/env python3
"""Apply the focused derivative-only anchor sampler repair; refuse unexpected source."""
from pathlib import Path
import argparse,json,hashlib,difflib
HELPER='''def distance_anchor_indices(full,heads,anchors):
    """Require exact disjoint stereochemical heads and unambiguous full anchors."""
    if len(heads)!=2 or len(anchors)!=2:
        raise ValidationError('Exactly two heads and anchors are required')
    if any(not 0<=i<h.GetNumAtoms() for h,i in zip(heads,anchors)):
        raise ValidationError('Anchor outside its exact heavy-atom head')
    pairs={(a[anchors[0]],b[anchors[1]]) for a,b in match_heads(full,heads)}
    if len(pairs)!=1:
        raise ValidationError('Ambiguous full-molecule anchor mapping')
    return next(iter(pairs))


'''

def repair(root):
    if root.name!='workflow' or not (root.parent/'SOURCE_MANIFEST.json').is_file():raise ValueError('Repair is restricted to the isolated workflow copy')
    base=root/'baseline_engine';c=base/'constrained_linker.py';s=c.read_text()
    if 'def distance_anchor_indices(' not in s:s=s.replace('def fit(mobile,reference):',HELPER+'def fit(mobile,reference):');c.write_text(s)
    for p in [base/'pyrosetta_adapter.py',base/'runtime/adapter/pyrosetta_adapter.py']:
        s=p.read_text();start=s.index('def sample_dist_free_conf(');end=s.index('def _make_docked_heads_mol',start);part=s[start:end]
        old_start=part.index('    # Map anchor atoms');old_end=part.index('    # Generate free conformers',old_start)
        part=part[:old_start]+'''    from constrained_linker import distance_anchor_indices
    anchor_A,anchor_B=distance_anchor_indices(mol,[HeadA,HeadB],Anchors)

'''+part[old_end:]
        part=part.replace('AllChem, rdFMCS, rdMolTransforms','AllChem, rdMolTransforms').replace('params.numThreads = 4','params.numThreads = 1')
        doc_start=part.index('    """');doc_end=part.index('    """',doc_start+7)
        part=part[:doc_start]+'''    """Inherited free ETKDGv3 distance histogram with exact stereo head mapping.

    This is a seed-dependent sampling heuristic, not an exhaustive physical
    linker-reach distribution or a probability estimate. Both-head constrained
    final linker construction is a separate stage.
'''+part[doc_end:]
        part=part.replace('    return (min(min(high_values), avg_index - min_margin),\n            max(max(high_values), avg_index + min_margin))','''    bounds=(min(min(high_values),avg_index-min_margin),
            max(max(high_values),avg_index+min_margin))
    from protac_qc import atomic_json,sha256
    atomic_json(Path('anchor_distance_sampling.json'),{
        'status':'SAMPLED_HEURISTIC_RANGE','requested_conformers':n,
        'generated_conformers':len(conf_ids),'num_threads':params.numThreads,
        'seed':params.randomSeed,'embedding_max_iterations_setting':params.maxIterations,
        'full_molecule_anchor_indices_zero_based':[anchor_A,anchor_B],
        'mapping':'Exact disjoint stereochemical complete heads; unique anchor pair',
        'histogram_bin_width_A':1.0,'histogram':histogram,'hist_threshold':hist_threshold,
        'min_margin_A':min_margin,'range_A':bounds,
        'source_sha256':{str(p):sha256(p) for p in [*Heads,Linkers]},
        'interpretation':'Finite free-conformer heuristic; no biological negative or population inferred'})
    return bounds''')
        if 'anchor_distance_sampling.json' not in part:raise ValueError('Cannot identify distance-range return')
        p.write_text(s[:start]+part+s[end:])
    (base/'runtime/adapter/constrained_linker.py').write_text(c.read_text())
    p=base/'runtime_manifest.json';m=json.loads(p.read_text());m['patched_runtime_sha256']={str(p.relative_to(base/'runtime')):hashlib.sha256(p.read_bytes()).hexdigest() for p in (base/'runtime').rglob('*') if p.is_file() and '__pycache__' not in p.parts};p.write_text(json.dumps(m,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workflow',type=Path,required=True);a=p.parse_args();repair(a.workflow)
