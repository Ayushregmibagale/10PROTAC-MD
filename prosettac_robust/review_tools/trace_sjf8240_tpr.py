#!/usr/bin/env python3
"""Inspect compiled production TPR graphs read-only; never load trajectories."""
import argparse,json,hashlib,warnings
from pathlib import Path
from collections import Counter
import MDAnalysis as mda
from trace_sjf8240_lineage import itp_graph,sha

def run(root,review):
    if review.resolve().is_relative_to(root.resolve()):raise ValueError('Audit output must be outside original project')
    ref=root/'ligand_params/SJF8240_3d.acpype/SJF8240_3d_GMX.itp';expected=itp_graph(ref)[0]
    metadata=json.loads((review/'md_graph_pose_lineage.json').read_text());rows=[]
    for i,row in enumerate(metadata):
        p=root/Path(row['parameter_path']).parent.parent/'prod/prod.tpr';result={'poi':row['poi'],'pose_id':row['pose_id'],'tpr':str(p),'exists':p.exists()}
        if p.exists():
            try:
                u=mda.Universe(str(p));a=u.select_atoms('resname UNL and mass 2.1:1000');idx=set(map(int,a.indices));names={int(x.index):str(x.name) for x in a}
                if len(a)!=79 or len(a.residues)!=1 or len(set(names.values()))!=79:raise ValueError('Ambiguous ligand selection')
                bonds=sorted(tuple(sorted((names[int(b.atoms[0].index)],names[int(b.atoms[1].index)]))) for b in a.bonds if all(int(x.index) in idx for x in b.atoms))
                key={'atoms':sorted((str(x.name),round(float(x.mass),2)) for x in a),'bonds':bonds}
                sig=hashlib.sha256(json.dumps(key,sort_keys=True).encode()).hexdigest();result.update(sha256=sha(p),heavy_atoms=len(a),heavy_connectivity_signature=sig,matches_legacy_connectivity=sig==expected)
            except Exception as e:result['error']=repr(e)
        rows.append(result)
        if (i+1)%20==0:print('TPR checked',i+1,'/',len(metadata),flush=True)
    result={'scope':'ACTUAL_COMPILED_PRODUCTION_TPR_HEAVY_CONNECTIVITY','rows':rows,'total':len(rows),'matches_legacy_connectivity':sum(x.get('matches_legacy_connectivity',False) for x in rows),'missing_or_unresolved':[x for x in rows if not x.get('matches_legacy_connectivity',False)],'by_poi':dict(Counter(x['poi'] for x in rows if x.get('matches_legacy_connectivity',False))),'limitations':['TPR graph proves compiled ligand connectivity; trajectory/TPR identity and original model execution require separate provenance.','Named heavy adjacency and masses do not encode bond order or stereochemistry.']}
    (review/'production_tpr_chemistry.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k not in ['rows','limitations','missing_or_unresolved']},flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--review',type=Path,required=True);a=p.parse_args();run(a.root,a.review)
