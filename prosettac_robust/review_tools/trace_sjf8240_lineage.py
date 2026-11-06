#!/usr/bin/env python3
"""Read original SJF8240 artifacts; write audit files only under --output.

Bond signatures establish named heavy-atom connectivity, not bond order/stereo.
No trajectories, original tables, graphs or models are modified or executed.
"""
import argparse,csv,hashlib,json
from pathlib import Path
from collections import Counter

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(p):
    if not p.exists():return []
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def itp_graph(p):
    section='';atoms={};bonds=[]
    for line in p.read_text().splitlines():
        line=line.split(';',1)[0].strip()
        if not line:continue
        if line.startswith('['):section=line.strip('[] ').strip();continue
        if line.startswith('#'):continue
        x=line.split()
        if section=='atoms':
            if len(x)<8:raise ValueError('Incomplete atom row: '+str(p))
            atoms[int(x[0])]={'name':x[4],'mass':float(x[7]),'type':x[1],'charge':float(x[6])}
        elif section=='bonds':bonds.append((int(x[0]),int(x[1])))
    heavy={i:a for i,a in atoms.items() if a['mass']>2.1}
    if not atoms or not bonds:raise ValueError('Missing topology atoms/bonds: '+str(p))
    if len({a['name'] for a in heavy.values()})!=len(heavy):raise ValueError('Nonunique heavy-atom names')
    key={'atoms':sorted((a['name'],round(a['mass'],2)) for a in heavy.values()),
         'bonds':sorted(tuple(sorted((heavy[i]['name'],heavy[j]['name']))) for i,j in bonds if i in heavy and j in heavy)}
    return hashlib.sha256(json.dumps(key,sort_keys=True).encode()).hexdigest(),len(heavy),len(atoms)
def audit(root,out):
    root=root.resolve();out=out.resolve()
    if out==root or out.is_relative_to(root):raise ValueError('Output must be outside original project')
    out.mkdir(parents=True,exist_ok=True)
    source=root/'ligand_params/SJF8240_3d.acpype/SJF8240_3d_GMX.itp'
    signature,nheavy,natoms=itp_graph(source);systems=[]
    for base in ['simulations','random_screen','kinome_md_screen']:
        for p in sorted((root/base).rglob('SJF8240.itp')):
            sig,nh,na=itp_graph(p);top=p.parent/'topol.top';text=top.read_text() if top.exists() else ''
            systems.append({'path':str(p.relative_to(root)),'sha256':sha(p),'heavy_connectivity_signature':sig,
                            'matches_legacy_parameter_connectivity':sig==signature,'heavy_atoms':nh,'all_atoms':na,
                            'topol_top_exists':top.exists(),'topol_top_sha256':sha(top) if top.exists() else None,
                            'topol_includes_sjf8240_itp':'#include "SJF8240.itp"' in text,
                            'directed_include_molecule_present':any(s.split(';')[0].split()[:1]==['SJF8240_3d'] for s in text.splitlines())})
    bypath={x['path']:x for x in systems};pose_rows=table(root/'analysis/gnn/md_contact_per_pose.tsv');linked=[]
    for x in pose_rows:
        rel=f"random_screen/{x['poi']}/{x['pose_id']}/topol/SJF8240.itp";hit=bypath.get(rel)
        linked.append(dict(x,parameter_path=rel,parameter_found=bool(hit),legacy_connectivity=hit['matches_legacy_parameter_connectivity'] if hit else None))
    evidence={}
    for rel in ['analysis/gnn/ensemble_graphs.pt','analysis/gnn/md_contact_graphs.pt','analysis/features/ml_table.tsv',
                'analysis/gnn/e2gnn_lopo_results.tsv','analysis/gnn/e2gnn_lopo_summary.tsv',
                'scripts/02_analysis/04a_build_ensemble_graphs.py','scripts/02_analysis/04c_build_md_contact_graphs.py',
                'scripts/02_analysis/04g_e2_scalar_gnn.py','random_screen/manifest.tsv']:
        p=root/rel
        if p.exists():evidence[rel]={'sha256':sha(p),'bytes':p.stat().st_size}
    summary={'scope':'READ_ONLY_HISTORICAL_SJF8240_CONNECTIVITY_AND_GNN_LINEAGE_AUDIT',
             'original_root':str(root),'legacy_parameter_source_sha256':sha(source),'legacy_heavy_connectivity_signature':signature,
             'topologies_examined':len(systems),'matching_legacy_connectivity':sum(x['matches_legacy_parameter_connectivity'] for x in systems),
             'connectivity_mismatches':[x['path'] for x in systems if not x['matches_legacy_parameter_connectivity']],
             'md_graph_pose_metadata_rows':len(linked),'md_graph_rows_with_matching_legacy_parameters':sum(x['legacy_connectivity'] is True for x in linked),
             'unresolved_md_graph_rows':[x for x in linked if not x['parameter_found']],
             'by_poi_md_graph_parameters':dict(Counter(x['poi'] for x in linked if x['legacy_connectivity'] is True)),
             'artifact_sha256':evidence,
             'limitations':['Text topologies and graph-builder metadata establish retained artifact lineage, not a cryptographic original-run receipt.',
                            'ITP signature checks named heavy-atom adjacency and masses, not stereochemistry/bond orders or actual TPR contents.',
                            'Static graph and GNN table rows need separate review; MD provenance does not certify their independence or accuracy.']}
    for name,value in [('lineage_summary.json',summary),('topology_connectivity_inventory.json',systems),('md_graph_pose_lineage.json',linked)]:
        (out/name).write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['artifact_sha256','limitations','unresolved_md_graph_rows']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();audit(a.root,a.output)
