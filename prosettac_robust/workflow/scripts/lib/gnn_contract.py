"""Publication derivative data contracts; no historical data is silently upgraded."""
import json,hashlib
from pathlib import Path
import numpy as np
from protac_qc import ValidationError,sha256

def verify_gnn_release(path,inputs,verified_inchikey):
    release=json.loads(Path(path).read_text())
    if release.get('status')!='RELEASED_FOR_EXPLORATORY_GNN':raise ValidationError('GNN input release is absent')
    if release.get('compound')!='SJF8240' or release.get('ligand_inchikey')!=verified_inchikey:
        raise ValidationError('GNN data chemistry differs from verified SJF8240')
    if release.get('feature_schema')!='SJF8240-periodic-contacts-v1':raise ValidationError('GNN feature schema not reviewed')
    for key,p in inputs.items():
        if release.get('input_sha256',{}).get(key)!=sha256(p):raise ValidationError('Changed or unreviewed GNN input: '+key)
    labels=release.get('experimental_labels',{})
    if len(labels)<2 or any(type(x) is not int or x not in (0,1) for x in labels.values()) or len(set(labels.values()))!=2:
        raise ValidationError('Explicit reviewed experimental binary labels required')
    groups=release.get('protein_groups',{})
    if set(groups)!=set(labels) or any(not isinstance(v,str) or not v for v in groups.values()):
        raise ValidationError('Reviewed protein identity groups required for every target')
    if len(set(groups.values()))!=len(groups):raise ValidationError('A protein target appears multiple times in LOPO cohort')
    provenance=release.get('label_provenance_file')
    if not provenance:raise ValidationError('Experimental label provenance is required')
    provenance=(Path(path).resolve().parent/provenance).resolve()
    if not provenance.is_file() or sha256(provenance)!=release.get('label_provenance_sha256'):
        raise ValidationError('Experimental label provenance is absent or changed')
    return release

def select_training_graphs(graphs,labels):
    selected={};excluded=[]
    for g in graphs:
        poi=str(g.poi)
        if poi in selected:raise ValidationError('Duplicate GNN POI: '+poi)
        if poi not in labels:
            excluded.append(poi);continue
        value=float(g.y.item())
        if bool(getattr(g,'artifact',False)) or value not in (0.,1.) or value!=labels[poi]:
            raise ValidationError('Graph artifact/label differs from reviewed cohort: '+poi)
        selected[poi]=g
    if set(selected)!=set(labels):raise ValidationError('Reviewed cohort has missing static graphs')
    return [selected[p] for p in sorted(labels)],sorted(excluded)

def normalize_fold_scalars(poi,table,train_pois,columns):
    if table.index.has_duplicates:raise ValidationError('Duplicate scalar POI identifiers')
    if len(set(train_pois))!=len(train_pois):raise ValidationError('Duplicate training POIs')
    values=table.loc[table.index.isin(train_pois),columns].to_numpy(dtype=float)
    if np.isinf(values).any():raise ValidationError('Infinite training scalars')
    # Training-only estimates; all-missing and singleton columns have explicit fallbacks.
    means=np.array([np.mean(x[np.isfinite(x)]) if np.isfinite(x).any() else 0. for x in values.T])
    scales=np.array([np.std(x[np.isfinite(x)],ddof=1) if np.isfinite(x).sum()>1 else 1. for x in values.T])
    scales[~np.isfinite(scales)|(scales==0)]=1.
    raw=table.loc[poi,columns].to_numpy(dtype=float) if poi in table.index else np.full(len(columns),np.nan)
    if np.isinf(raw).any():raise ValidationError('Infinite target scalars')
    norm=(raw-means)/scales;return np.nan_to_num(norm,nan=0.)

def fold_seed(seed,poi):
    return int.from_bytes(hashlib.sha256((str(seed)+'|'+poi).encode()).digest()[:4],'big')%2147483647


def periodic_minima(ligand,protein,box):
    from MDAnalysis.lib.distances import distance_array
    box=np.asarray(box,dtype=float)
    ligand=np.asarray(ligand,dtype=float);protein=np.asarray(protein,dtype=float)
    if box.shape!=(6,) or not np.isfinite(box).all() or np.any(box[:3]<=0) or np.any((box[3:]<=0)|(box[3:]>=180)):
        raise ValidationError('Valid triclinic periodic box required')
    if ligand.ndim!=2 or protein.ndim!=2 or ligand.shape[1:]!=(3,) or protein.shape[1:]!=(3,) or not len(ligand) or not len(protein) or not np.isfinite(ligand).all() or not np.isfinite(protein).all():
        raise ValidationError('Finite nonempty coordinate groups required')
    return distance_array(ligand,protein,box=box).min(axis=0)


def observed_joint_contacts(left,right):
    a=np.asarray(left,dtype=bool);b=np.asarray(right,dtype=bool)
    if a.ndim!=2 or b.ndim!=2 or not len(a) or len(a)!=len(b):raise ValidationError('Matched nonempty frame/contact matrices required')
    return a.astype(np.uint64).T@b.astype(np.uint64)/len(a)


def require_verified_ligand_graph(universe,smiles):
    """Reject the historical connectivity directly in the compiled MD topology."""
    import networkx as nx
    from rdkit import Chem
    expected=Chem.MolFromSmiles(smiles)
    if expected is None or Chem.MolToInchiKey(expected)!='SCUZPTNDZHEHLB-ZWUDBOSPSA-N':
        raise ValidationError('Registry chemistry is not verified SJF8240')
    ligand=universe.select_atoms('resname UNL PTC and mass 2.1:1000')
    if len(ligand)!=expected.GetNumHeavyAtoms() or len(ligand.residues)!=1:
        raise ValidationError('Compiled topology has an ambiguous or incomplete SJF8240 ligand')
    graph=nx.Graph();elements=[(6,12.01),(7,14.01),(8,16.),(9,19.),(16,32.06)]
    for a in ligand:
        z,mass=min(elements,key=lambda x:abs(x[1]-a.mass))
        if abs(mass-a.mass)>.2:raise ValidationError('Unknown heavy atom mass in ligand')
        graph.add_node(int(a.index),z=z)
    for b in ligand.bonds:
        i,j=map(int,b.indices)
        if i in graph and j in graph:graph.add_edge(i,j)
    target=nx.Graph();target.add_nodes_from((a.GetIdx(),{'z':a.GetAtomicNum()}) for a in expected.GetAtoms())
    target.add_edges_from((b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in expected.GetBonds())
    if not nx.is_isomorphic(graph,target,node_match=lambda a,b:a['z']==b['z']):
        raise ValidationError('Compiled ligand connectivity differs from verified SJF8240')
    return ligand


def configure_native_pyg():
    """Use native PyTorch CPU operators; exclude incompatible optional binaries.

    Process-local configuration only. Does not alter installed packages or files.
    This must be applied before importing torch_geometric.
    """
    import sys
    if 'torch_geometric' in sys.modules and not globals().get('_NATIVE_PYG_CONFIGURED',False):raise ValidationError('Configure native PyG before importing it')
    disabled=['pyg_lib','torch_scatter','torch_sparse','torch_cluster','torch_spline_conv','torch_frame','intel_extension_for_pytorch']
    for name in disabled:
        if name in sys.modules and sys.modules[name] is not None:raise ValidationError('Optional compiled GNN backend already loaded: '+name)
        sys.modules[name]=None
    globals()['_NATIVE_PYG_CONFIGURED']=True
    return disabled


def validate_feature_graph(graph,width,schema='SJF8240-periodic-contacts-v1'):
    if getattr(graph,'feature_schema',None)!=schema:raise ValidationError('Historical or incompatible graph feature schema')
    def array(x):return x.detach().cpu().numpy() if hasattr(x,'detach') else np.asarray(x)
    x=array(graph.x);ei=array(graph.edge_index);ea=array(graph.edge_attr)
    if x.ndim!=2 or x.shape[1]!=width or not len(x) or ei.ndim!=2 or ei.shape[0]!=2 or ea.shape!=(ei.shape[1],3) or not np.isfinite(x).all() or not np.isfinite(ea).all():
        raise ValidationError('Invalid or nonfinite GNN graph tensors')
    if not np.issubdtype(ei.dtype,np.integer) or (ei.size and (ei.min()<0 or ei.max()>=len(x))):
        raise ValidationError('GNN graph index outside node range')
