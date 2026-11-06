#!/usr/bin/env python3
"""Pack explicitly selected, hash-bound static graph packets for fresh training."""
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,sha256,atomic_json
from gnn_contract import configure_native_pyg
configure_native_pyg()
import torch
from torch_geometric.data import Data

def pack(folders,output):
    if output.exists():raise ValidationError('Output already exists')
    if any(output.resolve().is_relative_to(Path(p)) for p in ['/home/ayushregmibagale/protac-md','/mnt/d/PROTAC degradation Prediction']):raise ValidationError('Output points into original work')
    graphs=[];seen=set();sources={}
    for folder in folders:
        r=json.loads((folder/'receipt.json').read_text());p=folder/'graph.npz'
        if r.get('status')!='GRAPH_PACKET_REQUIRES_REVIEW' or r.get('feature_schema')!='SJF8240-periodic-contacts-v1' or r.get('graph_sha256')!=sha256(p):raise ValidationError('Unbound or incompatible static graph packet')
        poi=r['poi']
        if poi in seen:raise ValidationError('Duplicate protein target in graph packet')
        seen.add(poi)
        with np.load(p,allow_pickle=False) as a:
            x=a['x'];ei=a['edge_index'];ea=a['edge_attr'];label=int(a['label'][0])
            if x.ndim!=2 or x.shape[1]!=9 or not len(x) or ei.ndim!=2 or ei.shape[0]!=2 or ea.shape!=(ei.shape[1],3) or not np.isfinite(x).all() or not np.isfinite(ea).all():raise ValidationError('Invalid static graph tensor contract')
            if label not in (-1,0,1) or label!=r['label'] or (ei.size and (ei.min()<0 or ei.max()>=len(x))):raise ValidationError('Invalid graph labels/indices')
            graphs.append(Data(x=torch.tensor(x,dtype=torch.float32),edge_index=torch.tensor(ei,dtype=torch.long),edge_attr=torch.tensor(ea,dtype=torch.float32),y=torch.tensor([label],dtype=torch.long),poi=poi,artifact=False,feature_schema='SJF8240-periodic-contacts-v1'))
        sources[str(folder/'receipt.json')]=sha256(folder/'receipt.json')
    if not graphs:raise ValidationError('An explicit nonempty static graph selection is required')
    output.mkdir(parents=True);torch.save(graphs,output/'static_graphs.pt');atomic_json(output/'receipt.json',{'status':'GRAPH_PACKET_REQUIRES_REVIEW','feature_schema':'SJF8240-periodic-contacts-v1','poi':sorted(seen),'source_receipts_sha256':sources,'output_sha256':sha256(output/'static_graphs.pt'),'source_sha256':sha256(Path(__file__)),'automatic_training_release':False})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--graph',type=Path,action='append',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();pack(a.graph,a.output)
