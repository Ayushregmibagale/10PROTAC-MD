"""Positive controls for linker generation at deposited bound-head geometries.

No protein docking/MD and no native ligand coordinates used to seed embedding.
Native head coordinates are supplied because this isolates linker construction.
"""
import argparse,collections,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from rdkit import Chem

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    sys.path.insert(0,str(a.engine));from constrained_linker import constrained_linker,match_heads,stereo_key
    a.output.mkdir(parents=True,exist_ok=False);rows=[];negatives=[]
    for ref in ('7PI4_GSK215','5T35_MZ1'):
        source=a.inputs/ref
        full=Chem.MolFromSmiles((source/'full.smi').read_text().strip())
        heads=[Chem.SDMolSupplier(str(source/f),removeHs=True)[0] for f in ('vhl_head.sdf','poi_head.sdf')]
        if full is None or any(h is None for h in heads):raise ValueError('Invalid control graphs')
        matches=match_heads(full,heads);targets=[h.GetConformer().GetPositions() for h in heads]
        for seed in (42,7,123):
            start=time.monotonic()
            mol,report=constrained_linker(full,heads,targets,seed,attempts=20,embed_iterations=20)
            row={'reference':ref,'seed':seed,'seconds':time.monotonic()-start,'accepted':mol is not None,
                 'sampling':report,'rejection_counts':dict(collections.Counter(report['rejection_reasons'])),
                 'head_mapping_alternatives':len(matches),
                 'source_hashes':{f:sha(source/f) for f in ('full.smi','vhl_head.sdf','poi_head.sdf')}}
            if mol is not None:
                if stereo_key(mol)!=Chem.MolToInchiKey(full):raise ValueError('Output stereo mismatch')
                indices=sum(report['head_mappings'],[]);target=np.vstack(targets)
                displacement=np.linalg.norm(mol.GetConformer().GetPositions()[indices]-target,axis=1)
                row['max_fixed_head_displacement_A']=float(displacement.max())
                if displacement.max()>1e-5:raise ValueError('Native head constraints moved')
                with Chem.SDWriter(str(a.output/f'{ref}_seed{seed}.sdf')) as w:w.write(mol)
            rows.append(row)
            (a.output/'DIRECT_NATIVE_LINKER_RESULTS.json').write_text(json.dumps({'status':'IN_PROGRESS','rows':rows},indent=2)+'\n')
            print(json.dumps({k:row[k] for k in ('reference','seed','seconds','accepted','rejection_counts')}),flush=True)
        displaced=[targets[0],targets[1]+1000.]
        impossible,diagnostic=constrained_linker(full,heads,displaced,42,attempts=2,embed_iterations=2)
        if impossible is not None:raise ValueError('Physically incompatible head separation was accepted')
        negatives.append({'reference':ref,'translation_applied_to_target_head_A':[1000.,1000.,1000.],
                          'correctly_rejected':True,'sampling':diagnostic})
    result={'status':'COMPLETE','scope':'Native-geometry linker positive controls, not independent docking recovery',
       'engine_sha256':sha(a.engine/'constrained_linker.py'),'source_script_sha256':sha(Path(__file__)),
       'python':sys.version,'rdkit':Chem.rdBase.rdkitVersion,'rows':rows,'incompatible_geometry_controls':negatives,
       'summary':{ref:{'accepted':sum(x['accepted'] for x in rows if x['reference']==ref),'attempted_seeds':3} for ref in ('7PI4_GSK215','5T35_MZ1')},
       'limitations':['Exact deposited head geometry only; this is not a kinase screen or degradation prediction',
         'Failed stochastic embedding does not prove the known native molecule is physically impossible',
         'No native full-ligand coordinates used for initialization; fixed head coordinates intentionally isolate construction',
         'Generated full-ligand protein clash assessment is not included in this isolated construction test']}
    (a.output/'DIRECT_NATIVE_LINKER_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
