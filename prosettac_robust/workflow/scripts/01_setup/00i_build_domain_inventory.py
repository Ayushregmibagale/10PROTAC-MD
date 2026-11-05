#!/usr/bin/env python3
"""Define kinase domains from cached UniProt annotations and verified coordinates.

No last-DFG windows, gene/PDB-only joins, or silent tandem-domain deduplication.
Non-kinase folds and insufficient/low-confidence domains remain explicit exclusions.
"""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner,substitution_matrices
from Bio.SeqUtils import seq1
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import pdb_atoms,ValidationError,sha256,atomic_json

def entries(provenance):
    result={}
    for file in sorted((provenance/'uniprot_batches').glob('*.json')):
        for entry in json.loads(file.read_text())['results']:
            uid=entry['primaryAccession']
            if uid in result and result[uid]!=entry:raise ValidationError('Conflicting UniProt snapshots')
            result[uid]=entry
    return result

def map_sequence(atoms,sequence):
    ca=[a for a in atoms if a['name']=='CA']
    observed=''.join(seq1(a['resname'],undef_code='X') for a in ca)
    aligner=PairwiseAligner();aligner.mode='local'
    aligner.substitution_matrix=substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score=-10;aligner.extend_gap_score=-0.5
    alignment=aligner.align(sequence,observed)[0]
    pairs=[(i,j) for (a,b),(c,d) in zip(*alignment.aligned) for i,j in zip(range(a,b),range(c,d))]
    identity=sum(sequence[i]==observed[j] for i,j in pairs)/len(pairs)
    if identity<0.90:raise ValidationError(f'Source/UniProt identity only {identity:.3f}')
    mapping={(ca[j]['resnum'],ca[j]['icode']):i+1 for i,j in pairs}
    return mapping,identity

def build(root,source_root):
    provenance=root/'structure_provenance'
    proteins=entries(provenance);kinases={}
    for file in (provenance/'kinase_info_batches').glob('*.json'):
        for k in json.loads(file.read_text()):kinases[k['name']]=k
    resolved=pd.read_csv(provenance/'resolved_structure_inventory.tsv',sep='\t')
    source_rows=[]
    for _,row in resolved[resolved.status=='EXACT_COORDINATES_VERIFIED'].iterrows():
        k=kinases.get(row.klifs_name)
        if not k:continue
        cached=Path(row.cache_path)
        cached=cached if cached.is_absolute() else root/cached
        source_rows.append({'uid':k['uniprot'],'gene':k['HGNC'],'source':'KLIFS','structure_path':str(cached),
            'original_path':row.filepath,'chain':row.chain,'structure_id':'KLIFS'+str(int(row.klifs_id)),
            'klifs_name':row.klifs_name,'klifs_id':int(row.klifs_id),'dfg':row.dfg_state,
            'resolution':row.resolution,'sha256':row.verified_source_sha256})
    legacy=pd.read_csv(source_root/'structures/merged_manifest.tsv',sep='\t')
    for _,row in legacy.iterrows():
        if pd.notna(row.af2_path) and Path(row.af2_path).exists() and row.uniprot_id in proteins:
            source_rows.append({'uid':row.uniprot_id,'gene':row.hgnc_gene,'source':'AF2',
                'structure_path':row.af2_path,'chain':'A','structure_id':'AF2','dfg':'unassigned',
                'sha256':sha256(row.af2_path)})
    ready=[];excluded=[];dest=root/'domain_structures';dest.mkdir(exist_ok=True)
    for source in source_rows:
        uid=source['uid'];entry=proteins.get(uid)
        try:
            if not entry:raise ValidationError('No cached UniProt record')
            if 'reviewed (Swiss-Prot)' not in entry['entryType']:raise ValidationError('Unreviewed target')
            domains=[f for f in entry.get('features',[]) if f['type']=='Domain' and f.get('description','').startswith('Protein kinase')]
            if not domains:raise ValidationError('No annotated protein-kinase domain; do not force another fold into this screen')
            all_atoms=pdb_atoms(source['structure_path'])
            atoms=[a for a in all_atoms if a['record']=='ATOM']
            chains={a['chain'] for a in atoms}
            chain=source['chain']
            if chain not in chains:
                if len(chains)!=1:raise ValidationError('Source chain cannot be identified uniquely')
                chain=next(iter(chains))
            atoms=[a for a in atoms if a['chain']==chain]
            # Modified amino acids must not disappear from the domain silently.
            modified=[a for a in all_atoms if a['chain']==chain and a['resname'] in
                      {'SEP','TPO','PTR','MSE','CSO','CSD','KCX','MLY','HYP'}]
            if modified:
                raise ValidationError('Modified protein residues require an explicit chemical-state model: '+
                    ','.join(sorted({a['resname']+str(a['resnum']) for a in modified})))
            mapping,identity=map_sequence(atoms,entry['sequence']['value'])
            for index,domain in enumerate(domains,1):
                start=domain['location']['start']['value'];end=domain['location']['end']['value']
                indices={key:value for key,value in mapping.items() if start<=value<=end}
                coverage=len(indices)/(end-start+1)
                if len(indices)<100 or coverage<0.85:
                    excluded.append(dict(source,domain_index=index,status='DOMAIN_COVERAGE_REJECTED',coverage=coverage));continue
                if len(indices)>400:
                    excluded.append(dict(source,domain_index=index,status='DOMAIN_TOO_LARGE_FOR_KINASE_TRANSFER'));continue
                selected=[a for a in atoms if (a['resnum'],a['icode']) in indices]
                ca=[a for a in selected if a['name']=='CA']
                confidence=[float(a['line'][60:66]) for a in ca]
                if source['source']=='AF2' and (np.mean(confidence)<80 or np.mean(np.array(confidence)<70)>0.20):
                    excluded.append(dict(source,domain_index=index,status='AF2_DOMAIN_CONFIDENCE_REJECTED',mean_plddt=float(np.mean(confidence))));continue
                domain_id=f'kinase{index}_{start}_{end}'
                name=f'{source["gene"]}_{uid}__{domain_id}__{source["structure_id"]}'
                file=dest/(name+'.pdb')
                # Canonical UniProt residue numbers and a dedicated POI chain.
                lines=[]
                offset=start-1 if max(indices.values())>9999 else 0
                for atom in selected:
                    line=atom['line'];resnum=indices[(atom['resnum'],atom['icode'])]-offset
                    lines.append(line[:21]+'D'+f'{resnum:4d}'+' '+line[27:])
                file.write_text('\n'.join(lines)+'\nTER\nEND\n')
                row=dict(source,domain_id=domain_id,domain_index=index,uniprot_id=uid,
                    uniprot_start=start,uniprot_end=end,coverage=coverage,sequence_identity=identity,
                    structure_path=str(file.relative_to(root)),source_structure_path=source['structure_path'],
                    source_chain=chain,chain='D',res_start=min(indices.values())-offset,res_end=max(indices.values())-offset,
                    pdb_numbering_offset=offset,
                    run_name=name,status='DOMAIN_DEFINED',mean_plddt=float(np.mean(confidence)) if source['source']=='AF2' else None)
                row['domain_sha256']=sha256(file)
                atomic_json(file.with_suffix('.provenance.json'),row);ready.append(row)
        except (ValidationError,ValueError,OSError,KeyError) as exc:
            excluded.append(dict(source,status='SOURCE_REJECTED',reason=str(exc)))
    output=provenance/'domains';output.mkdir(exist_ok=True)
    pd.DataFrame(ready).to_csv(output/'domain_inventory.tsv',sep='\t',index=False)
    pd.DataFrame(excluded).to_csv(output/'exclusions.tsv',sep='\t',index=False)
    summary={'defined_domains':len(ready),'unique_uniprot':len({r['uid'] for r in ready}),
             'by_source':dict(Counter(r['source'] for r in ready)),
             'excluded':dict(Counter(r['status'] for r in excluded))}
    atomic_json(output/'summary.json',summary);print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--source-root',type=Path,required=True)
    a=p.parse_args();build(a.root,a.source_root)
