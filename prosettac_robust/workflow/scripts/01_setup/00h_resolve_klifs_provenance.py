#!/usr/bin/env python3
"""Recover exact KLIFS IDs by verifying source coordinates, never by name alone.

The original selection rule nominates one ID per downloaded file; it is accepted
only if fresh ID-specific coordinates have an identical atom fingerprint. The
download cache and full result table preserve failures and mismatches explicitly.
"""
import argparse
import concurrent.futures
import hashlib
import json
import sys
import urllib.request
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import atomic_json, sha256

def fingerprint(text):
    atoms=sorted((line[:6],line[12:27],line[30:54],line[76:78].strip())
                 for line in text.splitlines() if line.startswith(('ATOM  ','HETATM')))
    if not atoms:raise ValueError('No source atoms')
    return hashlib.sha256(json.dumps(atoms).encode()).hexdigest()

def resolve(source,output,workers=4):
    metadata=pd.read_csv(source/'structures/klifs_manifest.tsv',sep='\t')
    nominations=metadata.sort_values(['curated','missing_residues','resolution'],
                  ascending=[False,True,True]).drop_duplicates(['gene','dfg_state'])
    audit=pd.read_csv(source/'audits/klifs_metadata_audit.tsv',sep='\t')
    output=Path(output).resolve();cache=output/'coordinate_cache';cache.mkdir(parents=True,exist_ok=True)
    def one(row):
        candidates=nominations[(nominations.gene==row['klifs_name']) &
                    (nominations.pdb_id.str.lower()==row['pdb_id'].lower()) & (nominations.chain==row['chain'])]
        result={k:row[k] for k in ['klifs_name','pdb_id','chain','filepath']}
        result['original_sha256']=sha256(row['filepath'])
        if len(candidates)!=1:return dict(result,status='UNRESOLVED_NOMINATION')
        chosen=candidates.iloc[0];sid=int(chosen['klifs_id'])
        url=f'https://klifs.net/api/structure_get_pdb_complex?structure_ID={sid}'
        try:
            file=cache/f'{sid}.pdb'
            if not file.exists():
                request=urllib.request.Request(url,headers={'User-Agent':'10PROTAC-MD-provenance-audit/1.0'})
                with urllib.request.urlopen(request,timeout=30) as response:body=response.read()
                # Validate before caching; interrupted/error responses cannot become valid inputs.
                fingerprint(body.decode());file.write_bytes(body)
            expected=fingerprint(Path(row['filepath']).read_text())
            found=fingerprint(file.read_text())
            result.update(klifs_id=sid,dfg_state=chosen['dfg_state'],resolution=chosen['resolution'],
                          missing_residues=chosen['missing_residues'],quality_score=chosen['quality_score'],
                          source_url=url,verified_source_sha256=sha256(file),atom_fingerprint=found,
                          cache_path=str(file.relative_to(output.parent)))
            result['status']='EXACT_COORDINATES_VERIFIED' if expected==found else 'COORDINATE_MISMATCH'
        except Exception as exc:result.update(status='FETCH_OR_VALIDATION_FAILED',reason=str(exc))
        return result
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        rows=list(pool.map(one,audit.to_dict('records')))
    table=pd.DataFrame(rows);table.to_csv(output/'resolved_structure_inventory.tsv',sep='\t',index=False)
    summary=table.status.value_counts().to_dict();atomic_json(output/'summary.json',summary)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args()
    if not 1<=a.workers<=4:p.error('Use 1-4 concurrent source requests')
    resolve(a.source_root,a.output,a.workers)
