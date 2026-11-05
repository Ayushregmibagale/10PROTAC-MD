#!/usr/bin/env python3
"""Prepare a separate explicit-pocket comparison; never overwrite or release a screen."""
import argparse
from collections import Counter
import csv
import importlib.util
import json
from pathlib import Path

spec=importlib.util.spec_from_file_location('pocket_screen_placement',Path(__file__).with_name('01_place_warhead_batch.py'))
placement=importlib.util.module_from_spec(spec);spec.loader.exec_module(placement)
from protac_qc import ValidationError,atomic_json,sha256


def build(root):
    root=root.resolve();inventory=root/'structure_provenance/domains/domain_inventory.tsv'
    correspondence=root/'structure_provenance/pocket_correspondence.json'
    maps=json.loads(correspondence.read_text())
    if sha256(inventory)!=maps['inventory_sha256']:raise ValidationError('Pocket map inventory changed')
    rows=list(csv.DictReader(inventory.open(),delimiter='\t'));dest=root/'prepared_pocket_screen'
    if dest.exists():raise ValidationError('Comparison already exists; preserve it and use a fresh root')
    dest.mkdir();results=[]
    for compound,reference in maps['references'].items():
        for row in rows:
            result={'protac':compound,'run_name':row['run_name'],'gene':row['gene'],
                    'uniprot_id':row['uniprot_id'],'domain_id':row['domain_id'],'source':row['source']}
            target=maps['targets'].get(row['run_name'])
            if target is None:
                results.append(dict(result,status='CORRESPONDENCE_UNAVAILABLE',reason='No validated explicit pocket map for this domain'))
                continue
            pair=dict(row,protac=compound,alignment_method='klifs',pocket_correspondence={
                'target':target,'reference':reference,'correspondence_sha256':sha256(correspondence)})
            if compound=='SJF_delta':pair['reference_bundle']=reference['reference_bundle']
            try:
                info=placement.prepare(root,pair,True,dest)
                folder=dest/compound/row['run_name'];info['source_path']=row['structure_path']
                atomic_json(folder/'preparation.json',info)
                results.append(dict(result,status='CANDIDATE',pocket_rmsd_A=info['pocket_rmsd_A'],
                    preparation_sha256=sha256(folder/'preparation.json')))
            except (ValidationError,ValueError,OSError) as exc:
                results.append(dict(result,status='REJECTED',reason=str(exc)))
    report={'status':'COMPARISON_ONLY','production_ready':False,'old_screen_unchanged':True,
        'correspondence_sha256':sha256(correspondence),'pairs':results,
        'counts':dict(Counter(r['status'] for r in results)),
        'by_compound':{c:dict(Counter(r['status'] for r in results if r['protac']==c)) for c in maps['references']},
        'known_target_coverage':{c:[r['run_name'] for r in results if r['protac']==c and r['gene']==gene and r['status']=='CANDIDATE']
             for c,gene in [('GSK215','PTK2'),('SJF8240','MET'),('SJF1521','EGFR'),('JNJ1013','IRAK1'),('SJF_delta','MAPK13')]}}
    atomic_json(dest/'screen_preparation.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='pairs'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    build(p.parse_args().root)
