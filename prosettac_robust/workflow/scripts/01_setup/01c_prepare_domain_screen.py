#!/usr/bin/env python3
"""Prepare each explicit compound/domain/conformer without collapsing identities.

This generates inspectable candidates and exclusion records, never releases or
launches the campaign. Both experimental conformers and AF2 alternatives survive.
"""
import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path
spec=importlib.util.spec_from_file_location('placement',Path(__file__).with_name('01_place_warhead_batch.py'))
placement=importlib.util.module_from_spec(spec);spec.loader.exec_module(placement)
from protac_qc import atomic_json,sha256,ValidationError

def build(root):
    inventory=root/'structure_provenance/domains/domain_inventory.tsv'
    rows=list(csv.DictReader(inventory.open(),delimiter='\t'))
    output=root/'prepared_screen'
    if output.exists():raise ValidationError('Screen exists; use a fresh root for a new preparation')
    output.mkdir();results=[]
    for compound in ['GSK215','SJF8240','SJF1521','JNJ1013','SJF_delta']:
        for row in rows:
            pair=dict(row,protac=compound,alignment_method='ce')
            try:
                info=placement.prepare(root,pair,True,output)
                result={'status':'CANDIDATE','pocket_rmsd_A':info['pocket_rmsd_A'],
                        'preparation_sha256':sha256(output/compound/row['run_name']/'preparation.json')}
            except (ValidationError,ValueError,OSError) as exc:
                result={'status':'REJECTED','reason':str(exc)}
            results.append(dict(result,protac=compound,run_name=row['run_name'],gene=row['gene'],
                uniprot_id=row['uniprot_id'],domain_id=row['domain_id'],source=row['source']))
    from collections import Counter
    report={'status':'CANDIDATES_ONLY','production_ready':False,
            'inventory_sha256':sha256(inventory),'pairs':results,
            'counts':dict(Counter(r['status'] for r in results)),
            'by_compound':{c:dict(Counter(r['status'] for r in results if r['protac']==c))
                           for c in sorted({r['protac'] for r in results})}}
    atomic_json(output/'screen_preparation.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='pairs'},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    build(p.parse_args().root)
