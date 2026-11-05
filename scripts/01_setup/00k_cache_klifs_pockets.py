#!/usr/bin/env python3
"""Cache official KLIFS pocket maps for an explicit structure-ID download plan."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError, atomic_json, sha256


def validate(rows):
    if not isinstance(rows,list) or sorted(int(r['index']) for r in rows)!=list(range(1,86)):
        raise ValidationError('KLIFS response must contain each of 85 pocket positions once')
    if any(not isinstance(r['Xray_position'],str) or not r.get('KLIFS_position') for r in rows):
        raise ValidationError('Invalid KLIFS residue map')
    import re
    observed=[r['Xray_position'] for r in rows if re.fullmatch(r'-?\d+[A-Za-z]?',r['Xray_position'])]
    if len(observed)!=len(set(observed)):
        raise ValidationError('Observed pocket residues must be one-to-one')
    return rows


def fetch(root,plan):
    request=json.loads(plan.read_text())
    dest=root/'structure_provenance/pocket_maps';dest.mkdir(parents=True,exist_ok=True)
    def one(identifier):
        url=request['url_template'].format(id=int(identifier))
        if not url.startswith('https://klifs.net/api/interactions_match_residues?structure_ID='):
            raise ValidationError('Download plan must point to the official KLIFS residue-map API')
        target=dest/(str(int(identifier))+'.json')
        if target.exists():
            validate(json.loads(target.read_text()))
            receipt=target.with_suffix('.receipt.json')
            if not receipt.exists() or json.loads(receipt.read_text())['sha256']!=sha256(target):
                raise ValidationError('Missing or altered cached-map receipt')
            return {'id':identifier,'status':'CACHED','sha256':sha256(target)}
        error=None
        for retry in range(3):
            try:
                with urllib.request.urlopen(url,timeout=45) as response:data=response.read()
                validate(json.loads(data));target.write_bytes(data)
                atomic_json(target.with_suffix('.receipt.json'),{'url':url,'structure_id':int(identifier),
                    'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    'sha256':sha256(target)})
                return {'id':identifier,'status':'RETRIEVED','sha256':sha256(target)}
            except (OSError,ValueError,KeyError,ValidationError) as exc:
                error=str(exc)
                if retry<2:time.sleep(retry+1)
        return {'id':identifier,'status':'FAILED','reason':error}
    with ThreadPoolExecutor(max_workers=4) as workers:
        rows=list(workers.map(one,request['ids']))
    report={'plan_sha256':sha256(plan),'maps':rows,'failed':sum(r['status']=='FAILED' for r in rows)}
    atomic_json(dest/'download_receipt.json',report)
    print(json.dumps({'maps':len(rows),'failed':report['failed']},indent=2))
    if report['failed']:raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--plan',type=Path,required=True)
    a=p.parse_args();fetch(a.root,a.plan)
