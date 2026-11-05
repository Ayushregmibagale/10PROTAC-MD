#!/usr/bin/env python3
"""Map observed KLIFS pockets to canonical kinase domains without guessing gaps.

AF2 inherits pocket indices only from consistent experimental maps of the same
UniProt domain. Domains without such evidence remain explicitly unavailable.
"""
import argparse
import csv
import importlib.util
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError, atomic_json, pdb_atoms, sha256


def source_residue(value):
    match=re.fullmatch(r'(-?\d+)([A-Za-z]?)',str(value).strip())
    return (int(match[1]),match[2] or ' ') if match else None


def consensus(maps):
    """Never choose between contradictory canonical identities."""
    candidates=defaultdict(set)
    for mapping in maps:
        if len(mapping.values())!=len(set(mapping.values())):
            raise ValidationError('A canonical residue cannot represent two pocket positions')
        for position,residue in mapping.items():
            candidates[position].add(residue)
    conflicts={k:sorted(v) for k,v in candidates.items() if len(v)>1}
    if conflicts:
        raise ValidationError('Conflicting canonical pocket positions: '+json.dumps(conflicts,sort_keys=True))
    result={k:next(iter(v)) for k,v in candidates.items()}
    if len(result.values())!=len(set(result.values())):
        raise ValidationError('Consensus pocket positions are not one-to-one')
    return result


def build(root):
    root=root.resolve();here=Path(__file__).resolve().parent
    spec=importlib.util.spec_from_file_location('pocket_domain_mapping',here/'00i_build_domain_inventory.py')
    domains=importlib.util.module_from_spec(spec);spec.loader.exec_module(domains)
    provenance=root/'structure_provenance';entries=domains.entries(provenance)
    inventory=provenance/'domains/domain_inventory.tsv'
    rows=list(csv.DictReader(inventory.open(),delimiter='\t'))
    pocket_dir=provenance/'pocket_maps';targets={};experimental=defaultdict(list);errors=[]
    for row in rows:
        if row['source']!='KLIFS':continue
        try:
            identifier=int(float(row['klifs_id']));path=pocket_dir/f'{identifier}.json'
            receipt=json.loads(path.with_suffix('.receipt.json').read_text())
            if sha256(path)!=receipt['sha256']:raise ValidationError('Pocket map hash mismatch')
            source=root/row['source_structure_path']
            atoms=[a for a in pdb_atoms(source) if a['record']=='ATOM' and a['chain']==row['source_chain']]
            mapping,identity=domains.map_sequence(atoms,entries[row['uniprot_id']]['sequence']['value'])
            start,end=int(row['uniprot_start']),int(row['uniprot_end'])
            canonical={str(p['index']):mapping[key] for p in json.loads(path.read_text())
                       if (key:=source_residue(p['Xray_position'])) in mapping
                       and start<=mapping[key]<=end}
            ca={a['resnum'] for a in pdb_atoms(root/row['structure_path']) if a['name']=='CA'}
            offset=int(row['pdb_numbering_offset'])
            positions={k:v-offset for k,v in canonical.items() if v-offset in ca}
            targets[row['run_name']]={'status':'EXPLICIT_EXPERIMENTAL_MAP','positions':positions,
                'pocket_map_sha256':sha256(path),'domain_sha256':sha256(root/row['structure_path']),
                'map_structure_id':identifier,'canonical_positions':canonical}
            experimental[(row['uniprot_id'],row['domain_id'])].append((row['run_name'],canonical))
        except (OSError,ValueError,KeyError,ValidationError) as exc:
            errors.append({'run_name':row['run_name'],'reason':str(exc)})
    for row in rows:
        if row['source']!='AF2':continue
        sources=experimental.get((row['uniprot_id'],row['domain_id']),[])
        try:
            if not sources:raise ValidationError('No experimentally mapped pocket for this exact UniProt domain')
            canonical=consensus([mapping for _,mapping in sources]);offset=int(row['pdb_numbering_offset'])
            ca={a['resnum'] for a in pdb_atoms(root/row['structure_path']) if a['name']=='CA'}
            positions={k:v-offset for k,v in canonical.items() if v-offset in ca}
            targets[row['run_name']]={'status':'CANONICAL_MAP_FROM_SAME_DOMAIN_EXPERIMENTAL_EVIDENCE',
                'positions':positions,'domain_sha256':sha256(root/row['structure_path']),
                'evidence_run_names':[name for name,_ in sources], 'canonical_positions':canonical}
        except (OSError,ValueError,KeyError,ValidationError) as exc:
            errors.append({'run_name':row['run_name'],'reason':str(exc)})
    references={}
    for compound,identifier,bundle in [
        ('GSK215',13914,'reference_structures/validated/GSK215'),
        ('SJF8240',2998,'reference_structures/validated/SJF8240'),
        ('SJF1521',844,'modeled_reference_candidates/SJF1521'),
        ('JNJ1013',8700,'modeled_reference_candidates/JNJ1013'),
        ('SJF_delta',2949,'modeled_reference_candidates/SJF_delta_MAPK13_published_overlay')]:
        path=pocket_dir/f'{identifier}.json';protein=root/bundle/'protein.pdb'
        observed={(a['resnum'],a['icode']) for a in pdb_atoms(protein) if a['name']=='CA' and a['chain']=='D'}
        positions={str(p['index']):key[0] for p in json.loads(path.read_text())
                   if (key:=source_residue(p['Xray_position'])) in observed and not key[1].strip()}
        if len(positions)<68:raise ValidationError('Reference has insufficient explicit pocket positions: '+compound)
        references[compound]={'positions':positions,'reference_bundle':bundle,
            'protein_sha256':sha256(protein),'pocket_map_sha256':sha256(path),'structure_id':identifier}
    report={'schema_version':1,'status':'CORRESPONDENCES_ONLY_NOT_A_CAMPAIGN_RELEASE',
        'inventory_sha256':sha256(inventory),'targets':targets,'references':references,
        'unavailable':errors,'counts':{'domains':len(rows),'mapped':len(targets),'unavailable':len(errors)}}
    atomic_json(provenance/'pocket_correspondence.json',report)
    print(json.dumps(report['counts'],indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    build(parser.parse_args().root)
