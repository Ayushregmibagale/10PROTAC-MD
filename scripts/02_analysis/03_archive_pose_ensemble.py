#!/usr/bin/env python3
"""Losslessly archive baseline poses as shared PDB templates + integer coordinates.

Coordinates retain exact PDB precision (0.001 Angstrom), atom fields, non-atom
records and source SHA256. Restore validates the complete original bytes.
Never deletes source files. NumPy arrays use no pickle/object payloads.
"""
import argparse,gzip,hashlib,json,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError,sha256,atomic_json
NEGATIVE_ZERO=np.iinfo(np.int32).min

def encode(path):
    raw=path.read_bytes();text=raw.decode('utf-8')
    if '\r' in text:raise ValidationError('Only LF PDB records supported; source retained')
    lines=text.splitlines();template=[];other=[];coords=[]
    for index,line in enumerate(lines):
        if not line.startswith(('ATOM  ','HETATM')):other.append([index,line]);continue
        fields=[line[30:38],line[38:46],line[46:54]];values=[]
        for field in fields:
            value=round(float(field)*1000)
            if not -2147483647<=value<=2147483647:raise ValidationError('Coordinate outside int32 range')
            if value==0 and field.strip().startswith('-'):value=NEGATIVE_ZERO
            restored=-0.0 if value==NEGATIVE_ZERO else value/1000
            if f'{restored:8.3f}'!=field:raise ValidationError('Nonstandard coordinate precision/format; source retained')
            values.append(value)
        coords.append(values);template.append([index,line[:30],line[54:]])
    if not coords:raise ValidationError('No PDB atom records')
    return template,other,np.array(coords,dtype=np.int32),len(lines),text.endswith('\n'),hashlib.sha256(raw).hexdigest()

def reconstruct(template,other,xyz,nlines,newline,patches=None):
    if patches:
        template=list(template)
        for ordinal,index,prefix,suffix in patches:template[ordinal]=[index,prefix,suffix]
    lines=[None]*nlines
    if len(template)!=len(xyz):raise ValidationError('Archive atom/template count mismatch')
    for (index,prefix,suffix),point in zip(template,xyz):
        fields=''.join(f'{(-0.0 if int(v)==NEGATIVE_ZERO else int(v)/1000):8.3f}' for v in point)
        lines[index]=prefix+fields+suffix
    for index,line in other:lines[index]=line
    if any(line is None for line in lines):raise ValidationError('Missing archived PDB records')
    return ('\n'.join(lines)+('\n' if newline else '')).encode('utf-8')

def pack(source,destination,chunk_size=100):
    if chunk_size<1:raise ValidationError('Chunk size must be positive')
    candidates=[p for p in source.iterdir() if p.is_file() and
        re.fullmatch(r'(?:combined_\d+_\d+_\d{4}|pd\.\d+|pd\.\d+_docking_\d{4})\.pdb',p.name)]
    candidates=sorted(candidates)
    if not candidates:raise ValidationError('No baseline pose files')
    destination.mkdir(parents=True,exist_ok=False);templates=[];template_lookup={};batch=[];records=[];chunks=[];original=0
    metadata=gzip.open(destination/'poses.jsonl.gz','wt',encoding='utf-8')
    def flush():
        if not batch:return
        filename=f'coordinates_{len(chunks):06d}.npz'
        np.savez_compressed(destination/filename,xyz_milliAngstrom=np.stack(batch))
        chunks.append({'file':filename,'sha256':sha256(destination/filename),'frames':len(batch),'atoms':len(batch[0])})
        for frame,record in enumerate(records):
            record.update(chunk=filename,frame=frame);metadata.write(json.dumps(record,separators=(',',':'))+'\n')
        batch.clear();records.clear()
    for p in candidates:
        template,other,xyz,nlines,newline,digest=encode(p)
        # Histidine tautomer H names can vary after packing. Sparse atom-field
        # changes should not duplicate an entire protein template per pose.
        key=len(template)
        if key not in template_lookup:template_lookup[key]=len(templates);templates.append(template)
        template_id=template_lookup[key]
        patches=[[i,*row] for i,row in enumerate(template) if row!=templates[template_id][i]]
        if batch and (len(batch)>=chunk_size or xyz.shape!=batch[0].shape):flush()
        record={'name':p.name,'source_sha256':digest,'template':template_id,'atom_field_patches':patches,
                'other_records':other,'n_lines':nlines,'trailing_newline':newline}
        # Verify reconstructability before committing any frame.
        if hashlib.sha256(reconstruct(templates[template_id],other,xyz,nlines,newline,patches)).hexdigest()!=digest:
            raise ValidationError('Lossless reconstruction mismatch: '+p.name)
        batch.append(xyz);records.append(record);original+=p.stat().st_size
    flush();metadata.close()
    manifest={'schema_version':2,'format':'PDB templates plus chunked int32 coordinates, sparse field patches and gzip JSONL',
        'coordinate_units':'0.001 Angstrom','negative_zero_sentinel':int(NEGATIVE_ZERO),'source_directory':str(source.resolve()),
        'source_files':len(candidates),'original_bytes':original,'templates':templates,'chunks':chunks,
        'metadata_sha256':sha256(destination/'poses.jsonl.gz'),'deletes_source_files':False}
    with gzip.open(destination/'manifest.json.gz','wt',encoding='utf-8') as f:json.dump(manifest,f,separators=(',',':'))
    return {'poses':len(candidates),'original_bytes':original,
            'archive_bytes':sum(p.stat().st_size for p in destination.iterdir() if p.is_file()),
            'lossless_frame_reconstruction_verified':True,'source_deleted':False}

def restore(archive,destination):
    if (archive/'manifest.json.gz').exists():
        with gzip.open(archive/'manifest.json.gz','rt',encoding='utf-8') as f:manifest=json.load(f)
    else:manifest=json.loads((archive/'manifest.json').read_text())
    if sha256(archive/'poses.jsonl.gz')!=manifest['metadata_sha256']:raise ValidationError('Archive metadata changed')
    for chunk in manifest['chunks']:
        if Path(chunk['file']).name!=chunk['file'] or sha256(archive/chunk['file'])!=chunk['sha256']:
            raise ValidationError('Archive chunk changed or unsafe path')
    destination.mkdir(parents=True,exist_ok=False);count=0;loaded_name=None;loaded=None
    with gzip.open(archive/'poses.jsonl.gz','rt',encoding='utf-8') as f:
        for line in f:
            record=json.loads(line)
            if Path(record['name']).name!=record['name']:raise ValidationError('Unsafe archived filename')
            if record['chunk']!=loaded_name:
                if record['chunk'] not in {c['file'] for c in manifest['chunks']}:raise ValidationError('Unknown coordinate chunk')
                with np.load(archive/record['chunk'],allow_pickle=False) as data:loaded=data['xyz_milliAngstrom']
                loaded_name=record['chunk']
            raw=reconstruct(manifest['templates'][record['template']],record['other_records'],loaded[record['frame']],record['n_lines'],record['trailing_newline'],record.get('atom_field_patches'))
            if hashlib.sha256(raw).hexdigest()!=record['source_sha256']:raise ValidationError('Restored source SHA256 mismatch')
            path=destination/record['name']
            if path.exists():raise ValidationError('Duplicate archived filename')
            path.write_bytes(raw);count+=1
    if count!=manifest['source_files']:raise ValidationError('Archive pose count mismatch')
    return {'restored':count,'all_source_sha256_verified':True}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['pack','restore'])
    p.add_argument('--source',type=Path,required=True);p.add_argument('--destination',type=Path,required=True)
    p.add_argument('--chunk-size',type=int,default=100);a=p.parse_args()
    result=pack(a.source,a.destination,a.chunk_size) if a.mode=='pack' else restore(a.source,a.destination)
    print(json.dumps(result,indent=2))
