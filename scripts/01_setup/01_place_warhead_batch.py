#!/usr/bin/env python3
"""Prepare GSK215/SJF8240 candidates using explicit domains and bound head graphs.

TSV fields: protac,gene,uniprot_id,structure_path,chain,res_start,res_end.
Writes separate prepared_inputs; candidates remain blocked from production.
"""
import argparse
import csv
import json
import re
import shutil
import sys
import uuid
from pathlib import Path
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqUtils import seq1
from rdkit import Chem
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from protac_qc import ValidationError, pdb_atoms, min_distance, sha256, atomic_json, read_params

def kabsch(mobile, reference):
    """Row-vector convention: mobile @ R.T + t = reference."""
    mc,rc=mobile.mean(0),reference.mean(0)
    u,_,vt=np.linalg.svd((mobile-mc).T@(reference-rc))
    r=vt.T@np.diag([1,1,np.linalg.det(vt.T@u.T)])@u.T
    return r,rc-mc@r.T

def ce_path_indices(result):
    """Normalize Biopython 1.84 paired lists and 1.85 named CE results."""
    path=getattr(result,'path',result)
    if not isinstance(path,(tuple,list)) or len(path)!=2 or len(path[0])!=len(path[1]) or not len(path[0]):
        raise ValidationError('Malformed CE residue correspondence')
    return tuple(path[0]),tuple(path[1])

def domain_atoms(path, chain, start, end):
    atoms=[a for a in pdb_atoms(path) if a['record']=='ATOM' and a['chain']==chain
           and start <= a['resnum'] <= end]
    ca=[a for a in atoms if a['name']=='CA']
    if not 100<=len(ca)<=400:
        raise ValidationError(f'Explicit kinase domain must have 100–400 Cα atoms, got {len(ca)}')
    return atoms,ca

def align_domain(target, reference, head,method='sequence'):
    aligner=PairwiseAligner()
    aligner.substitution_matrix=substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score=-10; aligner.extend_gap_score=-0.5
    alignment=aligner.align(''.join(seq1(a['resname'],undef_code='X') for a in reference),
                            ''.join(seq1(a['resname'],undef_code='X') for a in target))[0]
    pairs=[(i,j) for (a,b),(c,d) in zip(*alignment.aligned) for i,j in zip(range(a,b),range(c,d))]
    identity=sum(reference[i]['resname']==target[j]['resname'] for i,j in pairs)/len(pairs)
    if method=='ce':
        # CE supplies a structural correspondence across divergent, independently
        # annotated kinase domains. Do not relax pocket/clash/coverage thresholds.
        from Bio.PDB.cealign import run_cealign
        paths=[ce_path_indices(p) for p in run_cealign([list(a['xyz']) for a in reference],[list(a['xyz']) for a in target],8,30)]
        if not paths:raise ValidationError('No CE domain correspondence')
        longest_length=max(len(p[0]) for p in paths)
        longest=[p for p in paths if len(p[0])==longest_length]
        def residual(path):
            a,b=path
            x=np.array([target[j]['xyz'] for j in b]);y=np.array([reference[i]['xyz'] for i in a])
            r,t=kabsch(x,y)
            return float(np.mean(np.sum((x@r.T+t-y)**2,axis=1)))
        best=min(longest,key=residual);pairs=list(zip(*best))
        identity=sum(reference[i]['resname']==target[j]['resname'] for i,j in pairs)/len(pairs)
    elif method!='sequence':raise ValidationError('Unknown domain alignment method')
    if (method=='sequence' and identity<0.25) or len(pairs)<100: raise ValidationError('Insufficient homologous domain correspondence')
    ref=np.array([reference[i]['xyz'] for i,j in pairs]); mob=np.array([target[j]['xyz'] for i,j in pairs])
    keep=np.ones(len(pairs),dtype=bool)
    for _ in range(8):
        r,t=kabsch(mob[keep],ref[keep]); distances=np.linalg.norm(mob@r.T+t-ref,axis=1)
        candidate=distances<=3.0
        if np.array_equal(candidate,keep): break
        if candidate.sum()<100 or candidate.mean()<0.6: raise ValidationError('Too little alignable domain core')
        keep=candidate
    r,t=kabsch(mob[keep],ref[keep]); distances=np.linalg.norm(mob@r.T+t-ref,axis=1)
    ref_pocket={i for i,a in enumerate(reference) if min_distance(np.array([a['xyz']]),head)<10.0}
    pocket_indices=[k for k,(i,j) in enumerate(pairs) if i in ref_pocket]
    if len(pocket_indices)<10 or len(pocket_indices)/len(ref_pocket)<0.8:
        raise ValidationError('Insufficient mapped pocket coverage')
    rmsd=float(np.sqrt(np.mean(distances[keep]**2)))
    pocket_rmsd=float(np.sqrt(np.mean(distances[pocket_indices]**2)))
    if rmsd>2.5 or pocket_rmsd>2.0: raise ValidationError(f'Poor domain/pocket fit: {rmsd:.2f}/{pocket_rmsd:.2f} Å')
    return r,t,{'identity':identity,'alignment_method':method,'core_rmsd_A':rmsd,'pocket_rmsd_A':pocket_rmsd,
                'pocket_coverage':len(pocket_indices)/len(ref_pocket),
                'aligned_residues':[{'reference':[reference[i]['resnum'],reference[i]['icode']],
                                     'target':[target[j]['resnum'],target[j]['icode']],
                                     'fit_core':bool(keep[k])} for k,(i,j) in enumerate(pairs)]}

def align_klifs_pocket(target,reference,head,evidence):
    """Fit explicit pocket identities; independently check the domain core.

    The original geometry thresholds remain: 80% of the head neighborhood,
    >=100 core residues, >=60% core, <=2.5 Å core and <=2 Å pocket RMSD.
    Missing/ambiguous pocket positions never become guessed correspondences.
    """
    def ca_index(atoms):
        index={}
        for i,a in enumerate(atoms):
            key=(a['resnum'],a['icode'].strip())
            if key in index:raise ValidationError('Ambiguous Cα identity in pocket alignment')
            index[key]=i
        return index
    ti,ri=ca_index(target),ca_index(reference)
    tp,rp=evidence['target']['positions'],evidence['reference']['positions']
    for mapping in [tp,rp]:
        if any(str(k) not in {str(i) for i in range(1,86)} for k in mapping):
            raise ValidationError('Invalid KLIFS position index')
        if len(mapping.values())!=len(set(mapping.values())):
            raise ValidationError('Pocket residue correspondence must be one-to-one')
    positions=sorted(set(tp)&set(rp),key=int)
    pairs=[(ri[(int(rp[k]),'')],ti[(int(tp[k]),'')]) for k in positions
           if (int(rp[k]),'') in ri and (int(tp[k]),'') in ti]
    if len(pairs)<68:raise ValidationError('Insufficient explicit 85-position pocket coverage')
    x=np.array([target[j]['xyz'] for i,j in pairs]);y=np.array([reference[i]['xyz'] for i,j in pairs])
    r,t=kabsch(x,y);error=np.linalg.norm(x@r.T+t-y,axis=1)
    rmsd=float(np.sqrt(np.mean(error**2)))
    neighborhood={i for i,a in enumerate(reference) if min_distance(np.array([a['xyz']]),head)<10.0}
    local=[k for k,(i,j) in enumerate(pairs) if i in neighborhood]
    coverage=len(local)/len(neighborhood) if neighborhood else 0
    if len(local)<10 or coverage<0.8:raise ValidationError('Insufficient mapped head-neighborhood coverage')
    pocket_rmsd=float(np.sqrt(np.mean(error[local]**2)))
    if rmsd>2.0 or pocket_rmsd>2.0:raise ValidationError(f'Poor explicit pocket fit: {rmsd:.2f}/{pocket_rmsd:.2f} Å')
    from Bio.PDB.cealign import run_cealign
    paths=run_cealign([list(a['xyz']) for a in reference],[list(a['xyz']) for a in target],8,30)
    checks=[]
    for path in paths:
        a,b=ce_path_indices(path)
        distances=np.linalg.norm(np.array([target[j]['xyz'] for j in b])@r.T+t-
                                 np.array([reference[i]['xyz'] for i in a]),axis=1)
        keep=distances<=3.0
        if keep.sum()>=100 and keep.mean()>=0.6:
            checks.append((float(np.sqrt(np.mean(distances[keep]**2))),int(keep.sum()),float(keep.mean())))
    if not checks:raise ValidationError('Insufficient independently alignable domain core after pocket fit')
    core_rmsd,core_atoms,core_fraction=min(checks)
    if core_rmsd>2.5:raise ValidationError('Domain core RMSD exceeds 2.5 Å after pocket fit')
    return r,t,{'alignment_method':'explicit_KLIFS_pocket','pocket_rmsd_A':pocket_rmsd,
        'whole_KLIFS_pocket_rmsd_A':rmsd,'pocket_coverage':coverage,'core_rmsd_A':core_rmsd,
        'core_atoms':core_atoms,'core_fraction':core_fraction,'mapped_KLIFS_positions':len(pairs),
        'pocket_correspondence_evidence':evidence,
        'aligned_residues':[{'reference':[reference[i]['resnum'],reference[i]['icode']],
                             'target':[target[j]['resnum'],target[j]['icode']]} for i,j in pairs]}


def prepare(root,row,allow_modeled=False,output_root=None):
    modeled=row['protac'] in {'SJF1521','JNJ1013','SJF_delta'}
    if row['protac'] not in {'GSK215','SJF8240'} and not (modeled and allow_modeled):
        raise ValidationError('No validated atom-mapped head for '+row['protac'])
    refdir=(root/'modeled_reference_candidates' if modeled else root/'reference_structures/validated')/row['protac']
    if row.get('reference_bundle'):
        candidate=(root/row['reference_bundle']).resolve()
        if not modeled or not candidate.is_relative_to((root/'modeled_reference_candidates').resolve()):
            raise ValidationError('Explicit reference bundle must be a modeled candidate under this project')
        refdir=candidate
    heads=json.loads((refdir/'heads.json').read_text())
    source=Path(row['structure_path'])
    source=(source if source.is_absolute() else root/source).resolve()
    atoms,ca=domain_atoms(source,row['chain'],int(row['res_start']),int(row['res_end']))
    _,refca=domain_atoms(refdir/'protein.pdb','D',-999,9999)
    mol=Chem.MolFromMolFile(str(refdir/'poi_head.sdf'),removeHs=False)
    if row.get('alignment_method')=='klifs':
        evidence=row['pocket_correspondence']
        if sha256(source)!=evidence['target']['domain_sha256'] or sha256(refdir/'protein.pdb')!=evidence['reference']['protein_sha256']:
            raise ValidationError('Pocket correspondence coordinate hash mismatch')
        r,t,meta=align_klifs_pocket(ca,refca,mol.GetConformer().GetPositions(),evidence)
    else:
        r,t,meta=align_domain(ca,refca,mol.GetConformer().GetPositions(),row.get('alignment_method','sequence'))
    placed=(mol.GetConformer().GetPositions()-t)@r
    distance=min_distance(placed,np.array([a['xyz'] for a in atoms if a['element'] not in ('H','D')]))
    if not 1.5 <= distance <= 4.5: raise ValidationError(f'Implausible head contacts/clash: {distance:.3f} Å')
    name=row.get('run_name') or row['gene']+'_'+row['uniprot_id']
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',name) or any('/' in v or '\\' in v or v in {'.','..'} for v in [row['gene'],row['uniprot_id']]):
        raise ValidationError('Unsafe run identifier')
    dest=(Path(output_root) if output_root is not None else root/'prepared_inputs')/row['protac']/name
    if dest.exists(): raise ValidationError('Candidate exists; use a new root for a new attempt')
    dest.mkdir(parents=True)
    for i,xyz in enumerate(placed): mol.GetConformer().SetAtomPosition(i,xyz)
    w=Chem.SDWriter(str(dest/'poi_head.sdf')); w.write(mol); w.close()
    (dest/'protein.pdb').write_text('\n'.join(a['line'][:21]+'D'+a['line'][22:] for a in atoms)+'\nTER\nEND\n')
    protein_lines=(refdir/'proteins.pdb').read_text().splitlines()
    (dest/'vbc.pdb').write_text('\n'.join(l for l in protein_lines if l.startswith('ATOM  ') and l[21] in 'ABC')+'\nEND\n')
    shutil.copy2(refdir/'vhl_head.sdf',dest/'vhl_head.sdf')
    shutil.copy2(root/'molecules'/row['protac']/'full.smi',dest/'full.smi')
    (dest/'Protac_params.txt').write_text('Structures: vbc.pdb protein.pdb\nChains: ABC D\nHeads: vhl_head.sdf poi_head.sdf\n'
        +f'Anchor atoms: {heads["vhl"]["anchor_one_based"]} {heads["poi"]["anchor_one_based"]}\nProtac: full.smi\nFull: True\nClusterName: LOCAL\n')
    read_params(dest/'Protac_params.txt')
    selection={'vhl':{'chains':['A'],'record':'ATOM'},'poi':{'chains':['D'],'record':'ATOM'},
               'protac':{'chains':['X'],'record':'HETATM','resname':'PTC'},
               'qc':{'min_nonbonded_heavy_distance_A':1.5,'units':'Angstrom',
                     'basis':'Same minimum contact distance as bound-head preparation; gross-overlap exclusion only'}}
    atomic_json(dest/'selection.json',selection)
    meta.update({'schema_version':1,'status':'PREPARATION_CANDIDATE','production_ready':False,
                 'protac':row['protac'],'gene':row['gene'],'uniprot_id':row['uniprot_id'],
                 'run_name':name,'domain_id':row.get('domain_id'),
                 'structure_source':row.get('source'),
                 'reference_kind':heads.get('method','native_bound_fragment'),
                 'reference_bundle':str(refdir.relative_to(root)),
                 'reference_model_evidence':heads.get('model_evidence'),
                 'contains_modeled_head':modeled,'head_bundle_sha256':sha256(refdir/'heads.json'),
                 'source_path':str(source),'source_sha256':sha256(source),'source_chain':row['chain'],
                 'domain':[int(row['res_start']),int(row['res_end'])], 'head_min_distance_A':distance,
                 'reference_sha256':sha256(refdir/'poi_head.sdf'),
                 'molecule_inchikey':Chem.MolToInchiKey(Chem.MolFromSmiles((dest/'full.smi').read_text().strip())),
                 'files':{p.name:sha256(p) for p in dest.iterdir() if p.is_file()}})
    atomic_json(dest/'preparation.json',meta)
    return meta

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--pair-table',type=Path,required=True)
    p.add_argument('--run',action='store_true',help='Write candidates; default lists explicit domains')
    p.add_argument('--allow-modeled-candidates',action='store_true',
                   help='Explicitly admit analog-model references for QC, never a production release')
    a=p.parse_args(); rows=list(csv.DictReader(a.pair_table.open(),delimiter='\t'))
    if not a.run: print(json.dumps(rows,indent=2)); return
    results=[]
    for row in rows:
        try:
            meta=prepare(a.root,row,a.allow_modeled_candidates)
            results.append({'protac':row['protac'],'gene':row['gene'],'uniprot_id':row['uniprot_id'],
                            'status':'CANDIDATE','rmsd':meta['pocket_rmsd_A']})
        except (ValidationError,ValueError,OSError) as exc:
            results.append({'protac':row['protac'],'gene':row['gene'],'uniprot_id':row['uniprot_id'],
                            'status':'REJECTED','reason':str(exc)})
    atomic_json(a.root/'audits/preparation_attempts'/(str(uuid.uuid4())+'.json'),
                {'pair_table_sha256':sha256(a.pair_table),'allow_modeled':a.allow_modeled_candidates,'results':results})
    summary=[]
    for path in sorted((a.root/'prepared_inputs').glob('*/*/preparation.json')):
        info=json.loads(path.read_text())
        summary.append({'protac':path.parent.parent.name,'run_name':path.parent.name,
                        'status':info['status'],'production_ready':info['production_ready'],
                        'contains_modeled_head':info.get('contains_modeled_head',False),
                        'pocket_rmsd_A':info['pocket_rmsd_A'],'metadata_sha256':sha256(path)})
    atomic_json(a.root/'prepared_inputs/preparation_summary.json',summary)
    print(json.dumps(results,indent=2))
    if any(r['status']=='REJECTED' for r in results): raise SystemExit(1)

if __name__=='__main__': main()
