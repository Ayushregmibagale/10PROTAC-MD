"""Verify canonical reference-chain identities against deposited entity sequences.

In particular, do not assume ELOB/ELOC have the same B/C ordering in different
VCB sources. Source metadata and observed-sequence alignment are both retained.
"""
import argparse,hashlib,json
from pathlib import Path
from Bio import Align
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def candidates(cif):
    d=MMCIF2Dict(str(cif));polys=dict(zip(d['_entity_poly.entity_id'],d['_entity_poly.pdbx_seq_one_letter_code_can']))
    out=[]
    for entity,db,accession,code in zip(d['_struct_ref.entity_id'],d['_struct_ref.db_name'],d['_struct_ref.pdbx_db_accession'],d['_struct_ref.db_code']):
        if db!='UNP' or entity not in polys:continue
        role={'P40337':'VHL','Q15370':'ELOB','Q15369':'ELOC'}.get(accession.split('-')[0],'POI')
        out.append({'entity_id':entity,'uniprot':accession,'database_code':code,'role':role,
                    'sequence':''.join(polys[entity].split())})
    return out
def match(observed,reference):
    aligner=Align.PairwiseAligner();aligner.mode='local';aligner.match_score=2.;aligner.mismatch_score=-1.
    aligner.open_gap_score=-8.;aligner.extend_gap_score=-.5
    a=aligner.align(reference,observed)[0];same=paired=0
    for (s,e),(q,r) in zip(*a.aligned):
        if e-s!=r-q:raise ValueError('Unexpected nonparallel alignment block')
        paired+=int(e-s);same+=sum(x==y for x,y in zip(reference[s:e],observed[q:r]))
    return {'score':float(a.score),'identity':same/paired if paired else 0.,
            'observed_coverage':paired/len(observed),'source_sequence_coverage':paired/len(reference),
            'aligned_residues':paired,'observed_residues':len(observed)}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--references',type=Path,required=True)
    p.add_argument('--sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    report={}
    for ref,pdb in (('5T35_MZ1','5T35'),('7PI4_GSK215','7PI4')):
        native=a.references/ref/'native_complex.pdb';source=a.sources/(pdb+'.cif');entities=candidates(source)
        model=PDBParser(QUIET=True).get_structure('p',str(native))[0];roles={}
        for chain in ('A','B','C','D'):
            observed=''.join(seq1(r.resname,custom_map={'MSE':'M'}) for r in model[chain] if r.id[0]==' ' and 'CA' in r)
            ranked=sorted([{**{k:v for k,v in e.items() if k!='sequence'},**match(observed,e['sequence'])} for e in entities],key=lambda r:r['score'],reverse=True)
            best=ranked[0]
            if best['identity']<.95 or best['observed_coverage']<.95:raise ValueError('Reference role not supported by observed sequence: '+ref+'/'+chain)
            tied=[x for x in ranked if abs(x['score']-best['score'])<1e-8]
            if len({x['role'] for x in tied})!=1:raise ValueError('Ambiguous role alignment')
            roles[chain]=best
        if roles['A']['role']!='VHL' or roles['D']['role']!='POI' or {roles[c]['role'] for c in ('B','C')}!={'ELOB','ELOC'}:
            raise ValueError('Canonical reference chain roles differ from required selections')
        report[ref]={'chain_roles':roles,'native_sha256':sha(native),'source_cif_sha256':sha(source),
                     'extended_role_chains':{role:[c for c,v in roles.items() if v['role']==role] for role in ('VHL','ELOB','ELOC','POI')}}
    a.output.mkdir(parents=True,exist_ok=True)
    result={'status':'REFERENCE_CHAIN_ROLES_VERIFIED','references':report,
       'interpretation':'Observed reference sequences match deposited entity/UniProt annotations; isoforms are retained explicitly',
       'limitations':['Reference controls only; every production input still needs explicit source/domain/role checks','No assumption of fixed B/C accessory ordering across sources']}
    (a.output/'CHAIN_ROLE_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({ref:d['extended_role_chains'] for ref,d in report.items()},indent=2))

if __name__=='__main__':main()
