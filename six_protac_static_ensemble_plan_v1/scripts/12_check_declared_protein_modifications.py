"""Guard against silently omitting HETATM protein modifications such as TPO.

Explicit residue/atom declarations are required; arbitrary HETATM records are
never guessed to be protein. This supplements, rather than certifies, pose QC.
"""
import argparse,json
from pathlib import Path
import numpy as np

STANDARD={'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
def inspect(path,chains,specs):
    expected={(s['chain'],int(s['resnum']),s.get('icode',' '),s['resname']):set(s['required_heavy_atoms']) for s in specs}
    if len(expected)!=len(specs):raise ValueError('Duplicate modification declaration')
    seen={};protein=[];allkeys=set()
    for line in path.read_text().splitlines():
        if not line.startswith(('ATOM  ','HETATM')) or line[16:17] not in (' ','A'):continue
        if line[21] not in chains:continue
        xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        if not np.isfinite(xyz).all():raise ValueError('Nonfinite protein coordinates')
        res=(line[21],int(line[22:26]),line[26:27],line[17:20].strip());atom=line[12:16].strip()
        key=(*res,atom)
        if key in allkeys:raise ValueError('Duplicate protein atom')
        allkeys.add(key)
        heavy=line[76:78].strip().upper() not in ('H','D')
        if res in expected:
            if heavy:seen.setdefault(res,set()).add(atom);protein.append(xyz)
        elif line.startswith('ATOM  ') and res[3] in STANDARD:
            if heavy:protein.append(xyz)
        else:raise ValueError('Undeclared protein-chain hetero/modified residue: '+str(res))
    for res,names in expected.items():
        if seen.get(res,set())!=names:raise ValueError('Missing/changed modification atom set: '+str(res))
    if not protein:raise ValueError('Empty declared protein selection')
    return {'status':'DECLARED_PROTEIN_MODIFICATIONS_PRESENT','protein_heavy_atoms_including_modifications':len(protein),
            'modifications':[{'chain':r[0],'resnum':r[1],'icode':r[2],'resname':r[3],'heavy_atoms':len(seen[r])} for r in expected],
            'scope':'Residue/atom presence and explicit role selection only; not a physical/biological validation'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pdb',type=Path,required=True);p.add_argument('--spec',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    spec=json.loads(a.spec.read_text());result=inspect(a.pdb,set(spec['protein_chains']),spec['protein_modifications'])
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
