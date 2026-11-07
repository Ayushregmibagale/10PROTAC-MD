"""Record actual topology protonation choices without claiming a pH assignment."""
from pathlib import Path
import sys,json,hashlib,importlib.util
ROOT=Path(__file__).resolve().parents[2];WORK=ROOT/'linker_method_review/md_handoff_DG_v1/system'
sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import atomic_json,sha256
spec=importlib.util.spec_from_file_location('audit_md_binding',ROOT/'workflow/scripts/01_setup/07_run_md_stage.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
meta=m.verify_preparation(WORK,'em');chains=[];states=[];total=0.;sources={}
for path in sorted(WORK.glob('topol_Protein_chain_*.itp')):
    section='';residues={};atoms={};bonds=[]
    for line in path.read_text().splitlines():
        line=line.split(';',1)[0].strip()
        if line.startswith('['):section=line.strip('[] ').lower();continue
        if not line or line.startswith('#'):continue
        fields=line.split()
        if section=='atoms' and len(fields)>=8:
            atom={'index':int(fields[0]),'residue':fields[2],'resname':fields[3],'name':fields[4],'charge':float(fields[6])}
            atoms[atom['index']]=atom;residues.setdefault((atom['residue'],atom['resname']),[]).append(atom)
        elif section=='bonds' and len(fields)>=2:bonds.append((int(fields[0]),int(fields[1])))
    chain=path.stem.rsplit('_',1)[-1];charge=sum(a['charge'] for a in atoms.values());total+=charge
    chains.append({'chain':chain,'topology_charge_e':charge,'disulfide_bonds':[{'left_residue':atoms[a]['residue'],'right_residue':atoms[b]['residue']} for a,b in bonds if atoms[a]['name']=='SG' and atoms[b]['name']=='SG']})
    for (number,name),group in residues.items():
        if name not in {'HIS','HID','HIE','HIP','ASP','GLU','LYS','ARG','CYS','CYX'}:continue
        names={a['name'] for a in group};state={'chain':chain,'topology_residue':number,'resname':name,'residue_charge_e':sum(a['charge'] for a in group),'hydrogen_names':sorted(n for n in names if n.startswith('H'))}
        if name in {'HIS','HID','HIE','HIP'}:state['histidine_state']='doubly protonated' if {'HD1','HE2'}<=names else 'delta' if 'HD1' in names else 'epsilon' if 'HE2' in names else 'UNRESOLVED'
        states.append(state)
    sources[path.name]=sha256(path)
counts={};section=''
for line in (WORK/'topol.top').read_text().splitlines():
    line=line.split(';',1)[0].strip()
    if line.startswith('['):section=line.strip('[] ').lower();continue
    if section=='molecules' and line and not line.startswith('#'):
        fields=line.split();counts[fields[0]]=int(fields[1])
net=total+meta['ligand_formal_charge']+counts.get('NA',0)-counts.get('CL',0)
assert abs(net)<0.01,net
out={'status':'TOPOLOGY_CHEMICAL_STATE_AUDIT_NOT_PH_VALIDATION','protein_chains':chains,'titratable_residues':states,'ligand_state':meta['ligand_microstate'],'ligand_formal_charge_e':meta['ligand_formal_charge'],'molecule_counts':counts,'system_charge_balance_e':net,'protein_state_assignment':'GROMACS pdb2gmx default standard residue states and geometry-based histidine choices; no explicit titration/pKa study performed','limitations':['Neutral supplier ligand is a reference microstate, not a demonstrated dominant solution state.','Histidine tautomer choices and terminal/domain protonation need scientific review before a production release.','This audit reports the topology used; it does not modify this control.'],'md_preparation_sha256':sha256(WORK/'md_preparation.json'),'protein_topology_sha256':sources,'audit_source_sha256':sha256(Path(__file__))}
atomic_json(ROOT/'linker_method_review/evidence/md_chemical_states.json',out)
print(json.dumps({'protein_charge_e':total,'histidines':sum('histidine_state' in s for s in states),'system_charge_balance_e':net,'counts':counts},indent=2))
