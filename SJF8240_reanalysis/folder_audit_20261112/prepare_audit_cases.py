from pathlib import Path
import sys,json
import MDAnalysis as mda
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'robust_analysis_v1'))
from contracts import *
from build_dynamic import sequence
PROJECT=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis');OLD=Path('/home/ayushregmibagale/protac-md');refs=PROJECT/'config/role_reference.json';labels=json.loads((PROJECT/'config/current_md_labels.json').read_text())['labels'];rolemap=json.loads(refs.read_text())['targets'];records=json.loads((PROJECT/'analysis/dynamic/completion.json').read_text())['records'];cases=[]
for r in records:
 if r['poi'] not in labels:continue
 target,pose=r['poi'],r['pose_id'];tpr=OLD/'random_screen'/target/pose/'prod/prod.tpr';xtc=OLD/'analysis/random_screen'/target/pose/'prod_nojump.xtc';u=mda.Universe(str(tpr));roles={}
 for role,chain,uid in [('elob','A','Q15370'),('vhl','C','P40337'),('poi','D',rolemap[target])]:
  selection='moltype Protein_chain_'+chain+' and name CA';roles[role]={'selection':selection,'uniprot':uid,'model_sequence':sequence(u.select_atoms(selection))}
 cases.append({'target_id':target,'pose_id':pose,'tpr':{'path':str(tpr),'sha256':sha(tpr)},'xtc':{'path':str(xtc),'sha256':sha(xtc)},'protein_roles':roles,'ligand_selection':'resname UNL','ligand_heavy_atoms':79})
c={'compound_id':'SJF8240','units':{'coordinates':'angstrom','time':'ns'},'raw_trajectory_time_unit':'ps','reference_sequences':{'path':str(refs),'sha256':sha(refs)},'target_uniprot':{t:rolemap[t] for t in labels},'expected_times_ns':[i/10 for i in range(201)],'duplicate_policy':'identical_only','contact_A':5.,'edge_frequency_min':.01,'cases':cases};atomic_json(ROOT/'real_extraction_cases.json',c);print(len(cases),'cases prepared')
