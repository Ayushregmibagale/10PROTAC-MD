"""Evaluate all saved global/local orientations; no generation or native-based selection."""
import argparse,csv,importlib.util,json
from pathlib import Path
import numpy as np

def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('runs','references','output'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 root=Path(__file__).parent;g=load('global_reader',root/'08_diagnose_global_coverage.py');h=load('recovery',root/'01_assess_reference_recovery.py')
 rows=[];summary=[]
 for ref in ('5T35_MZ1','7PI4_GSK215'):
  raw=a.runs/f'global_{ref}_seed42'/'engine_raw'/'Patchdock_Results'
  native=g.protein_roles(a.references/ref/'native_complex.pdb');rr=[]
  for i in range(1,501):
   row={'reference':ref,'global_index':i}
   for stage,name in [('global',f'pd.{i}.pdb'),('local',f'pd.{i}_docking_0001.pdb')]:
    source=raw/name
    if not source.exists():raise ValueError('Missing expected orientation: '+str(source))
    pose=g.protein_roles(source)
    if any(pose[k][0]!=native[k][0] for k in ('vhl','poi')):raise ValueError('Residue correspondence differs')
    r,t=h.transform(pose['vhl'][1],native['vhl'][1])
    rmsd=float(np.sqrt(np.mean(np.sum((pose['poi'][1]@r+t-native['poi'][1])**2,axis=1))))
    row[stage+'_RMSD_A']=rmsd
    row[stage+'_sha256']=h.sha(source)
   row['local_minus_global_RMSD_A']=row['local_RMSD_A']-row['global_RMSD_A'];rr.append(row);rows.append(row)
  summary.append({'reference':ref,'paired_orientations':len(rr),'best_global_RMSD_A':min(x['global_RMSD_A'] for x in rr),'best_local_RMSD_A':min(x['local_RMSD_A'] for x in rr),'global_RMSD_le5A':sum(x['global_RMSD_A']<=5 for x in rr),'local_RMSD_le5A':sum(x['local_RMSD_A']<=5 for x in rr),'improved_RMSD_count':sum(x['local_minus_global_RMSD_A']<0 for x in rr),'median_RMSD_change_A':float(np.median([x['local_minus_global_RMSD_A'] for x in rr]))})
 with (a.output/'REFINEMENT_TRACE.tsv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 report={'status':'PAIRED_REFINEMENT_TRACE_COMPLETE','summary':summary,'limitations':['Native used only for retrospective evaluation; no poses selected for further docking','RMSD-only diagnostic; not frozen conjunction gate or full-linker/steric certification','One seed, 500 global x1 local; cannot establish convergence'],'source_sha256':h.sha(Path(__file__))}
 (a.output/'REFINEMENT_TRACE.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
