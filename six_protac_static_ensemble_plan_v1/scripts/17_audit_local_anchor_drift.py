"""Read-only audit of linker-anchor drift across all global/local refinements."""
import argparse,csv,json,hashlib
from pathlib import Path
import numpy as np

def anchor_name(path,serial):
 rows=[x for x in path.read_text().splitlines() if x.startswith('HETATM') and int(x[6:11])==serial]
 if len(rows)!=1:raise ValueError('Ambiguous anchor serial')
 return rows[0][12:16].strip()
def distance(path,names):
 xyz=[]
 for res,name in zip(('PT0','PT1'),names):
  rows=[x for x in path.read_text().splitlines() if x.startswith(('ATOM  ','HETATM')) and x[17:20]==res and x[12:16].strip()==name]
  if len(rows)!=1:raise ValueError('Ambiguous/missing named anchor: '+str(path))
  xyz.append(np.array([float(rows[0][i:i+8]) for i in (30,38,46)]))
 return float(np.linalg.norm(xyz[0]-xyz[1]))
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('runs','output'):p.add_argument('--'+k,type=Path,required=True)
 a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True);allrows=[];summary=[]
 for ref in ('5T35_MZ1','7PI4_GSK215'):
  raw=a.runs/f'global_{ref}_seed42/engine_raw';params=dict(x.split(':',1) for x in (raw/'Protac_params.txt').read_text().splitlines() if ':' in x)
  serials=list(map(int,params['Anchor atoms'].split()));names=[anchor_name(raw/f'Init{i}.pdb',serial) for i,serial in enumerate(serials)]
  lo,hi=json.loads((raw/'anchor_distance_sampling.json').read_text())['range_A'];rows=[]
  for i in range(1,51):
   initial=distance(raw/f'Patchdock_Results/pd.{i}.pdb',names)
   for j in range(1,11):
    f=raw/f'Patchdock_Results/pd.{i}_docking_{j:04d}.pdb';d=distance(f,names)
    ledger=raw/f'Patchdock_Results/confs_{i}_{j}.sdf.sampling.json';status=json.loads(ledger.read_text())['status']
    row={'reference':ref,'global_index':i,'local_index':j,'initial_anchor_distance_A':initial,'local_anchor_distance_A':d,'outside_sampled_range':not(lo-.003<=d<=hi+.003),'linker_status':status};rows.append(row);allrows.append(row)
  summary.append({'reference':ref,'anchor_atom_names':names,'range_A':[lo,hi],'local_poses':len(rows),'outside_sampled_range':sum(x['outside_sampled_range'] for x in rows),'median_local_distance_A':float(np.median([x['local_anchor_distance_A'] for x in rows])),'max_local_distance_A':max(x['local_anchor_distance_A'] for x in rows),'feasible_linker_poses':sum(x['linker_status']=='FEASIBLE_CONFORMERS_GENERATED' for x in rows)})
 with (a.output/'LOCAL_ANCHOR_DRIFT.tsv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
 report={'status':'READ_ONLY_ANCHOR_DRIFT_AUDIT_COMPLETE','summary':summary,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limitations':['Sampled range is a finite free-conformer heuristic, not a physical impossibility bound','No native pose used for selection or restraint generation','Distance feasibility alone does not establish full linker closure or recovery']}
 (a.output/'LOCAL_ANCHOR_DRIFT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
