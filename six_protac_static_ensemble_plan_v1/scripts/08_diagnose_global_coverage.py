"""Matched-stream global coverage diagnostic; no local Rosetta docking or MD.

Regenerate up to500 global orientations and require exact first50 agreement
with the completed pilot. Native metrics are evaluation-only, never a selector.
"""
import argparse,contextlib,csv,hashlib,importlib.util,io,json,os,shutil,sys,time
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def protein_roles(path):
    # Global sampler concatenates two PDBs; strip intermediate END records for
    # analysis only, without changing any coordinates or the archived source.
    text='\n'.join(x for x in path.read_text().splitlines() if x.startswith(('ATOM  ','HETATM','TER')))+'\nEND\n'
    model=PDBParser(QUIET=True).get_structure('s',io.StringIO(text))[0];out={}
    for role,chain in (('vhl','A'),('poi','D')):
        rs=[r for r in model[chain] if r.id[0]==' ' and 'CA' in r]
        out[role]=([(r.id[1],r.id[2],r.resname) for r in rs],np.array([r['CA'].coord for r in rs],dtype=float),
                   [np.array([a.coord for a in r if a.element not in ('H','D')],dtype=float) for r in rs])
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine',type=Path,required=True);p.add_argument('--pilot',type=Path,required=True)
    p.add_argument('--references',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    for key in ('engine','pilot','references','output'):
        setattr(a,key,getattr(a,key).resolve())
    a.output.mkdir(parents=True,exist_ok=False)
    helper=load('recovery_metric',Path(__file__).with_name('01_assess_reference_recovery.py'))
    original_cwd=Path.cwd();results=[];allrows=[]
    try:
        for ref in ('7PI4_GSK215','5T35_MZ1'):
            source=a.pilot/'runs'/f'global_{ref}_seed42'/'engine_raw'
            output=a.output/ref;output.mkdir()
            for name in ('Init0.pdb','Init1.pdb'):shutil.copy2(source/name,output/name)
            interval=json.loads((source/'anchor_distance_sampling.json').read_text())['range_A']
            params={x.split(':',1)[0]:x.split(':',1)[1].strip() for x in (source/'Protac_params.txt').read_text().splitlines() if ':' in x}
            anchors=[int(x) for x in params['Anchor atoms'].split()]
            os.chdir(output);os.environ.update(PROTAC_ATTEMPT_ROOT=str(output),PROTAC_SEED='42',
                PROTAC_MAX_GLOBAL_TRIALS='20000',PROTAC_PROSETTAC=str(a.engine/'runtime/PRosettaC'))
            os.environ.pop('PROTAC_TASK_KEY',None)
            adapter=load('global_adapter_'+ref,a.engine/'runtime/adapter/pyrosetta_adapter.py')
            start=time.monotonic()
            with (output/'GLOBAL_GENERATION.log').open('w') as f,contextlib.redirect_stdout(f):
                count=adapter.global_dock(['Init0.pdb','Init1.pdb'],anchors,*interval,num_results=500)
            count=count or 0
            if count<50:raise ValueError('Did not regenerate the original50global poses')
            agreement=all(sha(output/'Patchdock_Results'/f'pd.{i}.pdb')==sha(source/'Patchdock_Results'/f'pd.{i}.pdb') for i in range(1,51))
            if not agreement:raise ValueError('Original50global poses differ: do not claim a matched-stream comparison')
            native=protein_roles(a.references/ref/'native_complex.pdb');rows=[]
            for i in range(1,count+1):
                pose=protein_roles(output/'Patchdock_Results'/f'pd.{i}.pdb')
                if any(pose[k][0]!=native[k][0] for k in ('vhl','poi')):raise ValueError('Residue correspondence differs')
                r,t=helper.transform(pose['vhl'][1],native['vhl'][1])
                rmsd=float(np.sqrt(np.mean(np.sum((pose['poi'][1]@r+t-native['poi'][1])**2,axis=1))))
                # Contact evaluation only needed for candidates within the fixed
                # RMSD gate; others cannot pass its conjunction by definition.
                fnat=helper.measure(native,pose)['native_residue_contact_fraction'] if rmsd<=5 else None
                row={'reference':ref,'global_index':i,'POI_CA_RMSD_after_VHL_fit_A':rmsd,
                     'native_residue_contact_fraction_if_RMSD_le5A':fnat,
                     'coarse_orientation_gate':rmsd<=5 and fnat is not None and fnat>=.5}
                rows.append(row);allrows.append(row)
            result={'reference':ref,'sampled_orientations':count,'seconds_generation_and_assessment':time.monotonic()-start,
                'first50_byte_identical_to_completed_pilot':True,
                'scenarios':{str(n):{'best_initial_POI_RMSD_A':min(x['POI_CA_RMSD_after_VHL_fit_A'] for x in rows[:n]),
                         'initial_orientations_passing_coarse_gate':sum(x['coarse_orientation_gate'] for x in rows[:n])}
                         for n in (50,100,200,500) if n<=count},
                'source_sha256':{x:sha(source/x) for x in ('Init0.pdb','Init1.pdb','anchor_distance_sampling.json')},
                'adapter_sha256':sha(a.engine/'runtime/adapter/pyrosetta_adapter.py')}
            results.append(result);print(json.dumps(result,indent=2),flush=True)
    finally:os.chdir(original_cwd)
    with (a.output/'GLOBAL_ORIENTATION_METRICS.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(allrows)
    report={'status':'MATCHED_STREAM_DIAGNOSTIC_COMPLETE','results':results,
       'limitations':['Global C-alpha clash/anchor sampler only, before local refinement or full PROTAC construction',
        'Native geometry used for evaluation only; no native-based pose selection is authorized',
        'One seed/two reference-bound systems; cannot establish a production budget or degradation accuracy',
        'More global coverage is not equivalent to physically valid final ternary ensembles']}
    (a.output/'GLOBAL_COVERAGE_DIAGNOSIS.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
