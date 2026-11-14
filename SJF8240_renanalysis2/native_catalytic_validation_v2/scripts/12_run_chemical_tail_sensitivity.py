"""Run/replay secondary modeled-tail comparison from already prepared frozen inputs."""
from pathlib import Path
import argparse,json,subprocess,sys,io,hashlib
import numpy as np
from Bio.PDB import MMCIFParser
ROOT=Path(__file__).resolve().parents[1];CHILD=ROOT/'chemical_tail_sensitivity'
def check_reference_build():
 old=np.load(ROOT/'references/8RX0_G76_fit75.npz');data={key:old[key] for key in old.files};z=np.load(ROOT/'references/8RX0_locally_relaxed_G76.npz');prior=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_renanalysis2/crl2_ubiquitination_v1');model=MMCIFParser(QUIET=True,auth_chains=False).get_structure('r',str(prior/'references/8RX0.cif'))[0];po=model['I'][(' ',75,' ')]['O'].coord;ii=np.flatnonzero(np.linalg.norm(data['xyz']-po,axis=1)<1e-4);assert len(ii)==1;data['xyz'][ii[0]]=z['G75_O'];data['xyz']=np.vstack([data['xyz'],z['G76_xyz'],z['E2_SG']]);data['elements']=np.r_[data['elements'],['N','C','C','O','S']];data['roles']=np.r_[data['roles'],['MODELED_G76']*4+['MODELED_C93_SG']];data['donor_C']=np.array([z['G76_xyz'][2]]);data['donor_tail_names']=np.array(['chemically_constrained_G76C_G75O_relaxed']);buffer=io.BytesIO();np.savez_compressed(buffer,**data);digest=hashlib.sha256(buffer.getvalue()).hexdigest();expected=json.loads((CHILD/'config/execution_manifest.json').read_text())['references/8RX0_relaxed_G76.npz'];assert digest==expected,'Reference rebuild differs from frozeninput';path=CHILD/'references/8RX0_relaxed_G76.npz'
 if not path.exists():path.write_bytes(buffer.getvalue())
 assert hashlib.sha256(path.read_bytes()).hexdigest()==expected
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--extract',action='store_true',help='Reprocess existingMD using preparedchemical-tail reference;not newMD');p.add_argument('--workers',type=int,default=4);a=p.parse_args();checks=json.loads((ROOT/'validation/locally_relaxed_tail_models.json').read_text());assert next(x['record']['validated_geometry_only'] for x in checks['models'] if x['reference']=='8RX0'),'No approved geometric-tail candidate'
 check_reference_build()
 if a.extract:subprocess.run([sys.executable,str(CHILD/'scripts/02_extract_chemical_tail_sensitivity.py'),'--workers',str(a.workers)],check=True)
 assert (CHILD/'completion.json').exists(),'Complete chemical-tail extraction first'
 for name in ['04_compare_chemical_tail_sensitivity.py','05_independent_verification.py']:subprocess.run([sys.executable,str(CHILD/'scripts'/name)],check=True)
if __name__=='__main__':main()
