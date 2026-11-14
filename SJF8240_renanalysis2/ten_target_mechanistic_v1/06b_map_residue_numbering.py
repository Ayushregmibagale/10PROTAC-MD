"""Annotate modeled residue IDs with sequence-validated canonical positions; no model changes."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1
ROOT=Path(__file__).resolve().parent;BASE=Path('/home/ayushregmibagale/10PROTAC-MD');MECH=BASE/'SJF8240_reanalysis/mechanistic_v3';REF=BASE/'SJF8240_reanalysis/config/role_reference.json'
def main():
 reference=json.loads(REF.read_text());rows=[];checks=[]
 for target in pd.read_csv(ROOT/'mechanistic/ten_target_summary.tsv',sep='\t').target:
  folders=ROOT/'extracted'/target
  first=sorted(folders.glob('*/base_packet.npz'))[0];z=np.load(first);ids=z['resids_D'];observed=''.join(seq1(str(x),custom_map={'HID':'H','HIE':'H','HIP':'H','CYX':'C','MSE':'M'},undef_code='X') for x in z['resnames_D']);uid=reference['targets'][target];canonical=reference['sequences'][uid];aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-3,open_gap_score=-5,extend_gap_score=-.5);alignments=aligner.align(canonical,observed)
  def mapping(a):return {j:i+1 for (a,b),(c,d) in zip(*a.aligned) for i,j in zip(range(a,b),range(c,d))}
  mapped=mapping(alignments[0]);identity=sum(canonical[i-1]==observed[j] for j,i in mapped.items())/len(mapped);assert identity>=.95 and len(mapped)/len(observed)>=.95,target+' insufficient sequence correspondence'
  ambiguous=set()
  try:
   other=mapping(alignments[1]);ambiguous={j for j in mapped if mapped.get(j)!=other.get(j)}
  except IndexError:pass
  for j,resid in enumerate(ids):rows.append({'target':target,'uniprot':uid,'topology_resid':int(resid),'canonical_resid':mapped.get(j) if j not in ambiguous else None,'amino_acid':observed[j],'canonical_amino_acid':canonical[mapped[j]-1] if j in mapped else None,'is_lysine':observed[j]=='K','sequence_matches':j in mapped and canonical[mapped[j]-1]==observed[j],'ambiguous':j in ambiguous})
  checks.append({'target':target,'uniprot':uid,'modeled_residues':len(ids),'canonical_first':min(mapped.values()),'canonical_last':max(mapped.values()),'full_canonical_length':len(canonical),'observed_sequence_identity':identity,'unmapped_modeled_residues':len(observed)-len(mapped),'sequence_mismatches':sum(canonical[i-1]!=observed[j] for j,i in mapped.items()),'ambiguous_positions':len(ambiguous)})
 df=pd.DataFrame(rows);df.to_csv(ROOT/'mechanistic/residue_numbering_map.tsv',sep='\t',index=False);df[~df.sequence_matches].to_csv(ROOT/'mechanistic/construct_sequence_differences.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(ROOT/'mechanistic/modeled_protein_coverage.tsv',sep='\t',index=False);keys=df[df.is_lysine & df.sequence_matches & ~df.ambiguous][['target','topology_resid','canonical_resid','uniprot']].rename(columns={'topology_resid':'lysine_resid'})
 for name in ['per_pose_lysine.tsv','per_lysine_equal_pose_summary.tsv','top5_lysines_each_target.tsv']:
  path=ROOT/'mechanistic'/name;a=pd.read_csv(path,sep='\t').drop(columns=['canonical_resid','uniprot','canonical_mapping_verified'],errors='ignore');oldrows=len(a);a=a.merge(keys,on=['target','lysine_resid'],how='left',validate='many_to_one');assert len(a)==oldrows;a['canonical_mapping_verified']=a.canonical_resid.notna();a.to_csv(path,sep='\t',index=False)
 (ROOT/'validation/residue_numbering.json').write_text(json.dumps({'status':'PASS','source_sha256':hashlib.sha256(REF.read_bytes()).hexdigest(),'checks':checks,'numbering_added_only_as_metadata':True},indent=2));print(pd.DataFrame(checks).to_string(index=False))
if __name__=='__main__':main()
