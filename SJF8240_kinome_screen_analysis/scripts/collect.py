"""Collect ranked candidates only after every inventoried pool finishes."""
from pathlib import Path
import json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def main():
 pools=json.loads((ROOT/'config/pools.json').read_text()); tables=[]; counts=[]
 for pool in pools:
  directory=ROOT/'results/pools'/pool['key']; receipt=directory/'completion.json'
  if not receipt.exists():raise SystemExit('Screen incomplete: '+pool['key'])
  r=json.loads(receipt.read_text());counts.append({'key':pool['key'],'status':r['status'],'error':r.get('error','')})
  if r['status']=='CONTEXT_REJECTED':continue
  if r['status']!='COMPLETE':raise SystemExit('Unexpected status: '+pool['key'])
  table=directory/'poses.tsv'
  if hashlib.sha256(table.read_bytes()).hexdigest()!=r['poses_sha256']:raise SystemExit('Changed results: '+str(table))
  if r['parsed']+r['rejected']!=r['expected']:raise SystemExit('Incomplete accounting: '+pool['key'])
  if r['parsed']:tables.append(pd.read_csv(table,sep='\t'))
 if not tables:raise SystemExit('No successfully parsed poses')
 df=pd.concat(tables,ignore_index=True)
 df['worst_pocket_min_A']=df[['ligand_VHL_pocket_min_A','ligand_POI_pocket_min_A']].max(axis=1)
 df=df.sort_values(['passes_full_screen','balanced_contact_atoms','worst_pocket_min_A','path'],ascending=[False,False,True,True],kind='stable')
 df.to_csv(ROOT/'results/all_static_ranked_poses.tsv',sep='\t',index=False)
 good=df[df.passes_full_screen.eq(True)].drop_duplicates('pdb_sha256')
 prospective=good[good.source_label.astype(str).eq('unknown')]
 for name,frame in [('top10_static_all_poses',good.head(10)),('top10_static_prospective_poses',prospective.head(10)),('top10_static_distinct_prospective_kinases',prospective.drop_duplicates('uid').head(10))]:
  for row in frame.itertuples():
   if hashlib.sha256(Path(row.path).read_bytes()).hexdigest()!=row.pdb_sha256:raise SystemExit('Selected source drift: '+row.path)
  frame.to_csv(ROOT/'results'/f'{name}.tsv',sep='\t',index=False)
 pd.DataFrame(counts).to_csv(ROOT/'results/pool_status.tsv',sep='\t',index=False)
 print(f'{len(df)} parsed poses; {len(good)} unique QC-passing poses; {prospective.uid.nunique()} prospective kinases. Missing/rejected pools are not negative results.')
if __name__=='__main__':main()
