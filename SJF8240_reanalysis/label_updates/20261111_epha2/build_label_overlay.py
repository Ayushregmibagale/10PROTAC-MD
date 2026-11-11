"""Record user-confirmed MD labels, preserving historical model/graph artifacts."""
from pathlib import Path
import json,hashlib,datetime,argparse
import numpy as np
SOURCE=Path('/home/ayushregmibagale/10PROTAC-MD/SJF8240_reanalysis')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(source,labels_path,destination):
 SOURCE=Path(source).resolve();LABELS=Path(labels_path).resolve();DEST=Path(destination).resolve()
 labels=json.loads(LABELS.read_text());mapping=labels['labels']
 assert len(mapping)==8 and {k for k,v in mapping.items() if v==0}=={'ABL1','AXL'} and mapping['EPHA2']==1
 assert sha(SOURCE/'config/protocol.json')==labels['old_protocol_sha256']
 if DEST.exists():raise ValueError('Preserve existing overlay; use a new version')
 DEST.mkdir(parents=True);graphs=[]
 for poi,label in mapping.items():
  for branch in ['static','dynamic']:
   for role in ['elob','vhl']:
    relative=Path('analysis/static')/poi/role/'graph.npz' if branch=='static' else Path('analysis/dynamic/targets')/poi/role/'graph.npz'
    path=SOURCE/relative;original_hash=sha(path)
    receipt=json.loads((SOURCE/'analysis/static'/poi/'receipt.json').read_text()) if branch=='static' else json.loads((path.parent/'receipt.json').read_text())
    expected=receipt['outputs'][role]['graph_sha256'] if branch=='static' else receipt['graph_sha256']
    assert original_hash==expected,'Historical graph hash drift'
    with np.load(path,allow_pickle=False) as loaded:arrays={k:loaded[k].copy() for k in loaded.files}
    assert arrays['label'].shape==(1,) and arrays['label'][0] in (0,1) and str(arrays['poi'][0])==poi
    old=int(arrays['label'][0]);output=None
    if old!=label:
     assert poi=='EPHA2' and old==0 and label==1
     arrays['label']=np.array([label],dtype=arrays['label'].dtype)
     output=DEST/'graphs'/poi/branch/role/'graph.npz';output.parent.mkdir(parents=True,exist_ok=False)
     np.savez_compressed(output,**arrays)
     with np.load(path,allow_pickle=False) as before,np.load(output,allow_pickle=False) as after:
      assert before.files==after.files
      for key in before.files:
       if key!='label':assert before[key].dtype==after[key].dtype and np.array_equal(before[key],after[key]),key
      assert after['label'][0]==1
    assert sha(path)==original_hash,'Historical source was changed'
    graphs.append({'poi':poi,'branch':branch,'role':role,'old_label':old,'current_label':label,'source_relative_path':str(relative),'source_sha256':original_hash,'overlay_relative_path':str(output.relative_to(DEST)) if output else None,'overlay_sha256':sha(output) if output else None,'feature_arrays_unchanged':True})
 manifest={'status':'CURRENT_MD_LABEL_OVERLAY_PREPARED_NOT_TRAINED','recorded_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'label_map_sha256':sha(LABELS),'builder_sha256':sha(__file__),'cohort':list(mapping),'positive_count':sum(mapping.values()),'negative_count':len(mapping)-sum(mapping.values()),'graph_references':graphs,'new_label_graphs':sum(r['overlay_relative_path'] is not None for r in graphs),'new_training_executed':False,'old_weights_predictions_metrics_modified':False,'original_docking_or_MD_modified':False,'required_next_analysis':'A separately versioned fit/evaluation using the corrected labels; do not score old weights as if trained with this correction.'}
 (DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps({k:v for k,v in manifest.items() if k!='graph_references'},indent=2))
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--source',default=str(SOURCE));parser.add_argument('--labels',required=True);parser.add_argument('--output',required=True);args=parser.parse_args()
 main(args.source,args.labels,args.output)
