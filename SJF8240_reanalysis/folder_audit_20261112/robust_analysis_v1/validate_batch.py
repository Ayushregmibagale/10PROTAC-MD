"""Validate every requested compound without fitting or silently skipping any."""
from pathlib import Path
import json,argparse
from contracts import *
def validate(path):
 path=Path(path).resolve();batch=json.loads(path.read_text());require(batch['expected_compound_ids'] and len(set(batch['expected_compound_ids']))==len(batch['expected_compound_ids']),'Unique expected compounds required');rows=[];keys=[]
 for row in batch['studies']:
  study=bound_file(path.parent,row);c,labels,graphs=validate_study(study);keys.append((c['compound_id'],c['model'],c['evaluation'],c['training']['loss']));rows.append({'compound_id':c['compound_id'],'model':c['model'],'evaluation':c['evaluation'],'labelled_targets':len(labels),'unknown_targets':len(c['targets'])-len(labels),'negatives':sum(v==0 for v in labels.values()),'folds_per_role':len(folds(labels,c['evaluation'])),'roles':list(c['roles']),'status':'VALIDATED_NOT_TRAINED'})
 require(len(keys)==len(set(keys)),'Duplicate compound/model study');require(set(batch['expected_compound_ids'])=={r['compound_id'] for r in rows},'Missing or unexpected compound');return rows
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--batch',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();atomic_json(a.output,{'status':'ALL_CONFIGURATIONS_VALIDATED_NO_FITTING','studies':validate(a.batch)})
