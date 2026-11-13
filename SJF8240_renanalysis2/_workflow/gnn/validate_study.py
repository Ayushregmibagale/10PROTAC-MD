from pathlib import Path
import argparse,json
from contracts import *
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);a=p.parse_args();c,labels,graphs=validate_study(a.study);print(json.dumps({'status':'VALIDATED_NOT_TRAINED','compound_id':c['compound_id'],'model':c['model'],'labelled_targets':len(labels),'prediction_only_targets':len(c['targets'])-len(labels),'negative_targets':sum(y==0 for y in labels.values()),'roles':list(c['roles']),'folds_per_role':len(folds(labels,c['evaluation']))}))
