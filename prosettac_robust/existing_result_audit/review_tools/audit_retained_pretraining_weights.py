"""Tensor-only inspection of retained MGDB pretraining checkpoints; no fitting."""
from pathlib import Path
import json,sys
import torch
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,atomic_json,ValidationError
project=Path('/home/ayushregmibagale/protac-md');output=ROOT/'existing_result_audit/retained_gnn_weights_v1'
if output.exists():raise ValidationError('Fresh evidence output required')
paths=[project/'gnn/models/pretrain_best.pt',project/'gnn/models/pretrain_final.pt',project/'gnn/models/training_log.tsv',project/'scripts/04_gnn/03_pretrain_gnn.py',project/'scripts/02_analysis/04f_finetune_from_mgdb.py',project/'scripts/02_analysis/04g_e2_scalar_gnn.py'];hashes={str(p):sha256(p) for p in paths};records=[]
for p in paths[:2]:
    record={'path':str(p),'sha256':hashes[str(p)],'bytes':p.stat().st_size}
    try:
        checkpoint=torch.load(p,map_location='cpu',weights_only=True);state=checkpoint['model_state']
        if not isinstance(state,dict) or not state or any(not isinstance(t,torch.Tensor) or not torch.isfinite(t).all().item() for t in state.values()):raise ValidationError('Invalid/nonfinite tensor checkpoint state')
        record.update(status='TENSOR_ONLY_CHECKPOINT_READABLE',epoch=checkpoint.get('epoch'),tensor_shapes={k:list(v.shape) for k,v in state.items()},tensor_entries=sum(t.numel() for t in state.values()),metadata_keys=sorted(checkpoint),args=checkpoint.get('args',{}))
    except Exception as exc:record.update(status='TENSOR_ONLY_INSPECTION_FAILED_CLOSED',reason=str(exc))
    records.append(record)
if any(sha256(p)!=hashes[str(p)] for p in paths):raise ValidationError('Checkpoint/source drift')
output.mkdir();atomic_json(output/'checkpoints.json',records)
report={'status':'RETAINED_PRETRAINING_WEIGHTS_AUDITED_NOT_HISTORICAL_E2_FOLD_REPLAY','checkpoints':len(records),'tensor_readable':sum(r['status']=='TENSOR_ONLY_CHECKPOINT_READABLE' for r in records),'declared_pretraining_model':'MGDB InterfaceGNN, 24-dimensional node inputs and 1-dimensional edge inputs; used by separate 04f adapter experiment','main_E2_model':'04g E2DualEncoderGNN has separate 9-feature static/3-feature dynamic encoders plus four E2 scalars; no pretrained-checkpoint load in that source','interpretation':'These saved pretrained weights exist and must not be described as missing, but they do not supply the missing trained per-fold states for the historical 17-row E2 GNN evaluation. No model was instantiated or fitted.','source_sha256':hashes,'audit_sha256':sha256(Path(__file__)),'records_sha256':sha256(output/'checkpoints.json'),'original_files_modified':False,'new_GNN_training':False}
atomic_json(output/'summary.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['source_sha256']},indent=2))
