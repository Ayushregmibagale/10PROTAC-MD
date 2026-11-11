"""Bounded MD-only GPU control; no study fitting or metric estimate."""
from pathlib import Path
import sys,ast,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from common import *
import torch.nn.functional as F
verify_frozen();p,graphs=load_inputs('vhl','cuda');require(torch.cuda.is_available(),'GPU not available');torch.set_num_threads(1);torch.use_deterministic_algorithms(True);torch.manual_seed(991)
original=ast.parse((ROOT/'provenance/04h_md_only_ablation_original.py').read_text());current=ast.parse((ROOT/'lib/model.py').read_text());extract=lambda t:next(n for n in t.body if isinstance(n,ast.ClassDef) and n.name=='MDEncoder')
require(ast.dump(extract(original),include_attributes=False)==ast.dump(extract(current),include_attributes=False),'Original MD-only architecture AST mismatch')
model=MDEncoder().cuda();opt=torch.optim.Adam(model.parameters(),lr=.001);started=time.monotonic()
for _ in range(20):
 model.train();opt.zero_grad();loss=F.binary_cross_entropy_with_logits(model(graphs['MET']),torch.ones(1,device='cuda'));require(bool(torch.isfinite(loss)),'Nonfinite control loss');loss.backward();opt.step()
score=predict(model,graphs['MET']);checkpoint=ROOT/'logs/preflight_weights.pt';torch.save({'model_state_dict':{k:v.detach().cpu().clone() for k,v in model.state_dict().items()}},checkpoint)
replay=MDEncoder().cuda();replay.load_state_dict(torch.load(checkpoint,map_location='cuda',weights_only=True)['model_state_dict']);error=abs(score-predict(replay,graphs['MET']));require(error<=1e-7,'GPU control replay mismatch')
report={'status':'MD_ONLY_GPU_PREFLIGHT_PASS_NOT_STUDY_RESULT','steps':20,'elapsed_s':time.monotonic()-started,'saved_state_replay_difference':error,'architecture_AST_matches_original_04h':True,'all_eight_evaluation_targets_have_MD':True,'EPHA2_label':p['labels']['EPHA2'],'torch':str(torch.__version__),'GPU':torch.cuda.get_device_name(0),'source_manifest_sha256':verify_frozen()};atomic_json(ROOT/'logs/preflight.json',report);print(json.dumps(report,indent=2))
