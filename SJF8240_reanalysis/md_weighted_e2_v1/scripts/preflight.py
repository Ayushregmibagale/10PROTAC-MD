"""Bounded MD-only GPU control; no study fitting or metric estimate."""
from pathlib import Path
import sys,ast,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'))
from common import *
import torch.nn.functional as F
verify_frozen();p,graphs=load_inputs('vhl','cuda');require(torch.cuda.is_available(),'GPU not available');torch.set_num_threads(1);torch.use_deterministic_algorithms(True);torch.manual_seed(991)
require(original_gat_backbone_matches(),'Original MD GAT backbone differs');apply_e2_features(graphs,p,p['cohort'][1:])
model=MDEncoder().cuda();opt=torch.optim.Adam(model.parameters(),lr=.001);started=time.monotonic()
for _ in range(20):
 model.train();opt.zero_grad();loss=F.binary_cross_entropy_with_logits(model(graphs['MET']),torch.ones(1,device='cuda'));require(bool(torch.isfinite(loss)),'Nonfinite control loss');loss.backward();opt.step()
score=predict(model,graphs['MET']);checkpoint=ROOT/'logs/preflight_weights.pt';torch.save({'model_state_dict':{k:v.detach().cpu().clone() for k,v in model.state_dict().items()}},checkpoint)
replay=MDEncoder().cuda();replay.load_state_dict(torch.load(checkpoint,map_location='cuda',weights_only=True)['model_state_dict']);error=abs(score-predict(replay,graphs['MET']));require(error<=1e-7,'GPU control replay mismatch')
report={'status':'MD_ONLY_GPU_PREFLIGHT_PASS_NOT_STUDY_RESULT','steps':20,'elapsed_s':time.monotonic()-started,'saved_state_replay_difference':error,'GAT_backbone_matches_original':True,'two_E2_graph_scalars_used':True,'all_eight_evaluation_targets_have_MD':True,'EPHA2_label':p['labels']['EPHA2'],'torch':str(torch.__version__),'GPU':torch.cuda.get_device_name(0),'source_manifest_sha256':verify_frozen()};atomic_json(ROOT/'logs/preflight.json',report);print(json.dumps(report,indent=2))
