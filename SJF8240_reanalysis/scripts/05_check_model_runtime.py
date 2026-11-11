"""Bounded execution/replay control using retained MET graphs; not study results."""
from pathlib import Path
import os,sys,importlib.util,time,json
os.environ['CUBLAS_WORKSPACE_CONFIG']=':4096:8'
r=Path(__file__).resolve().parents[1];spec=importlib.util.spec_from_file_location('training',r/'scripts/04_train_dual_lopo.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
torch=m.torch;device='cuda';torch.set_num_threads(1);torch.use_deterministic_algorithms(True);torch.manual_seed(991)
sg=m.load_graph(r/'analysis/static/MET/vhl/graph.npz','MET',1,9).to(device);dg=m.load_graph(r/'analysis/dynamic/targets/MET/vhl/graph.npz','MET',1,3).to(device)
model=m.DualEncoderGNN(hidden=64,heads=4,dropout=.3).to(device);opt=torch.optim.Adam(model.parameters(),lr=.001);t=time.monotonic()
for _ in range(20):
 model.train();opt.zero_grad();loss=m.F.binary_cross_entropy_with_logits(model(sg,dg),torch.ones(1,device=device));loss.backward();opt.step()
torch.cuda.synchronize();elapsed=time.monotonic()-t
first=m.predict(model,sg,dg);state={k:v.detach().cpu().clone() for k,v in model.state_dict().items()};replay=m.DualEncoderGNN(hidden=64,heads=4,dropout=.3).to(device);replay.load_state_dict(state);second=m.predict(replay,sg,dg)
if abs(first-second)>1e-7:raise RuntimeError('Model replay differs')
report={'status':'BOUNDED_GPU_CONTROL_PASS_NOT_STUDY_RESULT','training_steps':20,'elapsed_s':elapsed,'step_mean_s':elapsed/20,'score_replay_error':abs(first-second),'torch':str(torch.__version__),'gpu':torch.cuda.get_device_name(0),'peak_GPU_memory_bytes':torch.cuda.max_memory_allocated(),'deterministic_algorithms':True,'inputs':'Retained MET corrected graphs; not an experimental evaluation','new_MD_or_docking':False}
(r/'logs/gpu_runtime_control.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
