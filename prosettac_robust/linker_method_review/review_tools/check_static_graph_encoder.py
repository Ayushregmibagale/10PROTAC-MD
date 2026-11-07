"""Check the untrained static encoder on corrected-chemistry diagnostic geometry."""
from pathlib import Path
import sys,json,hashlib,importlib.util
ROOT=Path(__file__).resolve().parents[2]
LIB=ROOT/'workflow/scripts/lib';sys.path.insert(0,str(LIB))
from gnn_contract import configure_native_pyg,validate_feature_graph
configure_native_pyg()
import torch
folder=ROOT/'linker_method_review/static_graph_packet_DG_control/graphs.pt'
receipt=json.loads((folder/'receipt.json').read_text());packet=folder/'static_graphs.pt'
assert hashlib.sha256(packet.read_bytes()).hexdigest()==receipt['output_sha256']
graphs=torch.load(packet,map_location='cpu',weights_only=False)
assert len(graphs)==1 and int(graphs[0].y.item())==-1
validate_feature_graph(graphs[0],9)
source=ROOT/'workflow/scripts/04_gnn/01_sjf8240_e2_lopo.py'
spec=importlib.util.spec_from_file_location('diagnostic_encoder',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
torch.set_num_threads(1);torch.use_deterministic_algorithms(True);torch.manual_seed(7919)
model=m.EnsembleEncoder().eval()
with torch.no_grad():embedding=model(graphs[0])
assert tuple(embedding.shape)==(1,128) and torch.isfinite(embedding).all()
result={'status':'UNTRAINED_ENCODER_INTEGRATION_PASSED','input_label':-1,'nodes':int(graphs[0].x.shape[0]),'directed_edges':int(graphs[0].edge_index.shape[1]),'embedding_shape':list(embedding.shape),'finite_embedding':True,'seed':7919,'packet_sha256':receipt['output_sha256'],'encoder_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'torch_version':torch.__version__,'training_performed':False,'biological_prediction_made':False,'automatic_training_release':False}
(ROOT/'linker_method_review/evidence/corrected_static_graph_encoder.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
