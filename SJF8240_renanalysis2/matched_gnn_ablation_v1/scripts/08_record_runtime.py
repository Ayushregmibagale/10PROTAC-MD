"""Record runtime versions with reviewed process-local nativePyGguard."""
from pathlib import Path
import sys,platform
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'lib'));from runtime import configure_native_pyg
DISABLED=configure_native_pyg()
import torch,torch_geometric,numpy,scipy,pandas
from contracts import atomic_json
if __name__=='__main__':
 atomic_json(ROOT/'validation/runtime_versions.json',{'python':sys.version,'platform':platform.platform(),'torch':str(torch.__version__),'torch_geometric':str(torch_geometric.__version__),'numpy':numpy.__version__,'scipy':scipy.__version__,'pandas':pandas.__version__,'training_device':'CPU','training_threads_per_process':1,'training_workers':2,'optional_native_modules_disabled':DISABLED,'installed_packages_modified':False});print('Savedruntimeversions;CPUonly,noinstalledpackagesmodified.')
