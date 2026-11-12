"""Read-only delivery/cached-input audit. Run before reusing completed packets."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));from extract import require,verify_source,sha256,atomic_json

def main():
 verify_source();completion=json.loads((ROOT/'analysis/completion.json').read_text());require(len(completion['records'])==141,'Incomplete case receipts');bindings={}
 for r in completion['records']:
  require(r['frames']==201,'Incomplete frame count');p=ROOT/'analysis/poses'/r['poi']/r['pose']/'packet.npz';require(sha256(p)==r['packet_sha256'],'Packet drift')
  for path,digest in r['input_sha256'].items():
   require(path not in bindings or bindings[path]==digest,'Conflicting source bindings');bindings[path]=digest
 for path,digest in bindings.items():require(Path(path).is_file() and sha256(Path(path))==digest,'Missing or changed input '+path)
 for path,digest in json.loads((ROOT/'config/reporting_manifest.json').read_text()).items():require(sha256(ROOT/path)==digest,'Reporting source drift')
 atomic_json(ROOT/'validation/current_input_recheck.json',{'status':'PASS','case_packets':141,'unique_external_files_rehashed':len(bindings),'source_inputs_sha256':bindings,'required_before_cached_reuse':True});print('Input audit passed:',len(bindings),'external files and141 packets')
if __name__=='__main__':main()
