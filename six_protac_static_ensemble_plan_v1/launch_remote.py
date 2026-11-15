from pathlib import Path
import json,subprocess,sys
root=Path('/home/s1638308/ernie_static_reference_pilot_v1')
if (root/'STATUS.json').exists() or (root/'LAUNCH_RECEIPT.json').exists():
    raise SystemExit('Existing launch/attempt found. Inspect; do not relaunch blindly.')
with (root/'supervisor.log').open('x') as log:
    command=[sys.executable,str(root/'run_pilot.py'),'--root',str(root)]
    proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,cwd=root)
receipt={'pid':proc.pid,'command':command,'scope':'bounded static reference pilot','MD':False,
         'STOP_command':'touch /home/s1638308/ernie_static_reference_pilot_v1/STOP'}
(root/'LAUNCH_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
