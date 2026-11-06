#!/usr/bin/env python3
"""Launch the preserved local docking overlay; licensed dependencies stay private."""
import json
import os
from pathlib import Path
import sys

def runtime_command(root, arguments):
    private=root/'.runtime'
    config=json.loads((private/'config.json').read_text())
    if config.get('schema_version')!=1:
        raise ValueError('Unsupported local runtime configuration')
    python=Path(config['python'])
    relative=Path(config['overlay_relative'])
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('Overlay must stay inside the private runtime directory')
    overlay=private/relative
    if not python.is_file() or not (overlay/'pyrosetta/database').is_dir() or not (overlay/'rdkit').is_dir():
        raise ValueError('Preserved docking runtime is incomplete; install the private dependencies first')
    environment=os.environ.copy()
    environment['PYTHONPATH']=str(overlay)+(os.pathsep+environment['PYTHONPATH'] if environment.get('PYTHONPATH') else '')
    return [str(python),*arguments],environment

if __name__=='__main__':
    if len(sys.argv)<2:
        raise SystemExit('Usage: python3 docking_python.py <Python arguments or script>')
    root=Path(__file__).resolve().parents[3]
    try:
        command,environment=runtime_command(root,sys.argv[1:])
    except (OSError,ValueError,KeyError) as error:
        raise SystemExit('Local docking runtime unavailable: '+str(error))
    os.execve(command[0],command,environment)
