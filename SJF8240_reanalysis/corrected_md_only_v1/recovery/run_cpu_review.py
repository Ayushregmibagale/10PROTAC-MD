"""Run frozen CPU verifier without the Python3.10 GPU package overlay."""
from pathlib import Path
import subprocess,os,sys
root=Path(__file__).resolve().parents[1]
environment=os.environ.copy();environment.pop('PYTHONPATH',None)
result=subprocess.run(['/home/ayushregmibagale/miniconda3/bin/python3',str(root/'scripts/verify_and_report.py')],env=environment)
sys.exit(result.returncode)
