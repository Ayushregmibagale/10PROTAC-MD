"""Step03: extract lysine features. Run --help for options."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/"_workflow"))
from launcher import main
if __name__=="__main__":main(3)
