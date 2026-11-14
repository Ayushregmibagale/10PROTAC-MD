"""Kinome chemical-engagement review; see SCRIPT_GUIDE.md."""
from pathlib import Path
import runpy
if __name__ == "__main__":
 runpy.run_path(str(Path(__file__).resolve().parent / "engagement_v2/scripts" / 'write_report.py'), run_name="__main__")
