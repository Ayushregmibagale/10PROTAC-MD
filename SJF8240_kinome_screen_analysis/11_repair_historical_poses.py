"""Repair-stage entry point; see SCRIPT_GUIDE.md before running."""
from pathlib import Path
import runpy
if __name__ == "__main__":
 runpy.run_path(str(Path(__file__).resolve().parent / "repair_v3/scripts" / 'repair_legacy_poses.py'), run_name="__main__")
