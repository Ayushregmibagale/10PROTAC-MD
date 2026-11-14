"""Numbered entry point;see SCRIPT_GUIDE.md."""
from pathlib import Path
import runpy,sys,ast
if __name__=="__main__":
 if not __debug__:raise SystemExit("Run without optimization;checks are required.")
 script=Path(__file__).resolve().parent/"scripts"/'08_write_complete_report.py'
 if "--help" in sys.argv and '08_write_complete_report.py' not in ["02_extract_md_native_references.py","12_run_chemical_tail_sensitivity.py"]:
  print(ast.get_docstring(ast.parse(script.read_text())));raise SystemExit(0)
 runpy.run_path(str(script),run_name="__main__")
