"""Numbered entry point; see SCRIPT_GUIDE.md."""
from pathlib import Path
import runpy,sys,ast
if __name__ == "__main__":
 if not __debug__:raise SystemExit("Run without Python optimization; integrity checks are required.")
 script=Path(__file__).resolve().parent/"scripts"/'05_independent_verification.py'
 if "--help" in sys.argv and '05_independent_verification.py' not in ["00_fetch_references.py","02_extract_md_compatibility.py"]:
  print(ast.get_docstring(ast.parse(script.read_text())));raise SystemExit(0)
 runpy.run_path(str(script),run_name="__main__")
