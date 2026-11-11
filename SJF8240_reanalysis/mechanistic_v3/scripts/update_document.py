from pathlib import Path
import json,datetime
from docx import Document
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'docs/STEP_BY_STEP.md';s=p.read_text().split('\n## Live results and progress')[0];s+='\n## Live results and progress\n\nUpdated '+datetime.datetime.now(datetime.timezone.utc).isoformat()+'.\n\n'
for stage in ['baseline','temporal']:
 for role in ['elob','vhl']:
  d=ROOT/'analysis'/stage/role;c=d/'completion.json';v=d/'independent_cpu_verification.json';n=len(list(d.glob('*.pt')))
  s+=f'- {stage}, {role}: {n}/12 saved folds. '
  if c.exists():
   r=json.loads(c.read_text());s+=f'Matched-pair AUROC {r["matched_pair_AUROC"]:.4f}; balanced accuracy {r["balanced_accuracy"]:.4f}; sensitivity {r["sensitivity"]:.4f}; specificity {r["specificity"]:.4f}. '
  s+=('Independent CPU replay passed.\n' if v.exists() else 'CPU verification pending.\n')
n=len(list((ROOT/'analysis/lysines/poses').glob('*/*/receipt.json')));s+=f'\nLysine extraction: {n}/141 cases saved. '
if (ROOT/'analysis/lysines/completion.json').exists():s+='Complete: 5,781 sampled frames.\n'
else:s+='In progress.\n'
if (ROOT/'RESULTS.md').exists():s+='\n'+(ROOT/'RESULTS.md').read_text()

if (ROOT/'validation/final_review.json').exists():
 s=s.replace('New three-stage investigation is IN PROGRESS','New three-stage investigation is COMPLETE').replace('— pending','— completed')
p.write_text(s);d=Document();d.core_properties.title='SJF8240: step-by-step analysis and verified results'
table=None
for line in s.splitlines():
 if line.startswith('|'):
  cells=[x.strip() for x in line.strip().strip('|').split('|')]
  if all(set(x)<=set('-: ') for x in cells):continue
  if table is None:table=d.add_table(rows=0, cols=len(cells));table.style='Table Grid'
  row=table.add_row().cells
  for cell,value in zip(row,cells):cell.text=value
  continue
 table=None
 if line.startswith('# '):d.add_heading(line[2:],0)
 elif line.startswith('## '):d.add_heading(line[3:],1)
 elif line.startswith('### '):d.add_heading(line[4:],2)
 elif line:d.add_paragraph(line)

from docx.shared import Inches
for name in ['lysine_geometry.png','pair_margins.png']:
 figure=ROOT/'results/figures'/name
 if figure.exists():d.add_picture(str(figure),width=Inches(6.4))
d.save(ROOT/'docs/STEP_BY_STEP.docx');print(p)
