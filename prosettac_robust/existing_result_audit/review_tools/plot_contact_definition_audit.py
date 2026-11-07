"""Export descriptive contact-definition/cutoff figures from existing MD audits."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0,str(ROOT/'workflow/scripts/lib'))
from protac_qc import sha256,atomic_json
run=ROOT/'existing_result_audit/contact_core_v1';source=run/'per_target_definitions.tsv';rows=list(csv.DictReader(source.open(),delimiter='\t'));out=run/'figures';out.mkdir(exist_ok=False)
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(8.4,4.2));fig.subplots_adjust(left=.075,right=.985,bottom=.28,top=.8,wspace=.3)
pois=sorted({x['poi'] for x in rows});x=np.arange(len(pois));primary={r['poi']:r for r in rows if float(r['cutoff_A'])==5.}
keys=['legacy_named_groups_all_atom_fraction','true_VHL_POI_all_atom_fraction','true_VHL_POI_heavy_atom_fraction'];labels=['Elongin B + POI, all atoms','VHL + POI, all atoms','VHL + POI, heavy atoms'];colors=['#999999','#56B4E9','#0072B2']
for i,(key,label,color) in enumerate(zip(keys,labels,colors)):axes[0].bar(x+(i-1)*.24,[100*float(primary[p][key]) for p in pois],width=.23,label=label,color=color)
axes[0].set_xticks(x,[p.removesuffix('_full') for p in pois],rotation=45,ha='right');axes[0].set_ylim(0,102);axes[0].set_ylabel('Mean proximity occupancy (%)');axes[0].set_title('Contact definitions at 5 Å',loc='left',fontsize=9);axes[0].legend(frameon=False,loc='lower left',bbox_to_anchor=(0,1.08),fontsize=7)
palette=['#000000','#E69F00','#56B4E9','#009E73','#D55E00','#0072B2','#CC79A7','#7A6A00']
for poi,color in zip(pois,palette):
 group=sorted([r for r in rows if r['poi']==poi],key=lambda r:float(r['cutoff_A']));axes[1].plot([float(r['cutoff_A']) for r in group],[100*float(r['true_VHL_POI_heavy_atom_fraction']) for r in group],marker='o',ms=2.8,lw=1.15,color=color,label=poi.removesuffix('_full'))
axes[1].axvline(5,color='#777777',ls='--',lw=.8,zorder=0);axes[1].set_ylim(0,102);axes[1].set_xlim(3.4,6.1);axes[1].set_xlabel('Distance cutoff (Å)');axes[1].set_ylabel('Mean proximity occupancy (%)');axes[1].set_title('True VHL/POI, heavy atoms',loc='left',fontsize=9);axes[1].legend(frameon=False,ncol=4,loc='upper center',bbox_to_anchor=(.5,-.24),fontsize=7,columnspacing=.8)
for letter,ax in zip(['a','b'],axes):ax.text(-.15,1.05,letter,transform=ax.transAxes,fontweight='bold',fontsize=10);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig.text(.075,.035,'141 existing poses; 201 sampled frames per pose, 0–20 ns; equal pose weights. Descriptive occupancies, without uncertainty intervals.',fontsize=7)
for ext in ['pdf','svg','png']:fig.savefig(out/('contact_definition_sensitivity.'+ext),dpi=300,facecolor='white')
plt.close(fig)
atomic_json(out/'provenance.json',{'status':'DESCRIPTIVE_REVIEW_FIGURE','input_sha256':sha256(source),'aggregate_summary_sha256':sha256(run/'aggregate_summary.json'),'source_sha256':sha256(Path(__file__)),'outputs_sha256':{p.name:sha256(p) for p in out.iterdir() if p.is_file()},'original_results_modified':False,'new_MD_executed':False})
print(str(out))
