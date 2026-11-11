from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[1];dest=ROOT/'results/figures';dest.mkdir(exist_ok=True)
df=pd.read_csv(ROOT/'results/pose_stability_and_bridge_exposure.tsv',sep='\t');colors=['#bb6b45' if x==0 else '#367f9b' for x in df.label];fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained')
for ax,col,sd,title,xlabel in [(axes[0],'NZ_sasa_pose_mean_A2','NZ_sasa_pose_SD_A2','Lysine solvent exposure','Mean lysine NZ SASA (Å²)'),(axes[1],'joint_fraction_pose_mean','joint_fraction_pose_SD','Joint geometric occupancy','Exposed + reachable + bridge fraction')]:
 for i,row in df.iterrows():ax.errorbar(row[col],i,xerr=row[sd],fmt='o',color=colors[i],capsize=3)
 ax.set(yticks=np.arange(len(df)),yticklabels=df.target,xlabel=xlabel,title=title);ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
fig.legend(handles=[Patch(color='#bb6b45',label='Non-degraded (2 targets)'),Patch(color='#367f9b',label='Degraded (6 targets)')],loc='outside lower center',ncol=2);fig.suptitle('Existing MD: overlapping per-lysine geometric summaries\nBars show pose SD, not uncertainty across independent proteins',fontsize=12);fig.savefig(dest/'lysine_geometry.png',dpi=180);fig.savefig(dest/'lysine_geometry.svg');plt.close(fig)
p=ROOT/'results/pair_scores.tsv'
if p.exists():
 data=pd.read_csv(p,sep='\t');fig,axes=plt.subplots(1,2,figsize=(10,6),layout='constrained')
 for ax,role in zip(axes,['elob','vhl']):
  d=data[data.role==role];grid=d.pivot(index='pair',columns='stage',values='margin')[['baseline','temporal']];bound=max(abs(grid.to_numpy()).max(),1e-6);im=ax.imshow(grid,cmap='RdBu',vmin=-bound,vmax=bound,aspect='auto');ax.set(xticks=[0,1],xticklabels=['Original','Enriched'],yticks=np.arange(12),yticklabels=grid.index,title=role.upper());fig.colorbar(im,ax=ax,label='Held-out degrader score − non-degrader score')
 fig.suptitle('Matched held-out scores: blue = correctly ordered pair');fig.savefig(dest/'pair_margins.png',dpi=180);fig.savefig(dest/'pair_margins.svg');plt.close(fig)
