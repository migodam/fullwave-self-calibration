"""All-case descriptive stress plot; log axes expose scale differences."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
report=json.loads((ROOT/'results/shape_material_stress_audit.json').read_text())
rows=report['all_rows']
assert len(rows)==32 and all(r['evaluation_status']=='ok' for r in rows)
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
methods=('low','warm_raw','isotropic','rank1')
names=('Low only','Raw warm','Isotropic','Rank 1')
colors=('#777777','#b64b32','#2082a0','#5b4b96')
markers=('x','s','^','o')
for col,ref in enumerate((False,True)):
    for i,(method,name,color,marker) in enumerate(zip(methods,names,colors,markers)):
        group=[next(r for r in rows if r['seed']==seed and r['use_reference']==ref and r['method']==method)
               for seed in (8401,8402,8403,8404)]
        x=np.arange(4)+(i-1.5)*.12
        for ax,key in ((axes[0,col],'material_percent'),(axes[1,col],'pose_mm')):
            vals=[r[key] for r in group]
            assert min(vals)>0
            ax.scatter(x,vals,label=name,c=color,marker=marker,s=48)
            ax.set_yscale('log')
            ax.set_xticks(range(4),['Ellipsoid\neps 1.8','Ellipsoid\neps 3.5','Box\neps 1.8','Box\neps 3.5'])
            ax.grid(axis='y',which='both',alpha=.15)
            ax.spines[['right','top']].set_visible(False)
    axes[0,col].set_title('Reference present' if ref else 'Reference absent')
    axes[0,col].set_ylabel('Material error (%) — log scale')
    axes[1,col].set_ylabel('Receiver-position error (mm) — log scale')
    axes[0,col].legend(frameon=False,ncol=2,fontsize=9)
fig.suptitle('Shape/material stress test: all four registered development cases\n'
             'Known supports; one unknown material coefficient; lower error is better',fontsize=13)
path=ROOT/'figures/shape_material_stress.png'
fig.savefig(path,dpi=170)
print(path)
