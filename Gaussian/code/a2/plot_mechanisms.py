"""Plot executed mechanisms by physical setting, not opaque case indices."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'figures/a2';OUT.mkdir(parents=True,exist_ok=True)
rows=json.loads((ROOT/'runs/a2/certificate/results.json').read_text())['rows']
fig,axes=plt.subplots(1,3,figsize=(13,4),sharey=True,layout='constrained')
for ax,f in zip(axes,sorted({r['frequency_GHz'] for r in rows})):
    for material,color,label in zip(range(3),('C0','C1','C2'),('Gaussian','Overlapping Gaussians','Sharp rectangle')):
        rr=[r for r in rows if r['frequency_GHz']==f and r['material']==material]
        ax.semilogy([r['R'] for r in rr],[r['corrected_rel'] for r in rr],'-o',color=color,label=label)
        ax.semilogy([r['R'] for r in rr],[r['bound_rel'] for r in rr],'--',color=color)
    ax.set_title(f'{f:g} GHz');ax.set_xlabel('Primal current modes (R)');ax.grid(alpha=.2)
axes[0].set_ylabel('Relative scattered-field error / bound');axes[1].legend(fontsize=8)
fig.suptitle('Lossy discrete models: solid = actual error; dashed = residual bound\nMore current modes improve accuracy; the bound remains conservative',fontsize=11)
fig.savefig(OUT/'certificate_by_setting.png',dpi=170);plt.close(fig)
p=json.loads((ROOT/'runs/a2/ports/quadrature_ports_results.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for r in p['local_t_response_port_sweep']:
    axes[0].semilogy([s['requested_ports_per_component'] for s in r['ports']],
                     [s['scattered_field_relative_error'] for s in r['ports']],'-o',label=r['configuration'].replace('_',' '))
axes[0].set_xlabel('Response ports per Gaussian component');axes[0].set_ylabel('Relative scattered-field error')
axes[0].set_title('Three components; exact local T before compression');axes[0].legend(fontsize=8)
c=p['cylinder_analytic_control']
axes[1].loglog([r['n'] for r in c],[r['scattered_relative_error_to_cylinder_series'] for r in c],'-o')
axes[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator());axes[1].set_xticks([32,64,128],['32','64','128']);axes[1].set_xlabel('Grid cells per axis')
axes[1].set_ylabel('Error against analytic cylinder series');axes[1].set_title('Independent scalar forward check')
for ax in axes:ax.grid(alpha=.2)
fig.savefig(OUT/'ports_and_grid.png',dpi=170);plt.close(fig)
