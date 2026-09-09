"""Descriptive all-scene plot, no uncertainty bar or population claim."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
DELEGATED = ROOT.parents[1]/'delegated'
new = json.loads((DELEGATED/'a3_discrepancy_weighting/fits.json').read_text())
old = json.loads((DELEGATED/'a3_matched_frequency/fits.json').read_text())
assert len(new)==16 and all(r.get('evaluation_status')=='ok' for r in new)
labels = ['Low only', 'Low repeat', 'Raw warm', 'Isotropic', 'Rank 1', 'Rank 2']
keys = ['low', 'low_repeat', 'warm_raw', 'isotropic', 'rank1', 'rank2']
plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
fig, axs = plt.subplots(2, 2, figsize=(11,7), layout='constrained')
for col, ref in enumerate((False,True)):
    for seed, color, marker in ((8101,'#176087','o'),(8102,'#bd572b','s')):
        selected = []
        for key in keys:
            source = old if key in ('low','low_repeat') else new
            selected.append(next(r for r in source if r['seed']==seed and
                r['use_reference']==ref and r.get('choice',r.get('method'))==key))
        phase = [1000*r['sensor_low_band']['weighted_phase_rmse_rad'] for r in selected]
        material = [100*r['material_relative_error'] for r in selected]
        for ax, values in ((axs[0,col],phase),(axs[1,col],material)):
            ax.plot(np.arange(6),values,color=color,marker=marker,linewidth=1,label=f'Scene {seed}')
            ax.set_xticks(np.arange(6),labels,rotation=25,ha='right')
            ax.grid(axis='y',alpha=.2)
    for ax in axs[:,col]:
        maximum = max(float(np.max(line.get_ydata())) for line in ax.lines)
        ax.set_ylim(0, maximum*1.12)
    axs[0,col].set_title('With electronics reference' if ref else 'No electronics reference')
    axs[0,col].set_ylabel('Held-out low-band sensor phase (mrad)')
    axs[1,col].set_ylabel('Material relative error (%)')
    axs[0,col].legend(frameon=False)
fig.suptitle('Frozen coarse/fine discrepancy weighting: both development scenes\n'
             'Same low + high data for all four continuation methods; lower error is better',fontsize=13)
path = ROOT/'figures/discrepancy_weighting.png'
fig.savefig(path,dpi=170)
print(path)
