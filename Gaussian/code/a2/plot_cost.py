"""Presentation-only regeneration from the completed v2 benchmark."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
x=json.loads((ROOT/'runs/a2/matrix_free_v2/benchmark.json').read_text())
fig,ax=plt.subplots(figsize=(8,4.6),layout='constrained')
labels=['dense/dual GN','plain CG','rank 12','rank 32','rank 64']
for r in x['results']:
    ys=[r['dense_GN']['wall_seconds']]+[r['methods'][k]['wall_seconds_including_range'] for k in ['plain_cg','data_som_rank_12','data_som_rank_32','data_som_rank_64']]
    ax.semilogy(labels,ys,'o-',label=f"p={r['parameters']}")
ax.set(ylabel='Seconds including Jacobian / range setup',title='Fixed geometry, amplitude-only local GN; four residuals')
ax.legend();ax.grid(alpha=.2);fig.savefig(ROOT/'figures/a2/local_cost.png',dpi=160)
