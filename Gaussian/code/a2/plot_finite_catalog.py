"""Render decision-count summary from an already completed frozen catalog run."""
import json
from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];x=json.loads((ROOT/'runs/a2/finite_catalog/results.json').read_text());stages=x['stages'];fig,ax=plt.subplots(figsize=(7,3.5))
for case in x['cases']:
    vals=[sum(s['stages'][i]['certificate_decision']=='supported_new_structure' for s in case['seeds']) for i in range(len(stages))]
    ax.plot(stages,vals,'o-',label=case['truth'])
ax.set(xlabel='GMRES restart, one cycle',ylabel='certificate supported-new count / 10',ylim=(-.2,10.2),title='Frozen finite-catalog decision counts');ax.legend();fig.tight_layout();fig.savefig(ROOT/'figures/a2/catalog_decision_counts.png',dpi=160)
