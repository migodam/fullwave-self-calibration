"""Summarize common-frequency frozen-estimate predictions, not pooled-band metrics."""
import json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent
rows=json.loads((OUT/'results/cross_band_audit.json').read_text())
assert len(rows)==8
summary=[];paired=[]
for method in ['unified','unified_reference']:
    for k in [3,9]:
        for band in ['low','high']:
            group=[r for r in rows if r['method']==method and (max(r['frequencies'])>9)==(band=='high')]
            assert len(group)==2
            values=[next(p for p in r['per_frequency'] if p['k']==k) for r in group]
            summary.append(dict(method=method,k=k,band=band,
                sensor_phase_mean=float(np.mean([p['sensor']['weighted_phase_rmse_rad'] for p in values])),
                sensor_field_mean=float(np.mean([p['sensor']['field_relative_error'] for p in values])),
                structural_phase_mean=float(np.mean([p['structural']['weighted_phase_rmse_rad'] for p in values])),
                structural_field_mean=float(np.mean([p['structural']['field_relative_error'] for p in values]))))
        for seed in [6101,6102]:
            get=lambda high:next(p for r in rows if r['seed']==seed and r['method']==method and (max(r['frequencies'])>9)==high for p in r['per_frequency'] if p['k']==k)
            lo,hi=get(False),get(True)
            paired.append(dict(seed=seed,method=method,k=k,
                sensor_phase_change=hi['sensor']['weighted_phase_rmse_rad']-lo['sensor']['weighted_phase_rmse_rad'],
                structural_field_change=hi['structural']['field_relative_error']-lo['structural']['field_relative_error']))
result=dict(scene_count=2,summary=summary,paired=paired,
    scope='same predicted k but training frequency sets/noise differ; frozen estimates; exploratory, not a pure high-frequency causal experiment')
(OUT/'results/cross_band_summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
