"""Parent audit of matched acquisitions and all twelve fixed outcomes."""
import hashlib
import json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent
WORK=OUT.parents[1]/'delegated/a3_matched_frequency'
rows=json.loads((WORK/'fits.json').read_text())
assert len(rows)==12
assert len({(r['seed'],r['choice'],r['use_reference']) for r in rows})==12
source_hash=hashlib.sha256((WORK/'matched_frequency.py').read_bytes()).hexdigest()
assert all(r['implementation_sha256']==source_hash for r in rows)
assert json.loads((WORK/'checks.json').read_text())['passed']
assert json.loads((WORK/'reference_data_checks.json').read_text())['passed']
for seed in [8101,8102]:
    data=np.load(WORK/f'data_{seed}.npz')
    assert np.array_equal(data['low'],data['low_high'][:3])
    assert np.array_equal(data['low'],data['low_repeat'][:3])
    assert not np.array_equal(data['low'][1],data['low_repeat'][3])
    assert len({r['sharing']['reference_sha256'] for r in rows if r['seed']==seed})==1
    assert len({r['sharing']['low_observation_sha256'] for r in rows if r['seed']==seed})==1
summary=[]
for ref in [False,True]:
    for choice in ['low','low_high','low_repeat']:
        group=[r for r in rows if r['use_reference']==ref and r['choice']==choice]
        entry=dict(use_reference=ref,choice=choice,fit_count=2,
            optimizer_statuses=[r.get('optimizer_status') for r in group],
            evaluation_statuses=[r.get('evaluation_status') for r in group])
        if all(r.get('evaluation_status')=='ok' for r in group):
            entry.update(pose_mm=float(np.mean([r['pose_error_m'] for r in group]))*1000,
                material_percent=float(np.mean([r['material_relative_error'] for r in group]))*100,
                sensor_phase_rad=float(np.mean([r['sensor_low_band']['weighted_phase_rmse_rad'] for r in group])),
                sensor_field_percent=float(np.mean([r['sensor_low_band']['relative_field_error'] for r in group]))*100,
                structural_phase_rad=float(np.mean([r['structural_low_band']['weighted_phase_rmse_rad'] for r in group])),
                structural_field_percent=float(np.mean([r['structural_low_band']['relative_field_error'] for r in group]))*100,
                mean_solve_seconds=float(np.mean([r['solve_seconds'] for r in group])))
        summary.append(entry)
result=dict(data_sharing_and_code_integrity_passed=True,summary=summary,
    source_sha256=source_hash,scope='two development scenes; original low data retained; equal extra count high versus repeat; independent-code N64 ADDA versus N32 inverse')
(OUT/'results/matched_frequency_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
