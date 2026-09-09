"""Parent audit of the completed shared-model attribution control."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
WORK=ROOT.parents[1]/'delegated/a3_correct_model_frequency'
OLD=WORK.parent/'a3_matched_frequency'
WEIGHT=WORK.parent/'a3_discrepancy_weighting'
rows=json.loads((WORK/'fits.json').read_text())
assert len(rows)==12 and len({(r['seed'],r['choice'],r['use_reference']) for r in rows})==12
assert all(r.get('evaluation_status') in ('ok','failed','not_reached') for r in rows)
manifest=json.loads((WORK/'manifest.json').read_text())
assert all(r['provenance']==manifest for r in rows)
for path,digest in manifest.items():
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
assert json.loads((WORK/'checks.json').read_text())['passed']
alignment=[]
for seed in (8101,8102):
    saved=np.load(OLD/f'data_{seed}.npz')
    generated=np.load(WORK/f'data_{seed}.npz')
    meta=json.loads((WORK/f'data_{seed}.json').read_text())
    assert meta['status']=='ok' and meta['data_sha256']==hashlib.sha256((WORK/f'data_{seed}.npz').read_bytes()).hexdigest()
    for key in ('true','sigma','low_noise','extra_noise','reference','training_receivers','held_receivers'):
        assert np.array_equal(saved[key],generated[key]),key
    sigma=float(saved['sigma'])
    assert np.array_equal(generated['low'],generated['training_mean'][:3]+sigma*saved['low_noise'])
    assert np.array_equal(generated['low_high'][:3],generated['low'])
    assert np.array_equal(generated['low_repeat'][:3],generated['low'])
    assert np.array_equal(generated['low_high'][3],generated['training_mean'][3]+sigma*saved['extra_noise'])
    assert np.array_equal(generated['low_repeat'][3],generated['training_mean'][1]+sigma*saved['extra_noise'])
    # Post-hoc oracle diagnostic only. It does NOT enter any estimator or weight.
    independent_mean=saved['low_high'][3]-sigma*saved['extra_noise']
    error=generated['training_mean'][3]-independent_mean
    e=np.sqrt(2)/sigma*np.r_[error.real.ravel(),error.imag.ravel()]
    for ref in (False,True):
        mode=json.loads((WEIGHT/f'mode_{seed}_{int(ref)}.json').read_text())
        d,di=np.array(mode['delta']),np.array(mode['delta_i'])
        cosine=float(e@d/(np.linalg.norm(e)*np.linalg.norm(d)))
        second=float(e@di/(np.linalg.norm(e)*np.linalg.norm(di)))
        alignment.append(dict(seed=seed,use_reference=ref,
            true_numerical_error_white_norm=float(np.linalg.norm(e)),
            pilot_mode_white_norm=float(np.linalg.norm(d)), signed_cosine=cosine,
            rank1_error_energy_fraction=cosine**2,
            rank2_error_energy_fraction=cosine**2+second**2,
            scope='Post-hoc diagnostic compares low-pilot mode with numerical discrepancy at truth, not at pilot; truth was not used in weights.'))
summary=[]
for ref in (False,True):
    for choice in ('low','low_high','low_repeat'):
        group=[r for r in rows if r['use_reference']==ref and r['choice']==choice]
        record=dict(use_reference=ref,choice=choice,count=len(group),
            statuses=[r.get('optimizer_status') for r in group],
            evaluation_statuses=[r.get('evaluation_status') for r in group])
        if all(r.get('evaluation_status')=='ok' for r in group):
            record.update(pose_mm=1000*float(np.mean([r['pose_error_m'] for r in group])),
                material_percent=100*float(np.mean([r['material_relative_error'] for r in group])),
                sensor_phase_rad=float(np.mean([r['sensor_low_band']['weighted_phase_rmse_rad'] for r in group])),
                structural_field_percent=100*float(np.mean([r['structural_low_band']['relative_field_error'] for r in group])))
        summary.append(record)
result=dict(integrity_and_intervention_passed=True,summary=summary,
            posthoc_oracle_mode_alignment=alignment,
            scope='Two inspected shared-model controls; neither independent validation nor expectation-level information adjudication.')
(ROOT/'results/correct_model_frequency_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
