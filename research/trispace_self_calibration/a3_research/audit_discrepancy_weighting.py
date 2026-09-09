"""Parent frozen-input, endpoint, and paired-mechanism audit; requires all 16 rows."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
WORK = ROOT.parents[1] / 'delegated/a3_discrepancy_weighting'
MATCHED = WORK.parent / 'a3_matched_frequency'
rows = json.loads((WORK/'fits.json').read_text())
methods = ('warm_raw', 'isotropic', 'rank1', 'rank2')
assert len(rows) == 16
assert len({(r['seed'], r['use_reference'], r['method']) for r in rows}) == 16
assert all(r.get('evaluation_status') in ('ok', 'failed', 'not_reached') for r in rows)
manifest = json.loads((WORK/'manifest.json').read_text())
for path, digest in manifest.items():
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest, path
assert all(r['provenance'] == manifest for r in rows)
assert json.loads((WORK/'checks.json').read_text())['passed']
prior = {(r['seed'], r['use_reference'], r['choice']): r
         for r in json.loads((MATCHED/'fits.json').read_text())}
mode_audits = []
for seed in (8101, 8102):
    for ref in (False, True):
        mode = json.loads((WORK/f'mode_{seed}_{int(ref)}.json').read_text())
        assert mode['provenance'] == manifest
        if mode['status'] != 'ok':
            mode_audits.append(dict(seed=seed, use_reference=ref, status=mode['status']))
            continue
        d, di = np.array(mode['delta']), np.array(mode['delta_i'])
        assert d.shape == di.shape == (288,)
        assert np.array_equal(mode['pilot'], prior[seed, ref, 'low']['estimated'])
        assert np.allclose(di, np.r_[-d[144:], d[:144]], rtol=1e-12, atol=1e-12)
        assert np.isclose(d@d, mode['trace_inflation'])
        mode_audits.append(dict(seed=seed, use_reference=ref, status='ok',
            trace_inflation=float(d@d), isotropic_precision=1/(1+float(d@d)/len(d)),
            rank1_mode_precision=1/(1+float(d@d)),
            construction_seconds=mode['construction_seconds']))
summary = []
for ref in (False, True):
    for method in methods:
        group = [r for r in rows if r['use_reference']==ref and r['method']==method]
        good = [r for r in group if r.get('evaluation_status')=='ok']
        item = dict(use_reference=ref, method=method, count=len(group),
                    evaluated_count=len(good), optimizer_statuses=[r.get('optimizer_status') for r in group])
        # Do not silently drop failures from means.
        if len(good) == len(group):
            item.update(pose_mm=float(np.mean([r['pose_error_m'] for r in group]))*1000,
                material_percent=float(np.mean([r['material_relative_error'] for r in group]))*100,
                sensor_phase_rad=float(np.mean([r['sensor_low_band']['weighted_phase_rmse_rad'] for r in group])),
                sensor_field_percent=float(np.mean([r['sensor_low_band']['relative_field_error'] for r in group]))*100,
                structural_field_percent=float(np.mean([r['structural_low_band']['relative_field_error'] for r in group]))*100,
                high_raw_objective=float(np.mean([r['high_white_objective'] for r in group])),
                historical_pilot_plus_new_cost_sum=float(np.mean([
                    r['pilot_setup_seconds_historical']+r['pilot_solve_seconds_historical']+
                    (r['shared_mode_construction_seconds'] if method!='warm_raw' else 0)+
                    r['continuation_attempt_wall_seconds'] for r in group])))
        summary.append(item)
paired = []
for seed in (8101, 8102):
    for ref in (False, True):
        idx = {r['method']: r for r in rows if r['seed']==seed and r['use_reference']==ref}
        raw, cold = idx['warm_raw'], prior[seed,ref,'low_high']
        item = dict(seed=seed, use_reference=ref)
        if raw.get('evaluation_status')=='ok':
            item['warm_minus_cold_raw_objective'] = raw['objective']-cold['objective']
        a,b = idx['rank1'], idx['isotropic']
        if a.get('evaluation_status')==b.get('evaluation_status')=='ok':
            item['rank1_minus_isotropic'] = dict(
                pose_mm=1000*(a['pose_error_m']-b['pose_error_m']),
                material_percent=100*(a['material_relative_error']-b['material_relative_error']),
                sensor_phase_rad=a['sensor_low_band']['weighted_phase_rmse_rad']-b['sensor_low_band']['weighted_phase_rmse_rad'])
        paired.append(item)
result = dict(integrity_passed=True, modes=mode_audits, summary=summary, paired=paired,
              failures=[{k:r.get(k) for k in ('seed','use_reference','method','status','evaluation_status','failure_stage','error')}
                        for r in rows if r.get('evaluation_status')!='ok'],
              scope='Two previously inspected development scenes, no population inference or novelty claim; cost sum is not fresh end-to-end timing.')
(ROOT/'results/discrepancy_weighting_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
