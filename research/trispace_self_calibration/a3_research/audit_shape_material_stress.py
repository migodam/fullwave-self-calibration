"""Require all registered stress-test outcomes before paired reporting."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
WORK=ROOT/'results/shape_material_stress'
rows=json.loads((WORK/'fits.json').read_text())
expected={(s,r,m) for s in (8401,8402,8403,8404) for r in (False,True)
          for m in ('low','warm_raw','isotropic','rank1')}
assert {(r['seed'],r['use_reference'],r['method']) for r in rows}==expected and len(rows)==32
assert all(r.get('evaluation_status') in ('ok','failed','not_reached') for r in rows)
manifest=json.loads((WORK/'manifest.json').read_text())
assert all(r['provenance']==manifest for r in rows)
for path,digest in manifest.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
assert json.loads((WORK/'checks.json').read_text())['passed']
summaries=[]
references=[]
paired=[]
for seed in (8401,8402,8403,8404):
    data=json.loads((WORK/f'data_{seed}.json').read_text())
    assert data['status']=='ok'
    assert hashlib.sha256((WORK/f'data_{seed}.npz').read_bytes()).hexdigest()==data['data_sha256']
    references.append({k:data[k] for k in ('seed','shape','epsilon','high_reference_relative_difference',
                       'held_high_reference_relative_difference','reference_sensitivity_below_one_percent')})
    for ref in (False,True):
        idx={r['method']:r for r in rows if r['seed']==seed and r['use_reference']==ref}
        for method,row in idx.items():
            assert row['data_sha256']==data['data_sha256']
            item={k:row.get(k) for k in ('seed','shape','epsilon','use_reference','method','status','evaluation_status','failure_stage')}
            if row.get('evaluation_status')=='ok':
                item.update(pose_mm=1000*row['pose_error_m'],material_percent=100*row['material_relative_error'],
                    sensor_phase_rad=row['sensor_low_band']['weighted_phase_rmse_rad'],
                    sensor_field_percent=100*row['sensor_low_band']['relative_field_error'],
                    structural_field_percent=100*row['structural_low_band']['relative_field_error'])
            summaries.append(item)
        a,b=idx['rank1'],idx['isotropic']
        comparison=dict(seed=seed,shape=a['shape'],epsilon=a['epsilon'],use_reference=ref,
                        direction='rank1 minus isotropic; negative error difference is improvement')
        if a.get('evaluation_status')==b.get('evaluation_status')=='ok':
            comparison.update(pose_mm=1000*(a['pose_error_m']-b['pose_error_m']),
                material_percent=100*(a['material_relative_error']-b['material_relative_error']),
                sensor_phase_rad=a['sensor_low_band']['weighted_phase_rmse_rad']-b['sensor_low_band']['weighted_phase_rmse_rad'],
                structural_field_percent=100*(a['structural_low_band']['relative_field_error']-b['structural_low_band']['relative_field_error']))
        else:comparison['status']='not_comparable_due_to_failure'
        paired.append(comparison)
result=dict(integrity_passed=True,reference_sensitivity=references,all_rows=summaries,
            primary_paired=paired,scope='Four predefined known-support/material development cases, not iid population or continuum proof.')
(ROOT/'results/shape_material_stress_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
