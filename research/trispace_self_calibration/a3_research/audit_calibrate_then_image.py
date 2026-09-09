"""Audit all two-stage outcomes without treating oracle geometry as an algorithm."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/calibrate_then_image'
SOURCE=ROOT/'results/shape_material_stress'
rows=json.loads((OUT/'fits.json').read_text())
expected={(s,r,m) for s in (8401,8402,8403,8404) for r in (False,True)
          for m in ('warm_raw','isotropic','rank1','oracle_true_geometry')}
assert len(rows)==32 and {(r['seed'],r['use_reference'],r['geometry_source']) for r in rows}==expected
assert all(r.get('evaluation_status') in ('ok','failed','not_reached') for r in rows)
manifest=json.loads((OUT/'manifest.json').read_text())
for path,digest in manifest.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
assert all(r['provenance']==manifest for r in rows)
assert json.loads((OUT/'checks.json').read_text())['passed']
prior={(r['seed'],r['use_reference'],r['method']):r for r in json.loads((SOURCE/'fits.json').read_text())}
all_rows=[];paired=[]
for row in rows:
    seed,ref,source=row['seed'],row['use_reference'],row['geometry_source']
    low=prior[seed,ref,'low']
    base=np.array(low['estimated'])
    if source=='oracle_true_geometry':
        assert row['oracle']
        with np.load(SOURCE/f'data_{seed}.npz') as d:base[1:4]=d['true'][1:4]
    else:
        assert not row['oracle']
        base[1:4]=np.array(prior[seed,ref,source]['estimated'])[1:4]
    assert np.array_equal(base,row['start'])
    item={k:row.get(k) for k in ('seed','shape','epsilon','use_reference','geometry_source','oracle','status','evaluation_status','failure_stage')}
    if 'estimated' in row:assert np.array_equal(np.array(row['estimated'])[1:4],base[1:4])
    if row.get('evaluation_status')=='ok':
        item.update(pose_mm=1000*row['pose_error_m'],material_percent=100*row['material_relative_error'],
                    sensor_phase_rad=row['sensor_low_band']['weighted_phase_rmse_rad'],
                    structural_field_percent=100*row['structural_low_band']['relative_field_error'])
        if source=='rank1':
            direct=prior[seed,ref,'rank1']
            paired.append(dict(seed=seed,use_reference=ref,
                second_stage_material_percent=item['material_percent'],
                direct_rank1_material_percent=100*direct['material_relative_error'],
                low_only_material_percent=100*low['material_relative_error'],
                second_stage_structural_field_percent=item['structural_field_percent'],
                direct_rank1_structural_field_percent=100*direct['structural_low_band']['relative_field_error'],
                low_only_structural_field_percent=100*low['structural_low_band']['relative_field_error']))
    all_rows.append(item)
result=dict(integrity_and_fixed_geometry_passed=True,all_rows=all_rows,primary_paired=paired,
            scope='Four reused development cases; oracle is diagnostic, no calibrated confidence or general superiority.')
(ROOT/'results/calibrate_then_image_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
