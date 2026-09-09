"""Independent fresh execution of the bounded pipeline package without its writes.

Verify expected-fail controls explicitly: a status string alone is not evidence.
"""
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
source=ROOT/'research/delegated/a3_algebra/regenerate_checks.py'
spec=importlib.util.spec_from_file_location('a3_pipeline_checks',source)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
cases=mod.build_cases()
by={c['case']:c for c in cases}
checks={
    'no_unexpected_failure':all(c['status'] in ('pass','expected_fail') for c in cases),
    'zero_denominator_flag_verified':by['2G_near_zero_negative_control']['metrics']['near_zero_flag_detected'] is True,
    'no_division_at_zero':by['2G_near_zero_negative_control']['metrics']['blind_division_attempted'] is False,
}
control=by['4C_small_residual_control']['metrics']
checks['small_residual_large_output_verified']=control['residual_norm']<=1e-5 and control['remaining_output_error']>=.99
payload={'status':'bounded algebra replication only','source':str(source),
         'cases':mod._clean_jsonable(cases),'parent_assertions':checks,'passed':all(checks.values())}
(HERE/'results/pipeline_replication.json').write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps({'case_count':len(cases),'parent_assertions':checks,'passed':payload['passed']}))
assert payload['passed']
