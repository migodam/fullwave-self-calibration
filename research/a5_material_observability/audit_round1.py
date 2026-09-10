"""Persist independent A4 regression and round-one provenance checks."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
A4=ROOT/'public_release/research/a4_reliability_v1'
data=json.loads((HERE/'results/scale_gauge.json').read_text())
assert data['complete'] and len(data['rows'])==18 and len(data['modal_counterexamples'])==9
for path,digest in data['source_hashes'].items():
    assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
assert {(r['scale'],r['radius']) for r in data['rows']}=={
    (s,r) for s in (.01,.03,.1,.3,1.,2.) for r in (.2,.6,2.)}
assert max(r['gain_compensated_relative_error'] for r in data['modal_counterexamples'])<1e-12
assert all(r['linear_scaling_control']<1e-12 for r in data['rows'])
assert all(r['gains']['entry']['finite_pair_relative_residual']==0 for r in data['rows'])
env=os.environ.copy()
env.update(PYTEST_DISABLE_PLUGIN_AUTOLOAD='1',PYTHONDONTWRITEBYTECODE='1',
           OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
runs=[]
for args,name in [(['-m','pytest','-q','-p','no:cacheprovider'],'a4_regression'),
                  (['audit_package.py'],'a4_package')]:
    start=time.perf_counter()
    run=subprocess.run([sys.executable,*args],cwd=A4,env=env,capture_output=True,text=True,timeout=60)
    (HERE/'results'/f'{name}.txt').write_text(run.stdout+run.stderr)
    runs.append(dict(name=name,returncode=run.returncode,seconds=time.perf_counter()-start))
    assert run.returncode==0
assert '21 passed' in (HERE/'results/a4_regression.txt').read_text()
result=dict(passed=True,mechanism_cases=18,modal_counterexamples=9,runs=runs,
            source_hashes=data['source_hashes'],
            scope='Regression, stored-data and mechanism checks; not final material recovery or scientific novelty.')
(HERE/'results/round1_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
