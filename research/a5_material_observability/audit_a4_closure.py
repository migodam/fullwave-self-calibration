"""Audit frozen evidence and latest handoff without modifying old evidence."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
checks = {}
for filename in ('scale_gauge.json', 'interaction_ablation.json'):
    data = json.loads((HERE/'results'/filename).read_text())
    checks[filename+'_complete'] = data['complete']
    checks[filename+'_hashes'] = all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
                                    for path, digest in data['source_hashes'].items())
new = json.loads((HERE/'results/interaction_ablation.json').read_text())
checks['six_scales'] = [r['scale'] for r in new['rows']] == [.01,.03,.1,.3,1.,2.]
checks['numerical_checks'] = all(new['checks'].values())
checks['isolated_visible_all'] = all(r['methods']['isolated_sum']['relative_sensitivity'] > 0 for r in new['rows'])
checks['reported_direction_all'] = all(r['methods']['isolated_sum']['relative_sensitivity'] > r['methods']['interacting']['relative_sensitivity'] for r in new['rows'])
docs = [ROOT/'communication/A4_COMPLETE_RESEARCH_REPORT_ZH.md',
        ROOT/'Theory/Questions/Q5.md', HERE/'A4_CLOSURE_ADDENDUM_EN.md']
manifest = {}
for path in docs:
    content = path.read_text()
    checks[path.name+'_nonempty'] = len(content) > 1000
    checks[path.name+'_display_delimiters'] = sum(line.strip() == '$$' for line in content.splitlines()) % 2 == 0
    manifest[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
env = os.environ.copy()
env.update(PYTEST_DISABLE_PLUGIN_AUTOLOAD='1', PYTHONDONTWRITEBYTECODE='1',
           OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
runs = []
for args in (['-m','pytest','-q','-p','no:cacheprovider'], ['audit_package.py']):
    start = time.perf_counter()
    run = subprocess.run([sys.executable,*args], cwd=ROOT/'public_release/research/a4_reliability_v1',
                         env=env, capture_output=True, text=True, timeout=60)
    runs.append(dict(args=args, returncode=run.returncode, output=run.stdout+run.stderr,
                     wall_seconds=time.perf_counter()-start))
checks['a4_regression'] = runs[0]['returncode'] == 0 and '21 passed' in runs[0]['output']
checks['a4_package'] = runs[1]['returncode'] == 0
result = dict(passed=all(checks.values()), checks=checks, handoff_hashes=manifest, runs=runs,
              scope='Artifact/evidence/regression checks, not rendered math acceptance or TAP readiness.')
(HERE/'results/a4_closure_audit.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
assert result['passed']
