"""Read-only integrity/invariant audit of the shipped A4 evidence.

This does not validate continuum accuracy or the statistical error model.
It never reruns an experiment or changes a raw result.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent

def load(name):
    return json.loads((ROOT / name).read_text())

def main() -> int:
    manifest = load('provenance/frozen_manifest.json')
    checks = {}
    for name, expected in manifest['source_hashes'].items():
        checks['hash:' + name] = hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    branch = load('results/frozen_branch.json')['rows']
    cal = load('results/frozen_calibration.json')['rows']
    dda = load('results/frozen_dda.json')['rows']
    checks['counts'] = len(branch) == 144 and len(cal) == 72 and len(dda) == 96
    checks['branch_scene_units'] = len({r['scene_id'] for r in branch}) == 24
    checks['calibration_scene_units'] = len({r['scene'] for r in cal}) == 8
    checks['branch_flags'] = all(
        r['wrong_accepted'] == (r['accepted'] and not r['correct']) and
        r['correct_rejected'] == (r['correct'] and not r['accepted']) and
        r['wrong_selection_when_covered'] == (r['final_bank_covered'] and not r['correct'])
        for r in branch)
    checks['cell_acceptance_never_hides_unresolved'] = all(
        not r['accepted'] or r['coverage']['unresolved_cells'] == 0
        for r in branch if r['policy'] == 'coverage_aware')
    checks['acquisition_counts'] = all(
        r['training_complex_data'] == 27 + r['new_complex_acquisitions'] and
        r['validation_complex_data'] == 27 + r['new_complex_acquisitions'] and
        r['new_complex_acquisitions'] in (0, 9) for r in branch)
    checks['loss_matches_raw_parameters'] = all(np.isclose(
        r['scaled_task_loss'],
        (np.linalg.norm(np.array(r['estimate'][:3]) - r['true'][:3]) / .015)**2 +
        ((r['estimate'][3] - r['true'][3]) / .15)**2,
        rtol=1e-10, atol=1e-10) for r in cal)
    checks['charged_cost_sum'] = all(np.isclose(
        r['total_seconds'], sum(r[k] for k in ('fit_seconds', 'pilot_seconds', 'mode_seconds', 'screen_seconds')),
        atol=1e-10) for r in cal)
    fine = {r['scene']: r for r in cal if r['method'] == 'full_fine'}
    controllers = [r for r in cal if r['method'] == 'risk_controller']
    checks['negative_controller_result'] = all(
        r['selected_method'] == 'full_fine' and r['estimate'] == fine[r['scene']]['estimate'] and
        r['total_seconds'] > fine[r['scene']]['total_seconds'] for r in controllers)
    original_audit = load('results/A4_AUDIT.json')
    for name in ('frozen_branch_data.npz', 'frozen_calibration_data.npz'):
        p = ROOT / 'results' / name
        checks['raw_hash:' + name] = hashlib.sha256(p.read_bytes()).hexdigest() == original_audit['raw_npz'][name]['sha256']
        with np.load(p, allow_pickle=False) as arrays:
            checks['finite:' + name] = all(np.isfinite(arrays[k]).all() for k in arrays.files)
    correction = load('results/reference_v2_fresh.json')
    checks['fresh_reference_sample'] = len(correction) == 16
    checks['known_bug_preserved_and_disclosed'] = (
        'J[:-6' in (ROOT/'src/controls.py').read_text() and
        'J[:-6' in (ROOT/'docs/PROVENANCE_AND_CORRECTIONS.md').read_text())
    print(json.dumps({'all_passed': all(checks.values()), 'checks': checks,
                      'scope': 'file integrity and stored finite-model invariants, not physical certification'}, indent=2))
    return 0 if all(checks.values()) else 1

if __name__ == '__main__':
    sys.exit(main())
