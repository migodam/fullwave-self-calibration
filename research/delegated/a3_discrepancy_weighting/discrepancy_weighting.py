"""Registered frozen discrepancy-mode continuation; default is bounded checks.

All prior acquisition data, estimates and physics are imported read-only.
No rank/strength tuning, ADDA truth, or measurement residual enters the mode.
"""
from __future__ import annotations

import argparse
import copy
import fcntl
import gc
import hashlib
import json
import os
from pathlib import Path
import sys
import time

for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
MATCHED = HERE.parent/'a3_matched_frequency'
sys.path.insert(0, str(MATCHED))
import matched_frequency as mf

METHODS = ('warm_raw', 'isotropic', 'rank1', 'rank2')
PROTOCOL = mf.A3/'DISCREPANCY_WEIGHTING_PROTOCOL.md'


def provenance():
    files = [Path(__file__), PROTOCOL, MATCHED/'matched_frequency.py', MATCHED/'fits.json',
             *[MATCHED/f'data_{seed}.npz' for seed in (8101, 8102)],
             *[mf.A3/name for name in ('nonspherical3d.py', 'nonspherical_calibration.py', 'maxwell3d.py')],
             HERE.parent/'a3_maxwell_refine/tangent_fft.py', HERE.parent/'a3_maxwell_fft/maxwell_fft.py']
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def load_inputs():
    rows = json.loads((MATCHED/'fits.json').read_text())
    index = {(r['seed'], r['use_reference'], r['choice']): r for r in rows}
    assert len(rows) == len(index) == 12
    assert all(r.get('evaluation_status') == 'ok' and 'estimated' in r for r in rows)
    data = {}
    for seed in (8101, 8102):
        with np.load(MATCHED/f'data_{seed}.npz') as saved:
            arrays = {k: saved[k].copy() for k in saved.files}
        for value in arrays.values():
            value.setflags(write=False)
        assert np.array_equal(arrays['low_high'][:3], arrays['low'])
        assert np.array_equal(arrays['low_repeat'][:3], arrays['low'])
        assert np.array_equal(arrays['training_receivers'], mf.receivers())
        assert np.array_equal(arrays['held_receivers'], mf.receivers(17, 1.6))
        for use_ref in (False, True):
            for choice in ('low', 'low_high', 'low_repeat'):
                prior = index[seed, use_ref, choice]
                assert prior['sharing']['low_observation_sha256'] == mf.digest(arrays['low'])
                assert prior['sharing']['reference_sha256'] == mf.digest(arrays['reference'])
                assert prior['sharing']['sigma'] == float(arrays['sigma'])
        data[seed] = arrays
    return index, data


def mode_vectors(d, sigma):
    # Difference convention is fine minus coarse; simultaneous sign reversal
    # would leave either rank covariance unchanged.
    return np.sqrt(2)/sigma*mf.real(d), np.sqrt(2)/sigma*mf.real(1j*d)


def whitening(delta, delta_i, method):
    """Dense symmetric inverse square root; high block only (n=288)."""
    delta = np.asarray(delta, float)
    delta_i = np.asarray(delta_i, float)
    n = len(delta)
    eye = np.eye(n)
    energy = float(delta@delta)
    if method == 'warm_raw':
        return eye, eye
    if method == 'isotropic':
        c = 1+energy/n
        return eye/np.sqrt(c), eye*c
    if method == 'rank1':
        e = delta[:, None]
    elif method == 'rank2':
        e = np.stack((delta, delta_i), axis=1)/np.sqrt(2)
    else:
        raise ValueError(method)
    u, s, _ = np.linalg.svd(e, full_matrices=False)
    w = eye + (u*((1+s*s)**(-.5)-1))@u.T
    return w, eye+e@e.T


def high_rows(block_complex_size):
    # matched residual layout is [Re(all four blocks); Im(all four blocks); ref].
    b = block_complex_size
    return np.r_[np.arange(3*b, 4*b), np.arange(7*b, 8*b)]


def weighted_residual_jac(z, model, arrays, use_reference, w):
    adapter_data = dict(sigma=float(arrays['sigma']),
                        observations={'low_high': arrays['low_high']},
                        references={'low_high': arrays['reference']})
    r, j = mf.residual_jac(z, model, adapter_data, 'low_high', use_reference)
    idx = high_rows(arrays['low_high'][3].size)
    r[idx] = w@r[idx]
    j[idx] = w@j[idx]
    return r, j


def construct_mode(pilot, sigma, frequencies=(18.,), ns=(16, 32)):
    """Only pilot, noise scale, and inverse models are accepted, never data/truth."""
    start = time.perf_counter()
    predictions = []
    work = {}
    for n in ns:
        model = mf.ExplicitFrequencyModel(n, frequencies)
        pred, _ = model.field_jac(np.array(pilot, float))
        predictions.append(pred[0])
        work[str(n)] = model.ledger()
        del model
        gc.collect()
    delta, delta_i = mode_vectors(predictions[1]-predictions[0], sigma)
    return dict(delta=delta.tolist(), delta_i=delta_i.tolist(),
                trace_inflation=float(delta@delta), dimension=len(delta),
                construction_seconds=time.perf_counter()-start, construction_work=work,
                construction_grids=list(ns), difference='N32 minus N16 (fine minus coarse)',
                pilot=list(pilot), sigma=float(sigma), strength=1., frozen=True)


def checks():
    start = time.perf_counter()
    initial_hashes = provenance()
    prior, inputs = load_inputs()
    rng = np.random.default_rng(93124)
    d = rng.normal(size=(12, 3, 4))+1j*rng.normal(size=(12, 3, 4))
    delta, delta_i = mode_vectors(d, 3.)
    n = len(delta)
    results = {}
    matrices = {}
    for method in METHODS:
        w, c = whitening(delta, delta_i, method)
        err = float(np.linalg.norm(w@c@w.T-np.eye(n))/np.sqrt(n))
        assert err < 1e-11
        assert np.max(abs(w-w.T)) < 1e-13
        results[method] = dict(whitening_identity_error=err,
                              trace_inflation=float(np.trace(c-np.eye(n))))
        matrices[method] = w
    traces = np.array([results[m]['trace_inflation'] for m in METHODS[1:]])
    assert np.allclose(traces, delta@delta, rtol=1e-13, atol=1e-12)
    assert abs(delta@delta_i) < 1e-12 and np.isclose(delta@delta, delta_i@delta_i)
    for method in METHODS:
        w0, c0 = whitening(np.zeros(n), np.zeros(n), method)
        assert np.array_equal(w0, np.eye(n)) and np.array_equal(c0, np.eye(n))
    # Tagged real and imaginary arrays catch accidentally whitening a contiguous
    # slice spanning low imaginary rows, or mixing the electronics reference.
    tagged = np.arange(4*12*3*4).reshape(4, 12, 3, 4)
    tagged = tagged + 1j*(10000+tagged)
    idx = high_rows(tagged[3].size)
    assert np.array_equal(mf.real(tagged)[idx], mf.real(tagged[3]))
    z = np.array([2.3, .03, -.02, .04, .03]+[.01]*4+[.02]*4)
    model = mf.ExplicitFrequencyModel(8, mf.CHOICES['low_high'])
    arrays = inputs[8101]
    fd_errors = []
    unchanged = []
    # All 13 columns of full weighted residuals, both reference conditions.
    h = 1e-5
    for use_ref in (False, True):
        raw_r, raw_j = weighted_residual_jac(z, model, arrays, use_ref, matrices['warm_raw'])
        keep = np.ones(len(raw_r), bool)
        keep[idx] = False
        for method in METHODS:
            w = matrices[method]
            r, j = weighted_residual_jac(z, model, arrays, use_ref, w)
            assert np.array_equal(r[keep], raw_r[keep])
            assert np.array_equal(j[keep], raw_j[keep])
            errors = []
            for column in range(13):
                step = np.eye(13)[column]*h
                fd = (weighted_residual_jac(z+step, model, arrays, use_ref, w)[0]-
                      weighted_residual_jac(z-step, model, arrays, use_ref, w)[0])/(2*h)
                errors.append(float(np.linalg.norm(fd-j[:, column])/np.linalg.norm(fd)))
            assert max(errors) < 1e-7
            fd_errors.append(dict(method=method, use_reference=use_ref, column_errors=errors))
            unchanged.append(dict(method=method, use_reference=use_ref, low_and_reference_bit_identical=True))
    assert provenance() == initial_hashes
    result = dict(passed=True, whitening=results,
                  equal_trace_relative_spread=float(np.ptp(traces)/(delta@delta)),
                  high_real_imag_ordering=True, unchanged=unchanged,
                  fixed_weight_full_residual_column_errors=fd_errors,
                  checkpoint_recovery=mf.checkpoint_checks(), frozen_inputs_unchanged=True,
                  input_checks=dict(scenes=[8101,8102], all_12_previous_endpoints_available=True),
                  seconds=time.perf_counter()-start, provenance=initial_hashes,
                  scope='N8 fixed-weight residual Jacobian and synthetic covariance checks only; no N16/N32 modes or continuation fits')
    mf.write_json(HERE/'checks.json', result)
    print(json.dumps(dict(passed=True, max_fd_error=max(max(r['column_errors']) for r in fd_errors),
          max_whitening_error=max(r['whitening_identity_error'] for r in results.values()),
          seconds=result['seconds'], implementation_sha256=initial_hashes[str(Path(__file__))])), flush=True)


def mode_record(seed, use_reference, pilot, sigma, hashes):
    path = HERE/f'mode_{seed}_{int(use_reference)}.json'
    old = json.loads(path.read_text()) if path.exists() else None
    if old is not None:
        assert old['provenance'] == hashes, 'Mode inputs/source changed; preserve existing run.'
        assert old['pilot_sha256'] == mf.digest(np.array(pilot['estimated'], dtype=float))
        if old['status'] == 'ok':
            assert old['delta_sha256'] == mf.digest(np.array(old['delta'], dtype=float))
            assert old['delta_i_sha256'] == mf.digest(np.array(old['delta_i'], dtype=float))
        if old['status'] in ('ok', 'failed'):
            return old
        if old['status'] != 'started':
            raise RuntimeError('Unknown mode checkpoint status')
    record = dict(seed=seed, use_reference=use_reference, status='started', provenance=hashes,
                  pilot=pilot['estimated'], pilot_sha256=mf.digest(np.array(pilot['estimated'], dtype=float)),
                  pilot_setup_seconds_historical=pilot['setup_seconds'],
                  pilot_solve_seconds_historical=pilot['solve_seconds'])
    if old:
        record['attempt_history'] = mf.preserve_interruption(old, 'restart_interrupted_mode_construction')
    mf.write_json(path, record)
    start = time.perf_counter()
    try:
        record.update(construct_mode(pilot['estimated'], sigma), status='ok')
        record['delta_sha256'] = mf.digest(np.array(record['delta'], dtype=float))
        record['delta_i_sha256'] = mf.digest(np.array(record['delta_i'], dtype=float))
    except Exception as exc:
        record.update(status='failed', failure_stage='mode_construction',
                      error_type=type(exc).__name__, error=str(exc))
    record['attempt_wall_seconds'] = time.perf_counter()-start
    mf.write_json(path, record)
    return record


def frozen_evaluation(row, arrays, w, save):
    """Save each successful audit stage; resume only missing stages."""
    z = np.array(row['estimated'])
    # Reconstruct inexpensive parameter metrics even if a previous process
    # stopped immediately after the first endpoint checkpoint.
    row['current_stage'] = 'parameter_evaluation'
    truth = arrays['true']
    row.update(material_relative_error=float(abs(z[0]-truth[0])/truth[0]),
               pose_error_m=float(np.linalg.norm(z[1:4]-truth[1:4])),
               delay_error_m=float(abs(z[4]-truth[4])),
               low_band_electronic_response=mf.metrics(mf.electronics(z, mf.LOW), mf.electronics(truth, mf.LOW)))
    save()
    if row.get('low_evaluation_status') != 'ok':
        row['current_stage'] = 'low_band_evaluation'
        save()
        start = time.perf_counter()
        try:
            audit = {}
            audit['estimated'] = row['estimated']
            mf.complete_evaluation(audit, arrays['held_sensor_mean'], arrays['held_structural_mean'])
            row.update({k: v for k, v in audit.items() if k not in ('estimated', 'evaluation_status')})
            row['low_evaluation_status'] = 'ok'
        finally:
            row['low_evaluation_attempt_seconds'] = time.perf_counter()-start
        save()
    if row.get('high_evaluation_status') != 'ok':
        row['current_stage'] = 'high_measurement_evaluation'
        save()
        start = time.perf_counter()
        model = None
        try:
            model = mf.ExplicitFrequencyModel(32, (18.,))
            high, _ = model.field_jac(z)
            residual = np.sqrt(2)/float(arrays['sigma'])*mf.real(high[0]-arrays['low_high'][3])
            row.update(high_measurement_fit=mf.metrics(high[0], arrays['low_high'][3]),
                       high_white_residual_norm=float(np.linalg.norm(residual)),
                       high_surrogate_residual_norm=float(np.linalg.norm(w@residual)),
                       high_white_objective=float(residual@residual/2),
                       high_surrogate_objective=float(np.linalg.norm(w@residual)**2/2),
                       high_evaluation_work=model.ledger(), high_evaluation_status='ok')
        finally:
            row['high_evaluation_attempt_seconds'] = time.perf_counter()-start
            if model is not None:
                row['high_evaluation_work'] = model.ledger()
        save()
    row['evaluation_status'] = 'ok'
    row['current_stage'] = 'completed'
    save()


def fit_one(row, arrays, use_reference, w, save):
    start = time.perf_counter()
    model = None
    try:
        row['current_stage'] = 'continuation_setup'
        save()
        model = mf.ExplicitFrequencyModel(32, mf.CHOICES['low_high'])
        row['continuation_setup_seconds'] = time.perf_counter()-start
        row['current_stage'] = 'continuation_fit'
        save()
        cache = {}
        def ev(z):
            if cache.get('key') != z.tobytes():
                r, j = weighted_residual_jac(z, model, arrays, use_reference, w)
                cache.update(key=z.tobytes(), r=r, j=j)
            return cache['r'], cache['j']
        solve_start = time.perf_counter()
        try:
            opt = least_squares(lambda z: ev(z)[0], np.array(row['pilot']), jac=lambda z: ev(z)[1],
                bounds=(mf.LO, mf.HI), x_scale=mf.SCALE, max_nfev=35,
                ftol=1e-9, xtol=1e-9, gtol=1e-7)
        finally:
            row['continuation_solve_seconds'] = time.perf_counter()-solve_start
        row.update(estimated=opt.x.tolist(), status='converged' if opt.status>0 else 'iteration_limit',
                   optimizer_status=int(opt.status), optimizer_message=str(opt.message),
                   nfev=int(opt.nfev), njev=int(opt.njev), objective=float(opt.cost),
                   optimality=float(opt.optimality), active_mask=opt.active_mask.tolist(),
                   training_work=model.ledger(), evaluation_status='pending')
        # First freeze and persist every endpoint, before using evaluation truth.
        save()
    finally:
        row['continuation_attempt_wall_seconds'] = time.perf_counter()-start
        if model is not None:
            row['training_work'] = model.ledger()
        del model
        gc.collect()


def comparisons(rows, prior):
    output = []
    for row in rows:
        seed, use_ref = row['seed'], row['use_reference']
        baselines = [(f'matched_{choice}', prior[seed, use_ref, choice])
                     for choice in ('low', 'low_repeat', 'low_high')]
        baselines += [(other['method'], other) for other in rows
                      if other['seed']==seed and other['use_reference']==use_ref and other['method']!=row['method']]
        for name, base in baselines:
            contrast = dict(seed=seed, use_reference=use_ref, method=row['method'], baseline=name,
                            signed_difference='method minus baseline; negative error difference is improvement')
            if row.get('evaluation_status') == base.get('evaluation_status') == 'ok':
                contrast['differences'] = {field: {metric: row[field][metric]-base[field][metric]
                    for metric in ('relative_field_error', 'weighted_phase_rmse_rad')}
                    for field in ('sensor_low_band', 'structural_low_band')}
                contrast['pose_error_difference_m'] = row['pose_error_m']-base['pose_error_m']
                contrast['material_relative_error_difference'] = row['material_relative_error']-base['material_relative_error']
                if row['method']=='warm_raw' and name=='matched_low_high':
                    contrast['raw_endpoint_scaled_distance'] = float(np.linalg.norm(
                        (np.array(row['estimated'])-np.array(base['estimated']))/mf.SCALE))
                    contrast['raw_objective_difference'] = row['objective']-base['objective']
                    contrast['inspection_required_if_endpoint_differs'] = True
            else:
                contrast['status'] = 'not_comparable_due_to_failure'
            output.append(contrast)
    mf.write_json(HERE/'comparisons.json', output)


def benchmark():
    hashes = provenance()
    manifest = HERE/'manifest.json'
    if manifest.exists() and json.loads(manifest.read_text()) != hashes:
        raise RuntimeError('Source, protocol or frozen inputs changed; preserve this experiment and choose a new directory.')
    mf.write_json(manifest, hashes)
    checks()
    prior, inputs = load_inputs()
    path = HERE/'fits.json'
    rows = json.loads(path.read_text()) if path.exists() else []
    if any(r['provenance'] != hashes for r in rows):
        raise RuntimeError('Checkpoint source/input mismatch')
    for seed in (8101, 8102):
        arrays = inputs[seed]
        for use_ref in (False, True):
            pilot = prior[seed, use_ref, 'low']
            for method in METHODS:
                matches = [r for r in rows if r['seed']==seed and r['use_reference']==use_ref and r['method']==method]
                if len(matches)>1:
                    raise RuntimeError('Duplicate checkpoint configuration')
                previous = matches[0] if matches else None
                action = mf.checkpoint_action(previous)
                if action == 'skip_completed':
                    continue
                if action == 'resume_evaluation':
                    row = previous
                    row['attempt_history'] = mf.preserve_interruption(row, action)
                else:
                    row = dict(seed=seed, use_reference=use_ref, method=method, provenance=hashes,
                               status='started', current_stage='mode_construction', pilot=pilot['estimated'],
                               pilot_setup_seconds_historical=pilot['setup_seconds'],
                               pilot_solve_seconds_historical=pilot['solve_seconds'],
                               strength=1., max_nfev=35,
                               time_accounting='Historical pilot plus new costs is a cost sum, not fresh end-to-end wall time',
                               scope='Two development scenes; frozen discrepancy proxy, not known covariance or error enclosure')
                    if previous:
                        row['attempt_history'] = mf.preserve_interruption(previous, action)
                        rows[rows.index(previous)] = row
                    else:
                        rows.append(row)
                def save():
                    mf.write_json(path, rows)
                save()
                attempt_start = time.perf_counter()
                try:
                    row['current_stage'] = 'mode_cache_validation' if action=='resume_evaluation' else 'mode_construction'
                    mode = mode_record(seed, use_ref, pilot, float(arrays['sigma']), hashes)
                    if mode['status'] != 'ok':
                        raise RuntimeError('Shared mode construction failed; preserved in mode record')
                    row['mode_file'] = f'mode_{seed}_{int(use_ref)}.json'
                    if 'mode_delta_sha256' in row:
                        assert row['mode_delta_sha256'] == mode['delta_sha256'], 'Frozen mode changed'
                    row['mode_delta_sha256'] = mode['delta_sha256']
                    row['shared_mode_construction_seconds'] = mode['construction_seconds']
                    row['mode_construction_cost_shared_across_four_methods'] = True
                    w, _ = whitening(mode['delta'], mode['delta_i'], method)
                    save()
                    if action != 'resume_evaluation':
                        fit_one(row, arrays, use_ref, w, save)
                    row['evaluation_status'] = 'started'
                    frozen_evaluation(row, arrays, w, save)
                except Exception as exc:
                    row.update(failure_stage=row['current_stage'], error_type=type(exc).__name__, error=str(exc))
                    if 'estimated' in row:
                        row['evaluation_status'] = 'failed'
                    else:
                        row.update(status='failed', evaluation_status='not_reached')
                row['latest_attempt_wall_seconds'] = time.perf_counter()-attempt_start
                if row.get('attempt_history'):
                    row['uncheckpointed_interrupted_cost'] = 'unknown; preserved history, not zero'
                save()
                print(json.dumps({k: row.get(k) for k in ('seed','use_reference','method','status','evaluation_status','failure_stage')}), flush=True)
                gc.collect()
    comparisons(rows, prior)
    assert provenance() == hashes, 'Frozen input mutation detected'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', action='store_true', help='Run registered 16 fits only after parent approval')
    args = parser.parse_args()
    if args.run:
        with (HERE/'benchmark.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            benchmark()
    else:
        checks()
