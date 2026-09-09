"""Registered shared-N32-model frequency attribution control.

Default checks use synthetic means and saved noise; no N32 generation or fits.
Shared-model agreement is intentional and is not independent Maxwell validation.
"""
from __future__ import annotations

import argparse
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

PROTOCOL = mf.A3/'CORRECT_MODEL_FREQUENCY_PROTOCOL.md'


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance():
    files = [Path(__file__), PROTOCOL, MATCHED/'matched_frequency.py', MATCHED/'fits.json',
             *[MATCHED/f'data_{seed}.npz' for seed in (8101, 8102)],
             *[mf.A3/name for name in ('nonspherical3d.py', 'nonspherical_calibration.py', 'maxwell3d.py')],
             HERE.parent/'a3_maxwell_refine/tangent_fft.py', HERE.parent/'a3_maxwell_fft/maxwell_fft.py']
    return {str(p): file_hash(p) for p in files}


def read_npz(path):
    with np.load(path) as archive:
        arrays = {k: archive[k].copy() for k in archive.files}
    for value in arrays.values():
        value.setflags(write=False)
    return arrays


def atomic_npz(path, arrays):
    temporary = path.with_name(path.name+f'.{os.getpid()}.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def assemble(saved, training_mean, held_sensor, held_structural):
    """Only substitute field means; preserve saved sigma, noises and reference."""
    assert training_mean.shape == (4, 12, 3, 4)
    assert held_sensor.shape == held_structural.shape == (3, 17, 3, 4)
    sigma = float(saved['sigma'])
    low = training_mean[:3] + sigma*saved['low_noise']
    high = training_mean[3] + sigma*saved['extra_noise']
    repeat = training_mean[1] + sigma*saved['extra_noise']
    arrays = dict(true=saved['true'], sigma=saved['sigma'], low_noise=saved['low_noise'],
                  extra_noise=saved['extra_noise'], reference=saved['reference'],
                  training_receivers=saved['training_receivers'], held_receivers=saved['held_receivers'],
                  training_mean=training_mean, held_sensor_mean=held_sensor,
                  held_structural_mean=held_structural, low=low,
                  low_high=np.concatenate((low, high[None])),
                  low_repeat=np.concatenate((low, repeat[None])))
    for value in arrays.values():
        value.setflags(write=False)
    return arrays


def adapter_data(arrays):
    return dict(sigma=float(arrays['sigma']),
                observations={name: arrays[name] for name in mf.CHOICES},
                references={name: arrays['reference'] for name in mf.CHOICES},
                reference=arrays['reference'], extra_noise=arrays['extra_noise'])


def assert_intervention(saved, arrays):
    for name in ('true', 'sigma', 'low_noise', 'extra_noise', 'reference',
                 'training_receivers', 'held_receivers'):
        assert np.array_equal(saved[name], arrays[name]), name
    sigma = float(saved['sigma'])
    mean = arrays['training_mean']
    assert np.array_equal(arrays['low'], mean[:3]+sigma*saved['low_noise'])
    assert np.array_equal(arrays['low_high'][3], mean[3]+sigma*saved['extra_noise'])
    assert np.array_equal(arrays['low_repeat'][3], mean[1]+sigma*saved['extra_noise'])
    sharing = mf.assert_shared(adapter_data(arrays))
    return dict(**sharing, low_noise_sha256=mf.digest(arrays['low_noise']),
                saved_sigma_exact=True, saved_noise_exact=True,
                only_field_means_substituted=True)


def generate_means(saved):
    """N32 physics is called only by --run, never by checks()."""
    start = time.perf_counter()
    truth = saved['true']
    training_model = mf.ExplicitFrequencyModel(32, mf.CHOICES['low_high'], saved['training_receivers'])
    training_mean, _ = training_model.field_jac(truth)
    training_work = training_model.ledger()
    del training_model
    gc.collect()
    held_model = mf.ExplicitFrequencyModel(32, mf.LOW, saved['held_receivers'])
    held_sensor, _ = held_model.field_jac(truth)
    held_structural = held_sensor/mf.electronics(truth, mf.LOW)[:, None, None, :]
    held_work = held_model.ledger()
    del held_model
    gc.collect()
    return assemble(saved, training_mean, held_sensor, held_structural), dict(
        generation_seconds=time.perf_counter()-start,
        training_generation_work=training_work, held_generation_work=held_work)


def scene_data(seed, saved, hashes):
    path = HERE/f'data_{seed}.npz'
    record_path = HERE/f'data_{seed}.json'
    previous = json.loads(record_path.read_text()) if record_path.exists() else None
    if previous:
        assert previous['provenance'] == hashes, 'Input/source hash mismatch'
        if previous['status'] == 'ok':
            assert file_hash(path) == previous['data_sha256'], 'Generated data cache changed'
            arrays = read_npz(path)
            assert_intervention(saved, arrays)
            return previous, arrays
        if previous['status'] == 'failed':
            return previous, None
        assert previous['status'] == 'started', 'Unknown data-generation checkpoint'
    record = dict(seed=seed, status='started', provenance=hashes,
                  scope='Shared N32 generator/inverse; intentional correct-model control')
    if previous:
        record['attempt_history'] = mf.preserve_interruption(previous, 'restart_interrupted_mean_generation')
    mf.write_json(record_path, record)
    start = time.perf_counter()
    arrays = None
    try:
        arrays, timing = generate_means(saved)
        record.update(timing, intervention_checks=assert_intervention(saved, arrays))
        atomic_npz(path, arrays)
        record.update(status='ok', data_sha256=file_hash(path))
    except Exception as exc:
        record.update(status='failed', failure_stage='mean_generation',
                      error_type=type(exc).__name__, error=str(exc))
    record['generation_attempt_wall_seconds'] = time.perf_counter()-start
    mf.write_json(record_path, record)
    return record, arrays


def load_prior():
    rows = json.loads((MATCHED/'fits.json').read_text())
    index = {(r['seed'], r['choice'], r['use_reference']): r for r in rows}
    assert len(index) == len(rows) == 12
    assert all(r.get('evaluation_status') == 'ok' for r in rows)
    return index


def checks():
    start = time.perf_counter()
    hashes = provenance()
    prior = load_prior()
    records = []
    for seed in (8101, 8102):
        saved = read_npz(MATCHED/f'data_{seed}.npz')
        assert np.array_equal(saved['training_receivers'], mf.receivers())
        assert np.array_equal(saved['held_receivers'], mf.receivers(17, 1.6))
        for choice in mf.CHOICES:
            for use_ref in (False, True):
                original = prior[seed, choice, use_ref]
                assert original['sharing']['reference_sha256'] == mf.digest(saved['reference'])
                assert original['sharing']['sigma'] == float(saved['sigma'])
                assert original['start'] == mf.BASE.tolist()
        # Deterministic nonphysical shape fixtures; no solver or new RNG draws.
        t = np.arange(4*12*3*4).reshape(4,12,3,4)
        raw = 1e-3*(np.cos(t)+1j*np.sin(t))
        training = raw*mf.electronics(saved['true'], mf.CHOICES['low_high'])[:,None,None,:]
        h = np.arange(3*17*3*4).reshape(3,17,3,4)
        structural = 1e-3*(np.cos(h)+1j*np.sin(h))
        factor = mf.electronics(saved['true'], mf.LOW)[:,None,None,:]
        sensor = structural*factor
        arrays = assemble(saved, training, sensor, sensor/factor)
        record = dict(seed=seed, **assert_intervention(saved, arrays))
        electronic_error = float(np.linalg.norm(arrays['held_structural_mean']-structural)/np.linalg.norm(structural))
        assert electronic_error < 1e-14
        record['structural_electronics_roundtrip_error'] = electronic_error
        records.append(record)
    assert provenance() == hashes
    result = dict(passed=True, scenes=records, provenance=hashes,
                  checkpoint_recovery=mf.checkpoint_checks(),
                  physical_solver_calls=0, generated_N32_means=False, benchmark_fits=0,
                  seconds=time.perf_counter()-start,
                  scope='Fixture assembly/convention and frozen-source/data checks; existing validated physics imported unchanged')
    mf.write_json(HERE/'checks.json', result)
    print(json.dumps(dict(passed=True, physical_solver_calls=0,
          implementation_sha256=hashes[str(Path(__file__))], seconds=result['seconds'])), flush=True)


def fit(row, arrays, save):
    model = None
    start = time.perf_counter()
    row['current_stage'] = 'fit_setup'
    save()
    try:
        model = mf.ExplicitFrequencyModel(32, mf.CHOICES[row['choice']], arrays['training_receivers'])
        row['setup_seconds'] = time.perf_counter()-start
        row['current_stage'] = 'fit_solve'
        save()
        data = adapter_data(arrays)
        cache = {}
        def ev(z):
            if cache.get('key') != z.tobytes():
                r, j = mf.residual_jac(z, model, data, row['choice'], row['use_reference'])
                cache.update(key=z.tobytes(), r=r, j=j)
            return cache['r'], cache['j']
        solve_start = time.perf_counter()
        try:
            opt = least_squares(lambda z: ev(z)[0], mf.BASE.copy(), jac=lambda z: ev(z)[1],
                bounds=(mf.LO,mf.HI), x_scale=mf.SCALE, max_nfev=35,
                ftol=1e-9, xtol=1e-9, gtol=1e-7)
        finally:
            row['solve_seconds'] = time.perf_counter()-solve_start
        row.update(estimated=opt.x.tolist(), status='converged' if opt.status>0 else 'iteration_limit',
                   optimizer_status=int(opt.status), optimizer_message=str(opt.message),
                   nfev=int(opt.nfev), njev=int(opt.njev), objective=float(opt.cost),
                   optimality=float(opt.optimality), active_mask=opt.active_mask.tolist(),
                   reduced_chisquare=float(opt.fun@opt.fun/(len(opt.fun)-13)),
                   training_work=model.ledger(), evaluation_status='pending')
        save()  # Freeze before all offline metrics.
    finally:
        row['fit_attempt_wall_seconds'] = time.perf_counter()-start
        if model is not None:
            row['training_work'] = model.ledger()
        del model
        gc.collect()


def evaluate(row, arrays, save):
    row['current_stage'] = 'parameter_evaluation'
    row['evaluation_status'] = 'started'
    save()
    z = np.array(row['estimated'])
    truth = arrays['true']
    row.update(material_relative_error=float(abs(z[0]-truth[0])/truth[0]),
               pose_error_m=float(np.linalg.norm(z[1:4]-truth[1:4])),
               delay_error_m=float(abs(z[4]-truth[4])),
               low_band_electronic_response=mf.metrics(mf.electronics(z,mf.LOW),mf.electronics(truth,mf.LOW)))
    row['current_stage'] = 'low_band_evaluation'
    save()
    start = time.perf_counter()
    try:
        mf.complete_evaluation(row, arrays['held_sensor_mean'], arrays['held_structural_mean'])
    finally:
        row['evaluation_attempt_wall_seconds'] = time.perf_counter()-start
    row['current_stage'] = 'completed'
    save()


def compare(rows, prior):
    output = []
    for row in rows:
        baselines = [('independent_matched_same_choice', prior[row['seed'],row['choice'],row['use_reference']])]
        if row['choice'] == 'low_high':
            baselines += [('correct_model_'+r['choice'],r) for r in rows
                          if r['seed']==row['seed'] and r['use_reference']==row['use_reference'] and r['choice'] in ('low','low_repeat')]
        for label, base in baselines:
            item = dict(seed=row['seed'], choice=row['choice'], use_reference=row['use_reference'],
                        baseline=label, direction='correct-model row minus baseline; negative error difference is improvement')
            if row.get('evaluation_status') == base.get('evaluation_status') == 'ok':
                item['differences'] = {field:{metric:row[field][metric]-base[field][metric]
                    for metric in ('relative_field_error','weighted_phase_rmse_rad')}
                    for field in ('sensor_low_band','structural_low_band')}
                item.update(pose_error_difference_m=row['pose_error_m']-base['pose_error_m'],
                            material_relative_error_difference=row['material_relative_error']-base['material_relative_error'])
            else:
                item['status'] = 'not_comparable_due_to_failure'
            output.append(item)
    mf.write_json(HERE/'comparisons.json',output)


def benchmark():
    hashes = provenance()
    manifest = HERE/'manifest.json'
    if manifest.exists() and json.loads(manifest.read_text()) != hashes:
        raise RuntimeError('Source/protocol/frozen inputs changed; preserve existing experiment')
    mf.write_json(manifest,hashes)
    checks()
    prior = load_prior()
    path = HERE/'fits.json'
    rows = json.loads(path.read_text()) if path.exists() else []
    assert all(r['provenance']==hashes for r in rows), 'Checkpoint hash mismatch'
    for seed in (8101,8102):
        saved = read_npz(MATCHED/f'data_{seed}.npz')
        generation, arrays = scene_data(seed,saved,hashes)
        for choice in mf.CHOICES:
            for use_ref in (False,True):
                matches = [r for r in rows if (r['seed'],r['choice'],r['use_reference'])==(seed,choice,use_ref)]
                assert len(matches)<=1, 'Duplicate fit checkpoint'
                previous = matches[0] if matches else None
                action = mf.checkpoint_action(previous)
                if action == 'skip_completed':
                    continue
                if action == 'resume_evaluation':
                    row = previous
                    row['attempt_history'] = mf.preserve_interruption(row,action)
                else:
                    row = dict(seed=seed, choice=choice, use_reference=use_ref,
                               provenance=hashes, status='started', current_stage='mean_generation',
                               start=mf.BASE.tolist(), max_nfev=35, n=32,
                               frequencies=mf.CHOICES[choice],
                               scope='Intentional shared-N32-model control; not independent Maxwell/continuum validation')
                    if previous:
                        row['attempt_history'] = mf.preserve_interruption(previous,action)
                        rows[rows.index(previous)] = row
                    else:
                        rows.append(row)
                def save():
                    mf.write_json(path,rows)
                save()
                start = time.perf_counter()
                try:
                    if generation['status'] != 'ok':
                        raise RuntimeError('Shared mean generation failed; retained in data record')
                    if 'data_sha256' in row:
                        assert row['data_sha256']==generation['data_sha256'], 'Frozen data changed'
                    row.update(data_file=f'data_{seed}.npz',data_sha256=generation['data_sha256'],
                               generation_seconds_shared=generation['generation_seconds'],
                               generation_cost_shared_across_six_fits=True,
                               sigma=float(arrays['sigma']),intervention_checks=assert_intervention(saved,arrays))
                    save()
                    if action != 'resume_evaluation':
                        fit(row,arrays,save)
                    evaluate(row,arrays,save)
                except Exception as exc:
                    row.update(failure_stage=row['current_stage'],error_type=type(exc).__name__,error=str(exc))
                    if 'estimated' in row:
                        row['evaluation_status'] = 'failed'
                    else:
                        row.update(status='failed',evaluation_status='not_reached')
                row['latest_attempt_wall_seconds'] = time.perf_counter()-start
                if row.get('attempt_history'):
                    row['uncheckpointed_interrupted_cost'] = 'unknown, not counted as zero'
                save()
                print(json.dumps({k:row.get(k) for k in ('seed','choice','use_reference','status','evaluation_status','failure_stage')}),flush=True)
                gc.collect()
    compare(rows,prior)
    assert provenance()==hashes, 'Frozen source/input mutation detected'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run',action='store_true',help='Generate N32 means and run12fits after parent authorization')
    args = parser.parse_args()
    if args.run:
        # Exclude the other two experiment runners as well as duplicate local
        # runs. r+b obtains their existing locks without creating/editing files.
        with (MATCHED/'benchmark.lock').open('r+b') as matched_lock, \
             (HERE.parent/'a3_discrepancy_weighting/benchmark.lock').open('r+b') as discrepancy_lock, \
             (HERE/'benchmark.lock').open('a') as own_lock:
            for lock in (matched_lock,discrepancy_lock,own_lock):
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            benchmark()
    else:
        checks()
