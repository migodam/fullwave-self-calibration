"""Registered 12-fit acquisition control; existing A3 physics imported read-only.

Default command runs bounded checks. Full benchmark requires explicit --run.
No global KS mutation; repeated acquisitions share physics, never observations.
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

# Set before importing NumPy/SciPy (also set these in the invoking environment).
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
A3 = HERE.parents[1] / "trispace_self_calibration" / "a3_research"
sys.path.insert(0, str(A3))
from nonspherical3d import ShapeVIE, reference, independent_radiation
from nonspherical_calibration import metrics
from maxwell3d import receivers, dipole_kernel

LOW = (3., 6., 9.)
CHOICES = {"low": LOW, "low_high": LOW + (18.,), "low_repeat": LOW + (6.,)}
LOSS = .03
REF_SIGMA = .02
LO = np.array([1.2] + [-.25]*3 + [-.2] + [-.5]*4 + [-1.]*4)
HI = np.array([5.] + [.25]*3 + [.2] + [.5]*4 + [1.]*4)
SCALE = np.array([1.] + [.1]*12)
BASE = np.array([2.] + [0.]*12)


def write_json(path, value):
    # An interrupted replacement leaves the previous complete checkpoint intact.
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    with temporary.open('w') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def digest(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def real(a):
    return np.r_[a.real.ravel(), a.imag.ravel()]


def electronics(z, frequencies):
    k = np.asarray(frequencies)
    return np.exp(z[5:9][None, :] + 1j*z[9:13][None, :] + 1j*k[:, None]*z[4])


def electronics_jac(z, frequencies=LOW):
    y = electronics(z, frequencies)
    j = np.zeros(y.shape + (13,), complex)
    j[..., 4] = 1j*np.asarray(frequencies)[:, None]*y
    for t in range(4):
        j[:, t, 5+t] = y[:, t]
        j[:, t, 9+t] = 1j*y[:, t]
    return y, j


class ExplicitFrequencyModel:
    """Same 13 parameters as Model; one physical model per unique explicit k."""

    def __init__(self, n, frequencies, rx=None):
        self.frequencies = tuple(float(k) for k in frequencies)
        self.rx = receivers() if rx is None else np.asarray(rx)
        self.models = {}
        self.setup_seconds = {}
        for k in dict.fromkeys(self.frequencies):
            start = time.perf_counter()
            self.models[k] = ShapeVIE('ellipsoid', n, k)
            self.setup_seconds[str(k)] = time.perf_counter() - start
        self.cache_material = None
        self.currents = {}
        self.calls = 0
        self.field_seconds = {str(k): 0. for k in self.models}

    def field_jac(self, z):
        self.calls += 1
        if self.cache_material != z[0]:
            self.currents = {k: m.currents([z[0]+1j*LOSS], True)
                             for k, m in self.models.items()}
            self.cache_material = z[0]
        outputs = {}
        for k, m in self.models.items():
            start = time.perf_counter()
            p, dp = self.currents[k]
            rx = self.rx + z[1:4]
            op = dipole_kernel(rx, m.points, k)
            raw = (op @ p).reshape(len(rx), 3, 4)
            j = np.zeros(raw.shape + (13,), complex)
            j[..., 0] = (op @ dp[:, :, 0]).reshape(raw.shape)
            for d in range(3):
                h = np.eye(3)[d]*1e-5
                dop = (dipole_kernel(rx+h, m.points, k)-dipole_kernel(rx-h, m.points, k))/(2e-5)
                j[..., 1+d] = (dop @ p).reshape(raw.shape)
            factor = electronics(z, (k,))[0]
            y = raw*factor
            j *= factor[None, None, :, None]
            j[..., 4] = 1j*k*y
            for t in range(4):
                j[:, :, t, 5+t] = y[:, :, t]
                j[:, :, t, 9+t] = 1j*y[:, :, t]
            outputs[k] = (y, j)
            self.field_seconds[str(k)] += time.perf_counter() - start
        return (np.stack([outputs[k][0] for k in self.frequencies]),
                np.stack([outputs[k][1] for k in self.frequencies]))

    def ledger(self):
        return dict(calls=self.calls, per_frequency={str(k): dict(
            setup_seconds=self.setup_seconds[str(k)],
            field_projection_seconds=self.field_seconds[str(k)],
            work=copy.deepcopy(m.work), fft_work=copy.deepcopy(m.fft_work),
            max_true_relative_residual=max(m.true_relative_residuals, default=0.),
            physical_frequency=k, acquisition_multiplicity=self.frequencies.count(k))
            for k, m in self.models.items()})


def scene(seed):
    rng = np.random.default_rng(seed)
    z = np.zeros(13)
    z[0] = 2.5
    z[1:4] = rng.normal(size=3)
    z[1:4] *= .09 / np.linalg.norm(z[1:4])
    z[4] = .06
    z[5:9] = rng.normal(0, .06, 4)
    z[9:13] = rng.normal(0, .12, 4)
    return z, rng


def cached_sources():
    """Preflight all cache records before calling existing validated reader."""
    for k in LOW + (18.,):
        folder = A3 / 'results/adda' / f'ellipsoid_n64_k{k:g}_2.5_0.03_ldr'
        required = [folder/'shape.dat']
        for group in range(2):
            required.append(folder/f'prop{group}_status.json')
            required.extend(folder/f'prop{group}'/name for name in
                            ('log', 'DipPol-X', 'DipPol-Y', 'IncBeam-X', 'IncBeam-Y'))
        missing = [str(p) for p in required if not p.is_file()]
        if missing:
            raise FileNotFoundError('Read-only N64 reference cache incomplete: ' + ', '.join(missing))
    start = time.perf_counter()
    records = {k: reference('ellipsoid', 64, k, 2.5+1j*LOSS) for k in LOW+(18.,)}
    return records, time.perf_counter()-start


def reference_means(sources, z, rx, frequencies):
    raw = np.stack([independent_radiation(*sources[k][:2], rx+z[1:4], k)
                    for k in frequencies])
    return raw*electronics(z, frequencies)[:, None, None, :], raw


def make_data(seed, low_mean, high_mean):
    """Generate low once; separate extra draw reused between k18 and repeat k6."""
    z, rng = scene(seed)
    sigma = float(np.linalg.norm(low_mean)/np.sqrt(low_mean.size)*10**(-30/20))
    low_noise = (rng.normal(size=low_mean.shape)+1j*rng.normal(size=low_mean.shape))/np.sqrt(2)
    extra_noise = (rng.normal(size=high_mean.shape)+1j*rng.normal(size=high_mean.shape))/np.sqrt(2)
    low_y = low_mean + sigma*low_noise
    high_y = high_mean + sigma*extra_noise
    repeat_y = low_mean[1] + sigma*extra_noise
    ref_rng = np.random.default_rng(seed+10000)
    ref_mean = electronics(z, LOW)
    ref_noise = (ref_rng.normal(size=ref_mean.shape)+1j*ref_rng.normal(size=ref_mean.shape))/np.sqrt(2)
    ref = ref_mean + REF_SIGMA*ref_noise
    observations = {'low': low_y, 'low_high': np.concatenate((low_y, high_y[None])),
                    'low_repeat': np.concatenate((low_y, repeat_y[None]))}
    references = {name: ref for name in CHOICES}
    for a in [*observations.values(), ref, low_noise, extra_noise]:
        a.setflags(write=False)
    return dict(true=z, sigma=sigma, observations=observations, references=references,
                low_noise=low_noise, extra_noise=extra_noise, reference=ref)


def assert_shared(data):
    obs = data['observations']
    for name in CHOICES:
        assert np.array_equal(obs[name][:3], obs['low'])
        assert data['references'][name] is data['reference']
        assert np.array_equal(data['references'][name], data['reference'])
    assert not np.array_equal(obs['low_repeat'][3], obs['low'][1])
    return dict(low_observation_sha256=digest(obs['low']),
                reference_sha256=digest(data['reference']),
                extra_standardized_noise_sha256=digest(data['extra_noise']),
                sigma=data['sigma'], low_exactly_shared=True,
                reference_exactly_shared=True, repeat_is_new_measurement=True)


def residual_jac(z, model, data, choice, use_reference):
    pred, j = model.field_jac(z)
    r = np.sqrt(2)/data['sigma']*real(pred-data['observations'][choice])
    jj = j.reshape(-1, 13)
    jr = np.sqrt(2)/data['sigma']*np.r_[jj.real, jj.imag]
    if use_reference:
        rp, rj = electronics_jac(z, LOW)
        r = np.r_[r, np.sqrt(2)/REF_SIGMA*real(rp-data['references'][choice])]
        rj = rj.reshape(-1, 13)
        jr = np.r_[jr, np.sqrt(2)/REF_SIGMA*np.r_[rj.real, rj.imag]]
    return r, jr


def checks():
    start = time.perf_counter()
    import nonspherical_calibration as original
    original_ks = original.KS.copy()
    z = np.array([2.3, .03, -.02, .04, .03]+[.01]*4+[.02]*4)
    model = ExplicitFrequencyModel(8, (3., 6., 9., 18., 6.))
    y, j = model.field_jac(z)
    assert len(model.models) == 4 and np.array_equal(y[1], y[4])
    assert np.array_equal(j[1], j[4])
    fd = []
    h = 1e-5
    for d in range(13):
        step = np.eye(13)[d]*h
        fd.append((model.field_jac(z+step)[0]-model.field_jac(z-step)[0])/(2*h))
    fd = np.stack(fd, -1)
    errors = [float(np.linalg.norm(fd[..., d]-j[..., d])/np.linalg.norm(fd[..., d])) for d in range(13)]
    assert max(errors) < 1e-7
    rp, rj = electronics_jac(z)
    ref_errors = []
    for d in range(13):
        step = np.eye(13)[d]*h
        f = (electronics(z+step, LOW)-electronics(z-step, LOW))/(2*h)
        ref_errors.append(float(np.linalg.norm(f-rj[..., d])/max(np.linalg.norm(f), 1.)))
    assert max(ref_errors) < 1e-7
    sharing = []
    for seed in (8101, 8102):
        data = make_data(seed, y[:3], y[3])
        sharing.append(dict(seed=seed, **assert_shared(data)))
        # Audit the common standardized extra noise via actual generated differences.
        extra_high = (data['observations']['low_high'][3]-y[3])/data['sigma']
        extra_repeat = (data['observations']['low_repeat'][3]-y[1])/data['sigma']
        assert np.allclose(extra_high, data['extra_noise'], rtol=0, atol=1e-13)
        assert np.allclose(extra_repeat, data['extra_noise'], rtol=0, atol=1e-13)
    # Full real/imag residual and reference layout directional derivative.
    m = ExplicitFrequencyModel(8, CHOICES['low_repeat'])
    data = make_data(8101, y[:3], y[3])
    direction = np.random.default_rng(92001).normal(size=13)*SCALE
    residual_errors = []
    for use_ref in (False, True):
        r, jr = residual_jac(z, m, data, 'low_repeat', use_ref)
        f = (residual_jac(z+h*direction, m, data, 'low_repeat', use_ref)[0]-
             residual_jac(z-h*direction, m, data, 'low_repeat', use_ref)[0])/(2*h)
        error = float(np.linalg.norm(f-jr@direction)/np.linalg.norm(f))
        assert error < 1e-7
        assert len(r) == 2*4*12*3*4+(24 if use_ref else 0)
        residual_errors.append(error)
    # Original model equivalence on its original band, without altering its KS.
    old = original.Model(8)
    old_y, old_j = old.field_jac(z)
    eq = ExplicitFrequencyModel(8, tuple(original.KS))
    eq_y, eq_j = eq.field_jac(z)
    equivalence = float(np.linalg.norm(eq_y-old_y)/np.linalg.norm(old_y))
    jac_equivalence = float(np.linalg.norm(eq_j-old_j)/np.linalg.norm(old_j))
    assert equivalence < 1e-13 and jac_equivalence < 1e-13
    assert np.array_equal(original.KS, original_ks)
    checkpoint_results = checkpoint_checks()
    result = dict(passed=True, n=8, jacobian_column_errors=errors,
                  electronics_jacobian_errors=ref_errors,
                  real_residual_directional_errors=residual_errors,
                  original_model_field_error=equivalence,
                  original_model_jacobian_error=jac_equivalence,
                  global_KS_unchanged=True, sharing=sharing,
                  checkpoint_recovery=checkpoint_results,
                  repeated_physics_reused=True, seconds=time.perf_counter()-start,
                  scope='Small N8 derivative/assembly controls with synthetic means; no benchmark fits',
                  implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    write_json(HERE/'checks.json', result)
    print(json.dumps(result), flush=True)


def evaluate_frozen(z):
    start = time.perf_counter()
    held = receivers(17, 1.6)
    evaluator = ExplicitFrequencyModel(32, LOW, rx=held)
    pred, _ = evaluator.field_jac(z)
    result = dict(prediction=pred, structural=pred/electronics(z, LOW)[:, None, None, :],
                  evaluation_seconds=time.perf_counter()-start, evaluation_work=evaluator.ledger())
    return result


def checkpoint_action(row):
    """Classify an existing attempt after the exclusive run lock is acquired."""
    if row is None:
        return 'new'
    if row.get('evaluation_status') in ('ok', 'failed', 'not_reached') or row.get('status') == 'failed':
        return 'skip_completed'
    if 'estimated' in row:
        if row.get('evaluation_status') not in (None, 'pending', 'started'):
            raise RuntimeError('Unrecognized evaluation checkpoint state; preserve and inspect it.')
        return 'resume_evaluation'
    if row.get('status') == 'started' and 'fit_and_evaluation_seconds' not in row:
        return 'restart_interrupted_fit'
    raise RuntimeError('Unrecognized pre-estimate checkpoint state; preserve and inspect it.')


def preserve_interruption(row, action):
    # Store exact prior attempt metadata, excluding recursively nested history.
    previous = copy.deepcopy({k: v for k, v in row.items() if k != 'attempt_history'})
    history = copy.deepcopy(row.get('attempt_history', []))
    history.append(dict(recovery_action=action, detected_at_unix=time.time(),
                        interruption='incomplete_checkpoint_from_previous_exclusive_run',
                        previous_attempt=previous, uncheckpointed_work_cost='unknown'))
    return history


def complete_evaluation(row, held_sensor, held_struct, evaluator=evaluate_frozen):
    """Only reads the frozen endpoint; has no data or optimizer dependency."""
    frozen = np.array(row['estimated'], dtype=float)
    audit = evaluator(frozen)
    row.update(sensor_low_band=metrics(audit['prediction'], held_sensor),
        structural_low_band=metrics(audit['structural'], held_struct),
        evaluation_seconds=audit['evaluation_seconds'], evaluation_work=audit['evaluation_work'],
        per_low_frequency={str(k): dict(sensor=metrics(audit['prediction'][i], held_sensor[i]),
            structural=metrics(audit['structural'][i], held_struct[i])) for i, k in enumerate(LOW)},
        evaluation_status='ok')


def checkpoint_checks():
    """Bounded crash-state checks; no optimizer or N32 physics is invoked."""
    assert checkpoint_action(None) == 'new'
    interrupted = dict(status='started', seed=8101, setup_seconds=.4)
    assert checkpoint_action(interrupted) == 'restart_interrupted_fit'
    history = preserve_interruption(interrupted, 'restart_interrupted_fit')
    assert history[0]['previous_attempt'] == interrupted
    snapshot = copy.deepcopy(interrupted)
    history[0]['previous_attempt']['setup_seconds'] = 999
    assert interrupted == snapshot
    frozen = dict(status='converged', estimated=BASE.tolist(), nfev=7,
                  solve_seconds=1.2, setup_seconds=.4, evaluation_status='started')
    assert checkpoint_action(frozen) == 'resume_evaluation'
    frozen['attempt_history'] = preserve_interruption(frozen, 'resume_evaluation')
    calls = []
    field = np.ones((3, 17, 3, 4), complex)
    def fixture_evaluator(z):
        calls.append(z.copy())
        return dict(prediction=field, structural=field, evaluation_seconds=.1, evaluation_work={})
    complete_evaluation(frozen, field, field, evaluator=fixture_evaluator)
    assert len(calls) == 1 and np.array_equal(calls[0], BASE)
    assert frozen['estimated'] == BASE.tolist() and frozen['nfev'] == 7
    assert frozen['solve_seconds'] == 1.2 and frozen['setup_seconds'] == .4
    assert checkpoint_action(frozen) == 'skip_completed'
    for failure in (dict(status='failed', evaluation_status='not_reached'),
                    dict(status='converged', estimated=BASE.tolist(), evaluation_status='failed')):
        assert checkpoint_action(failure) == 'skip_completed'
    return dict(passed=True, interrupted_metadata_preserved=True,
                frozen_evaluation_resumed_without_refit=True, completed_failures_not_retried=True)


def benchmark():
    # Recheck current implementation first; never use stale pass artifacts.
    checks()
    sources, cache_seconds = cached_sources()
    write_json(HERE/'reference_cache_audit.json', dict(
        cache_read_seconds=cache_seconds,
        per_frequency={str(k): v[2] for k, v in sources.items()},
        historical_generation_cost_is_not_current_read_cost=True))
    dest = HERE/'fits.json'
    rows = json.loads(dest.read_text()) if dest.exists() else []
    code_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if any(r['implementation_sha256'] != code_hash for r in rows):
        raise RuntimeError('Existing fits use different code; preserve them and choose a fresh output folder.')
    for seed in (8101, 8102):
        start = time.perf_counter()
        true, _ = scene(seed)
        mean, _ = reference_means(sources, true, receivers(), LOW+(18.,))
        data = make_data(seed, mean[:3], mean[3])
        sharing = assert_shared(data)
        held_sensor, held_struct = reference_means(sources, true, receivers(17, 1.6), LOW)
        prep_seconds = time.perf_counter()-start
        np.savez_compressed(HERE/f'data_{seed}.npz', true=true,
            sigma=data['sigma'], low=data['observations']['low'],
            low_high=data['observations']['low_high'], low_repeat=data['observations']['low_repeat'],
            reference=data['reference'], low_noise=data['low_noise'], extra_noise=data['extra_noise'],
            held_sensor_mean=held_sensor, held_structural_mean=held_struct,
            training_receivers=receivers(), held_receivers=receivers(17, 1.6))
        write_json(HERE/f'data_{seed}.json', dict(seed=seed, **sharing,
                   preparation_seconds=prep_seconds, frequencies=CHOICES,
                   reference_frequencies=LOW, reference_sigma=REF_SIGMA))
        for choice, frequencies in CHOICES.items():
            for use_reference in (False, True):
                matching = [r for r in rows if r['seed']==seed and r['choice']==choice and r['use_reference']==use_reference]
                if len(matching) > 1:
                    raise RuntimeError('Duplicate fit checkpoints; preserve and inspect them.')
                existing = matching[0] if matching else None
                action = checkpoint_action(existing)
                if action == 'skip_completed':
                    continue
                if action == 'resume_evaluation':
                    # Never re-enter least_squares for an already frozen endpoint.
                    existing['attempt_history'] = preserve_interruption(existing, action)
                    existing['evaluation_status'] = 'started'
                    write_json(dest, rows)
                    resume_start = time.perf_counter()
                    try:
                        complete_evaluation(existing, held_sensor, held_struct)
                    except Exception as exc:
                        existing.update(evaluation_status='failed', error_type=type(exc).__name__, error=str(exc))
                    existing['resumed_evaluation_wall_seconds'] = time.perf_counter()-resume_start
                    existing['uncheckpointed_interrupted_work_cost'] = 'unknown; history retained, not counted as zero'
                    write_json(dest, rows)
                    print(json.dumps(dict(seed=seed, choice=choice, use_reference=use_reference,
                          recovery=action, evaluation_status=existing['evaluation_status'])), flush=True)
                    gc.collect()
                    continue
                start = time.perf_counter()
                row = dict(seed=seed, choice=choice, use_reference=use_reference,
                           frequencies=frequencies, n=32, true=true.tolist(),
                           start=BASE.tolist(), max_nfev=35, implementation_sha256=code_hash,
                           sharing=sharing, status='started',
                           scope='Two development scenes, N64 ADDA LDR versus N32 CM+RR; no population inference')
                if action == 'restart_interrupted_fit':
                    row['attempt_history'] = preserve_interruption(existing, action)
                    rows[rows.index(existing)] = row
                else:
                    rows.append(row)
                write_json(dest, rows)
                model = None
                try:
                    model = ExplicitFrequencyModel(32, frequencies)
                    row['setup_seconds'] = time.perf_counter()-start
                    cache = {}
                    def ev(v):
                        if cache.get('key') != v.tobytes():
                            r, j = residual_jac(v, model, data, choice, use_reference)
                            cache.update(key=v.tobytes(), r=r, j=j)
                        return cache['r'], cache['j']
                    solve_start = time.perf_counter()
                    opt = least_squares(lambda v: ev(v)[0], BASE.copy(), jac=lambda v: ev(v)[1],
                        bounds=(LO, HI), x_scale=SCALE, max_nfev=35,
                        ftol=1e-9, xtol=1e-9, gtol=1e-7)
                    z = opt.x.copy()
                    row.update(status='converged' if opt.status>0 else 'iteration_limit',
                        solve_seconds=time.perf_counter()-solve_start,
                        optimizer_status=int(opt.status), optimizer_message=str(opt.message),
                        nfev=int(opt.nfev), njev=int(opt.njev), estimated=z.tolist(),
                        objective=float(opt.cost), optimality=float(opt.optimality),
                        active_mask=opt.active_mask.tolist(), training_work=model.ledger(),
                        reduced_chisquare=float(opt.fun@opt.fun/(len(opt.fun)-13)),
                        material_relative_error=float(abs(z[0]-true[0])/true[0]),
                        pose_error_m=float(np.linalg.norm(z[1:4]-true[1:4])),
                        delay_error_m=float(abs(z[4]-true[4])),
                        log_amplitude_rmse=float(np.sqrt(np.mean((z[5:9]-true[5:9])**2))),
                        gain_phase_rmse_rad=float(np.sqrt(np.mean(np.angle(np.exp(1j*(z[9:13]-true[9:13])))**2))),
                        low_band_electronic_response=metrics(electronics(z, LOW), electronics(true, LOW)),
                        evaluation_status='pending')
                    # Persist the estimate before any independent evaluation.
                    write_json(dest, rows)
                    del model
                    model = None
                    gc.collect()
                    complete_evaluation(row, held_sensor, held_struct)
                except Exception as exc:
                    row.update(status='failed' if 'estimated' not in row else row['status'],
                               error_type=type(exc).__name__, error=str(exc),
                               evaluation_status='failed' if 'estimated' in row else 'not_reached')
                    if model is not None:
                        row['partial_training_work'] = model.ledger()
                row['fit_and_evaluation_seconds'] = time.perf_counter()-start
                write_json(dest, rows)
                print(json.dumps({k: row.get(k) for k in ('seed','choice','use_reference','status','nfev','sensor_low_band','error')}), flush=True)
                del model
                gc.collect()
    comparisons = []
    for seed in (8101, 8102):
        for use_ref in (False, True):
            selected = {r['choice']: r for r in rows if r['seed']==seed and r['use_reference']==use_ref}
            for baseline in ('low', 'low_repeat'):
                high, base = selected['low_high'], selected[baseline]
                row = dict(seed=seed, use_reference=use_ref, contrast='low_high minus '+baseline,
                           negative_error_difference_means_high_better=True)
                if all(r.get('evaluation_status')=='ok' for r in (high, base)):
                    row['differences'] = {field: {metric: high[field][metric]-base[field][metric]
                        for metric in ('relative_field_error', 'weighted_phase_rmse_rad')}
                        for field in ('sensor_low_band', 'structural_low_band')}
                else:
                    row['status'] = 'not_comparable_due_to_failure'
                comparisons.append(row)
    write_json(HERE/'comparisons.json', comparisons)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', action='store_true', help='Run registered 12 fits (after parent authorization).')
    parser.add_argument('--check-reference-data', action='store_true', help='Read existing N64 caches and check actual shared data only.')
    args = parser.parse_args()
    if args.run:
        # The OS releases this lock on process termination, so an incomplete
        # checkpoint cannot belong to another active compliant runner.
        with (HERE/'benchmark.lock').open('a') as run_lock:
            fcntl.flock(run_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            benchmark()
    elif args.check_reference_data:
        sources, seconds = cached_sources()
        records = []
        for seed in (8101, 8102):
            true, _ = scene(seed)
            mean, _ = reference_means(sources, true, receivers(), LOW+(18.,))
            records.append(dict(seed=seed, **assert_shared(make_data(seed, mean[:3], mean[3]))))
        write_json(HERE/'reference_data_checks.json', dict(passed=True, cache_read_seconds=seconds, scenes=records,
                   reference_metadata={str(k): values[2] for k, values in sources.items()},
                   implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
        print(json.dumps(records), flush=True)
    else:
        checks()
