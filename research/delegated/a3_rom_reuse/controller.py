"""Parent takeover: budget-safe fixed-chart reuse, development only.

Current coefficients remain physically constrained: chart rank is a numerical
approximation dimension, not an independent nuisance-current dimension.
Offline exact derivatives NEVER enter cache admission or update decisions.
"""
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parent / 'a3_rom'))
from common import A2, FREQ_IDS, POSE_SCALE, realify, exact_loss, select_rows, write_json
from common import make_scene, scene_arrays, pose_init, lever_metric_error
from wave import fixed_chart_evaluate, scaled_reduced_jacobian
from solver import _gn_trial_step, DAMPING_SCHEDULE, accept_trial
sys.path.insert(0, str(OUT))
from cost import AuditedLedger, build_chart_charged


class BudgetStop(RuntimeError):
    pass


def bounds(q):
    return [(0., 2.5)] * q + [(-.7*s, .7*s) for s in POSE_SCALE]


def pack(a, x):
    return np.r_[a, x * POSE_SCALE]


def unpack(z):
    return z[:-3], z[-3:] / POSE_SCALE


def required_rhs(model, kind, ids=FREQ_IDS):
    states = len(ids)*model.n_pose*model.n_tx
    return states * ({'forward': 1, 'jacobian': 1+model.n_alpha+3, 'adjoint': 2}[kind])


def chart_guard(red, fw, tolerance=1e-3):
    # Normalize the output discrepancy by scattered, not dominant direct,
    # fields; this avoids hiding poor scattering prediction behind direct data.
    output = float(np.linalg.norm(red['total']-fw['total']) /
                   max(np.linalg.norm(fw['scattered']), 1e-14))
    state = float(red['max_state_res_rel'])
    rank_ok = all(v['rank_ok'] for v in red['per_freq'].values())
    return dict(passed=bool(rank_ok and np.isfinite(state) and state <= tolerance and output <= tolerance),
                state_residual=state, scattered_output_error=output, rank_ok=rank_ok)


def solve(model, data, method, maxiter=24, budget=8000, ranks=None, reuse=True):
    q = model.n_alpha
    led = AuditedLedger(model.n_cells)
    z = pack(np.full(q, .08), data['x0'])
    box = np.array(bounds(q))
    events, snapshots = [], []
    chart = None
    status = 'iteration_limit'
    start = time.perf_counter()
    selected = select_rows(data['y'], FREQ_IDS)

    def forward(zv, jac=False):
        cost = required_rhs(model, 'jacobian' if jac else 'forward')
        if led.full_rhs_columns + cost > budget:
            raise BudgetStop()
        a, x = unpack(zv)
        fw = model.forward(a, x, list(FREQ_IDS), jacobian=jac)
        assert fw['work']['rhs_solves_total'] == cost
        led.charge_model_work(fw['work'])
        if not jac:
            led.charge_acceptance()
        return fw

    cache = {}
    accepted = z.copy()
    if method == 'direct_adjoint':
        def fg(zz):
            key = zz.tobytes()
            if key in cache:
                return cache[key]
            cost = required_rhs(model, 'adjoint')
            if led.full_rhs_columns + cost > budget:
                raise BudgetStop()
            a, x = unpack(zz)
            out = model.adjoint_gradient(a, x, selected, data['sigma'], list(FREQ_IDS))
            assert out['work']['rhs_solves_total'] == cost
            led.charge_model_work(out['work'])
            cache[key] = (out['loss'], np.r_[out['grad_alpha'], out['grad_x']/POSE_SCALE])
            return cache[key]
        def cb(zz):
            nonlocal accepted
            accepted = zz.copy()
        try:
            res = minimize(fg, z, jac=True, method='L-BFGS-B', bounds=bounds(q), callback=cb,
                           options=dict(maxiter=maxiter, maxls=24, ftol=1e-12, gtol=1e-7))
            z = res.x
            status = str(res.message)
        except BudgetStop:
            z = accepted
            status = 'budget_stopped_before_call'
    else:
        try:
            fw = forward(z)
            for it in range(maxiter):
                a, x = unpack(z)
                L = exact_loss(fw['total']-selected, data['sigma'])
                ev = dict(iteration=it, loss=L, reused=False, rebuilt=False, fallback=False)
                red = None
                if method != 'direct_gn':
                    if chart is not None and reuse:
                        red = fixed_chart_evaluate(model, a, x, chart, data['y'], data['sigma'], ledger=led)
                        ev['reuse_guard'] = chart_guard(red, fw)
                        ev['reused'] = ev['reuse_guard']['passed']
                    if not ev['reused']:
                        chart = build_chart_charged(model, a, x, FREQ_IDS, method, 128,
                                                   data['y'], led, rank_by_freq=ranks)
                        ev['rebuilt'] = True
                        red = fixed_chart_evaluate(model, a, x, chart, data['y'], data['sigma'], ledger=led)
                    ev['guard'] = chart_guard(red, fw)
                    ev['ranks'] = {str(k): int(v.shape[1]) for k,v in chart.items()}
                    if not ev['guard']['passed']:
                        ev['fallback'] = True
                        red = None
                if red is None:
                    fw = forward(z, jac=True)
                    J = scaled_reduced_jacobian(fw['A'], fw['B'], data['sigma'])
                    r = realify(fw['total']-selected, data['sigma'])
                else:
                    J = scaled_reduced_jacobian(red['A'], red['B'], data['sigma'])
                    r = realify(red['total']-selected, data['sigma'])
                accepted_step = False
                # Reduced and direct updates share damping, clipping and exact acceptance.
                for damp in DAMPING_SCHEDULE:
                    step, _ = _gn_trial_step(J, r, damp, led)
                    trial = np.clip(z+step, box[:,0], box[:,1])
                    if np.linalg.norm(trial-z) < 1e-12:
                        continue
                    ft = forward(trial)
                    Lt = exact_loss(ft['total']-selected, data['sigma'])
                    if accept_trial(L, Lt, np.linalg.norm(trial-z)):
                        if red is not None:
                            # Immutable cache reference; all audits deferred until timer stops.
                            snapshots.append((z.copy(), chart, ev))
                        ev.update(accepted_reduced=red is not None, trial_loss=Lt, damping=damp)
                        z, fw = trial, ft
                        accepted_step = True
                        break
                events.append(ev)
                if not accepted_step:
                    status = 'no_improving_trial'
                    break
        except BudgetStop:
            status = 'budget_stopped_before_call'
    online_seconds = time.perf_counter()-start
    online = led.snapshot()
    assert online['full_rhs_columns'] <= budget
    # Separate, fully charged audits; no audit selection affects optimizer.
    audit = AuditedLedger(model.n_cells)
    audit_start = time.perf_counter()
    for zz, cc, ev in snapshots:
        aa, xx = unpack(zz)
        rr = fixed_chart_evaluate(model, aa, xx, cc, data['y'], data['sigma'], freq_ids=(3,), ledger=audit)
        ff = model.forward(aa, xx, [3], jacobian=True)
        audit.charge_model_work(ff['work'])
        ev['offline_high_frequency'] = {
            'map_jac_relative_error': float(np.linalg.norm(rr['A']-ff['A'])/max(np.linalg.norm(ff['A']),1e-14)),
            'pose_jac_relative_error': float(np.linalg.norm(rr['B']-ff['B'])/max(np.linalg.norm(ff['B']),1e-14)),
            'state_residual': rr['max_state_res_rel']}
    a, x = unpack(z)
    ff = model.forward(a,x,list(FREQ_IDS),jacobian=True)
    audit.charge_model_work(ff['work'])
    phase = np.angle(ff['scattered'] * data['true_scattered'].conj())
    weights = np.abs(data['true_scattered'])**2
    return dict(method=method, q=q, grid=model.config.N, status=status, alpha=a, pose=x,
                pose_error=lever_metric_error(x, data['x_true']), material_rmse=float(np.linalg.norm(a-data['alpha_true'])/np.sqrt(q)),
                final_loss=exact_loss(ff['total']-selected,data['sigma']),
                phase_rms=float(np.sqrt(np.sum(weights*phase**2)/weights.sum())),
                online_seconds=online_seconds, audit_seconds=time.perf_counter()-audit_start,
                online_work=online, offline_work=audit.snapshot(), events=events,
                reduced_steps=sum(e.get('accepted_reduced',False) for e in events),
                reuse_hits=sum(e.get('reused',False) for e in events),
                rebuilds=sum(e.get('rebuilt',False) for e in events),
                fallback_steps=sum(e.get('fallback',False) for e in events),
                budget_rhs=budget, ranks_requested=ranks, reuse_enabled=reuse,
                scope='development; timings not formal; same-discretization 2D solver comparison')


def make_case(n, q, seed):
    if q == 9:
        model = A2.Model(A2.Config(N=n, kmax=12.))
        d = scene_arrays(make_scene(seed))
        d['x0'] = pose_init(seed)
    else:
        pts, h = A2.make_grid(n)
        axis = np.linspace(-.65, .65, int(np.sqrt(q)))
        centers = np.array([(i,j) for i in axis for j in axis])
        mb = np.exp(-np.sum((pts[:,None,:]-centers[None,:,:])**2,axis=2)/(2*.16**2))
        model = A2.Model(A2.Config(N=n, kmax=12., material_basis=mb))
        rng = np.random.default_rng(seed)
        a = np.full(q, .08)
        a[rng.choice(q,6,replace=False)] += rng.uniform(.25,.7,6)
        clean = model.forward(a,np.zeros(3),jacobian=False)['total']
        sigma = float(np.sqrt(np.mean(np.abs(clean)**2)/1000))
        d = dict(alpha_true=a, x_true=np.zeros(3), sigma=sigma,
                 y=clean+sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape)),
                 x0=np.array([.05,-.03,.025]))
    d['true_scattered'] = model.forward(d['alpha_true'],d['x_true'],list(FREQ_IDS),jacobian=False)['scattered']
    return d


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--grid',type=int,default=16)
    p.add_argument('--q',type=int,default=9)
    p.add_argument('--seeds',type=int,nargs='+',default=[4101,4102,4103,4104])
    p.add_argument('--methods',nargs='+',default=['direct_gn','direct_adjoint','generic_task','twofold_task'])
    p.add_argument('--iterations',type=int,default=24)
    p.add_argument('--budget',type=int,default=8000)
    p.add_argument('--rank',type=int,default=128)
    args = p.parse_args()
    dest = OUT / f'reuse_n{args.grid}_q{args.q}_r{args.rank}.json'
    rows = json.loads(dest.read_text()) if dest.exists() else []
    for seed in args.seeds:
        d = make_case(args.grid,args.q,seed)
        for method in args.methods:
            if any(r['seed']==seed and r['method']==method for r in rows):
                continue
            start = time.perf_counter()
            if args.q == 9:
                model = A2.Model(A2.Config(N=args.grid, kmax=12.))
            else:
                pts,_ = A2.make_grid(args.grid)
                axis=np.linspace(-.65,.65,int(np.sqrt(args.q)))
                centers=np.array([(i,j) for i in axis for j in axis])
                mb=np.exp(-np.sum((pts[:,None,:]-centers[None,:,:])**2,axis=2)/(2*.16**2))
                model=A2.Model(A2.Config(N=args.grid,kmax=12.,material_basis=mb))
            setup = time.perf_counter()-start
            row=solve(model,d,method,args.iterations,args.budget,{f:args.rank for f in FREQ_IDS})
            row.update(seed=seed,model_setup_seconds=setup,cold_seconds=setup+row['online_seconds'])
            rows.append(row)
            write_json(dest,rows)
            print(json.dumps({k:row[k] for k in ['seed','method','q','final_loss','pose_error','online_seconds','reduced_steps','reuse_hits','fallback_steps']}),flush=True)


if __name__ == '__main__':
    main()
