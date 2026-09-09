"""Matched development seed ablation; frozen prior modules imported read-only."""
import argparse
import gc
import json
from pathlib import Path
import sys
import time
import numpy as np

OUT = Path(__file__).resolve().parent
REUSE = OUT.parents[1] / 'delegated/a3_rom_reuse'
sys.path.insert(0, str(REUSE))
import controller as ctl
import cost
from basis import (sensing_chart, block_krylov_chart, _reduced_states_on_U,
                   _append_orth)

ORIGINAL = cost.build_chart_charged
METHODS = ['direct_gn', 'direct_adjoint', 'generic_task', 'sensing_task',
           'twofold_task', 'krylov_task']


class SeedLedger(cost.AuditedLedger):
    def snapshot(self):
        result = super().snapshot()
        result['seed_events'] = getattr(self, 'seed_events', [])
        return result


def build(model, alpha, x, freq_ids, method, rank, y_full, ledger,
          rank_by_freq=None, include_tangent_rhs=True, **kwargs):
    if method not in ('sensing_task', 'krylov_task'):
        return ORIGINAL(model, alpha, x, freq_ids, method, rank, y_full, ledger,
                        rank_by_freq=rank_by_freq,
                        include_tangent_rhs=include_tangent_rhs, **kwargs)
    geom = model._effective_geometry(x)
    charts = {}
    for fi in sorted(set(freq_ids)):
        target = int((rank_by_freq or {}).get(fi, rank))
        initial = min(64, target)
        if method == 'sensing_task':
            u = sensing_chart(model, geom, fi, initial, ledger)
        else:
            u = block_krylov_chart(model, alpha, geom, fi, initial, ledger)
        seed_rank = u.shape[1]
        yr = list(np.asarray(y_full)[fi*72:(fi+1)*72].reshape(-1, model.n_rx))
        for _ in range(14):
            if u.shape[1] >= target:
                break
            records = _reduced_states_on_U(model, alpha, geom, fi, u, ledger)
            directions = cost._task_directions_charged(
                model, records, yr, include_tangent_rhs, ledger, fi)
            if not directions:
                break
            nxt = _append_orth(u, np.column_stack(directions), ledger, max_rank=target)
            if nxt.shape[1] <= u.shape[1]:
                break
            u = nxt
        if not hasattr(ledger, 'seed_events'):
            ledger.seed_events = []
        ledger.seed_events.append(dict(method=method, fi=fi, seed_rank=seed_rank,
                                       final_rank=u.shape[1]))
        charts[fi] = u.copy()
    return charts


def model_for(n, q):
    if q == 9:
        return ctl.A2.Model(ctl.A2.Config(N=n, kmax=12.))
    pts, _ = ctl.A2.make_grid(n)
    axis = np.linspace(-.65, .65, int(np.sqrt(q)))
    centers = np.array([(i,j) for i in axis for j in axis])
    mb = np.exp(-np.sum((pts[:,None,:]-centers[None,:,:])**2,axis=2)/(2*.16**2))
    return ctl.A2.Model(ctl.A2.Config(N=n, kmax=12., material_basis=mb))


def check():
    model = model_for(16,9)
    data = ctl.make_case(16,9,4101)
    a, x = np.full(9,.08), data['x0']
    checks = []
    for method in METHODS[2:]:
        ledger = cost.AuditedLedger(model.n_cells)
        chart = build(model,a,x,[3],method,128,data['y'],ledger)
        u = chart[3]
        orth = float(np.linalg.norm(u.conj().T@u-np.eye(u.shape[1])))
        assert orth < 1e-10
        if method in ('generic_task','twofold_task'):
            old = ORIGINAL(model,a,x,[3],method,128,data['y'],cost.AuditedLedger(model.n_cells))
            assert np.array_equal(u,old[3])
        red = ctl.fixed_chart_evaluate(model,a,x,chart,data['y'],data['sigma'],freq_ids=(3,),ledger=ledger)
        # Differentiate this fixed chart, not a newly rebuilt chart.
        h = 1e-5
        plus = ctl.fixed_chart_evaluate(model,a,x+np.array([h,0,0]),chart,data['y'],data['sigma'],freq_ids=(3,))
        minus = ctl.fixed_chart_evaluate(model,a,x-np.array([h,0,0]),chart,data['y'],data['sigma'],freq_ids=(3,))
        fd = (plus['total']-minus['total'])/(2*h)
        err = float(np.linalg.norm(fd-red['B'][:,0])/np.linalg.norm(fd))
        assert err < 1e-6, (method,err)
        checks.append(dict(method=method,rank=u.shape[1],orthogonality_error=orth,
                           fixed_chart_pose_derivative_error=err,
                           seed_events=getattr(ledger,'seed_events',[])))
    ctl.write_json(OUT/'results/rom_seed_checks.json',dict(passed=True,checks=checks))
    print(json.dumps(checks),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--test',action='store_true')
    p.add_argument('--q',type=int,choices=[9,49],default=9)
    args = p.parse_args()
    if args.test:
        check()
    else:
        # Replace only the dispatch binding in this process, never frozen files.
        ctl.build_chart_charged = build
        ctl.AuditedLedger = SeedLedger
        n = 16 if args.q == 9 else 32
        seeds = [4101,4102,4103,4104] if args.q == 9 else [4201,4202]
        dest = OUT/f'results/rom_seed_q{args.q}.json'
        rows = json.loads(dest.read_text()) if dest.exists() else []
        for seed in seeds:
            data = ctl.make_case(n,args.q,seed)
            for method in METHODS:
                if any(r['seed']==seed and r['method']==method for r in rows):
                    continue
                start = time.perf_counter()
                model = model_for(n,args.q)
                setup = time.perf_counter()-start
                row = ctl.solve(model,data,method,24,8000,{f:128 for f in ctl.FREQ_IDS})
                row.update(seed=seed,model_setup_seconds=setup,
                    cold_seconds=setup+row['online_seconds'],
                    scope='serial DEVELOPMENT seed comparison; not equal wall budgets or final testing')
                rows.append(row)
                ctl.write_json(dest,rows)
                print(json.dumps({k:row[k] for k in ['seed','method','cold_seconds',
                    'pose_error','material_rmse','final_loss','reduced_steps','fallback_steps']}),flush=True)
                del model
                gc.collect()
