"""Validate deadline bookkeeping before any new heldout scene is opened."""
import gc
import json
from pathlib import Path
import time
import numpy as np
import timed_rom as timed

OUT=Path(__file__).resolve().parent


def main():
    sm=timed.seedmod;rows=[]
    data=sm.ctl.make_case(16,9,4101)
    for budget in [0.,10.]:
        for method in sm.METHODS:
            start=time.perf_counter();model=sm.model_for(16,9);setup=time.perf_counter()-start
            row=timed.solve(model,data,method,budget,setup,maxiter=200)
            row['seed']=4101
            if budget==0:
                assert not row['events'] and row['status']=='decision_deadline'
                assert np.allclose(row['alpha'],.08) and np.allclose(row['pose'],data['x0'])
                assert row['online_work']['full_rhs_columns']==0
            else:
                assert row['events'] and all(e['seconds']<budget for e in row['events'])
            rows.append(row)
            sm.ctl.write_json(OUT/'results/timed_rom_development.json',rows)
            print(json.dumps({k:row[k] for k in ['method','decision_budget_seconds','status','actual_consumed_seconds','final_loss','reduced_updates']}),flush=True)
            del model;gc.collect()
    old=json.loads((OUT/'results/rom_seed_q9.json').read_text())
    for row in rows:
        if row['decision_budget_seconds']==10 and row['method']!='direct_adjoint':
            previous=next(r for r in old if r['seed']==4101 and r['method']==row['method'])
            assert abs(previous['final_loss']-row['final_loss'])<1e-7
    sm.ctl.write_json(OUT/'results/timed_rom_checks.json',dict(passed=True,
        checks=['zero-deadline accepts no update','positive-deadline timestamps valid',
                'GN/ROM unrestricted endpoints match prior development controller'],
        scope='one inspected development scene; direct adjoint allowed 200 iterations rather than 24'))


if __name__=='__main__':main()
