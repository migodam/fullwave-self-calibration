"""Anytime decision-budget controller. Never accept an iterate after deadline.

Atomic linear algebra cannot be preempted: record overrun as real consumed
runtime. Compare accepted states available by the decision deadline, not an
untrue claim of exactly equal consumed wall time. Offline audits are separate.
"""
import time
import numpy as np
from scipy.optimize import minimize
import rom_seed_comparison as seedmod

ctl=seedmod.ctl


class Deadline(RuntimeError):pass


def solve(model,data,method,seconds,setup_seconds=0.,maxiter=200):
    ledger=seedmod.SeedLedger(model.n_cells)
    start=time.perf_counter()-setup_seconds
    events=[];z=ctl.pack(np.full(model.n_alpha,.08),data['x0'])
    accepted=z.copy();chart=None;status='iteration_limit'
    selected=ctl.select_rows(data['y'],ctl.FREQ_IDS)
    bounds=np.array(ctl.bounds(model.n_alpha))

    def elapsed():return time.perf_counter()-start
    def check():
        if elapsed()>=seconds:raise Deadline()
    def record(zz,loss,reduced=False):
        nonlocal accepted
        # Check after physical solves and immediately before publishing state.
        check();t=elapsed()
        if t>=seconds:raise Deadline()
        accepted=zz.copy()
        events.append(dict(seconds=t,z=zz.tolist(),loss=float(loss),reduced=reduced))
    def forward(zz,jac=False):
        check();a,x=ctl.unpack(zz)
        out=model.forward(a,x,list(ctl.FREQ_IDS),jacobian=jac)
        ledger.charge_model_work(out['work'])
        if not jac:ledger.charge_acceptance()
        check();return out

    fallback=0;rebuild=0;reuse=0
    try:
        if method=='direct_adjoint':
            cache={}
            def fg(zz):
                check();key=zz.tobytes()
                if cache.get('key')!=key:
                    a,x=ctl.unpack(zz)
                    out=model.adjoint_gradient(a,x,selected,data['sigma'],list(ctl.FREQ_IDS))
                    ledger.charge_model_work(out['work'])
                    cache.update(key=key,value=(out['loss'],np.r_[out['grad_alpha'],out['grad_x']/ctl.POSE_SCALE]))
                check();return cache['value']
            def callback(zz):
                loss,_=fg(zz);record(zz,loss)
            opt=minimize(fg,z,jac=True,method='L-BFGS-B',bounds=ctl.bounds(model.n_alpha),
                callback=callback,options=dict(maxiter=maxiter,maxls=24,ftol=1e-12,gtol=1e-7))
            # A converged final point may duplicate the last accepted callback.
            record(opt.x,opt.fun);status=str(opt.message)
        else:
            fw=forward(z)
            record(z,ctl.exact_loss(fw['total']-selected,data['sigma']))
            for it in range(maxiter):
                check();a,x=ctl.unpack(z)
                loss=ctl.exact_loss(fw['total']-selected,data['sigma']);red=None
                if method!='direct_gn':
                    if chart is not None:
                        red=ctl.fixed_chart_evaluate(model,a,x,chart,data['y'],data['sigma'],ledger=ledger)
                        check()
                        if ctl.chart_guard(red,fw)['passed']:reuse+=1
                        else:red=None
                    if red is None:
                        chart=seedmod.build(model,a,x,ctl.FREQ_IDS,method,128,data['y'],ledger)
                        rebuild+=1;check()
                        red=ctl.fixed_chart_evaluate(model,a,x,chart,data['y'],data['sigma'],ledger=ledger)
                        check()
                        if not ctl.chart_guard(red,fw)['passed']:red=None;fallback+=1
                if red is None:
                    fw=forward(z,True)
                    jac=ctl.scaled_reduced_jacobian(fw['A'],fw['B'],data['sigma'])
                    residual=ctl.realify(fw['total']-selected,data['sigma'])
                else:
                    jac=ctl.scaled_reduced_jacobian(red['A'],red['B'],data['sigma'])
                    residual=ctl.realify(red['total']-selected,data['sigma'])
                moved=False
                for damp in ctl.DAMPING_SCHEDULE:
                    check();step,_=ctl._gn_trial_step(jac,residual,damp,ledger)
                    trial=np.clip(z+step,bounds[:,0],bounds[:,1]);check()
                    if np.linalg.norm(trial-z)<1e-12:continue
                    ft=forward(trial);lt=ctl.exact_loss(ft['total']-selected,data['sigma'])
                    if ctl.accept_trial(loss,lt,np.linalg.norm(trial-z)):
                        record(trial,lt,red is not None)
                        z,fw=trial,ft;moved=True;break
                if not moved:status='no_improving_trial';break
    except Deadline:
        status='decision_deadline'
    online=elapsed()
    assert all(e['seconds']<seconds for e in events)
    work=ledger.snapshot()
    # Frozen decision state only; truth is accessed after all online decisions.
    audit_start=time.perf_counter();a,x=ctl.unpack(accepted)
    f=model.forward(a,x,list(ctl.FREQ_IDS),jacobian=False)
    return dict(method=method,status=status,decision_budget_seconds=seconds,
        actual_consumed_seconds=online,atomic_overrun_seconds=max(0.,online-seconds),
        setup_seconds=setup_seconds,events=events,alpha=a.tolist(),pose=x.tolist(),
        final_loss=ctl.exact_loss(f['total']-selected,data['sigma']),
        pose_error_m=ctl.lever_metric_error(x,data['x_true']),
        material_rmse=float(np.linalg.norm(a-data['alpha_true'])/np.sqrt(model.n_alpha)),
        reduced_updates=sum(e['reduced'] for e in events),rebuilds=rebuild,reuse_hits=reuse,
        fallback_steps=fallback,online_work=work,offline_work=f['work'],
        offline_seconds=time.perf_counter()-audit_start,
        scope='anytime accepted states; atomic overrun disclosed, not exact equal consumed wall time')
