"""Independent fresh validation after adding data and rebuilding a failed bank.

Uses saved first-pass fitted parameters, not truth, as starts. The first-pass
validation batch may become new training data; its noise is then not reused
for verification. All policies use the same cost/cardinality and fresh k9
verification. Four scenes remain development, not 800 independent scenes.
"""
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2
from branch_acquisition import OUT,A2,pack,unpack,realify,bounds,POSE_SCALE,rows_for,TRAIN_RX,lever_metric_error


def refit(model, z0, trainidx, ytrain, sigma, max_nfev=45):
    # Training never includes frequency index 2 (reserved fresh validation).
    ids=[0,1,3]
    full_rows=np.concatenate([np.arange(fi*72,(fi+1)*72) for fi in ids])
    lookup={int(v):i for i,v in enumerate(full_rows)}
    idx=np.array([lookup[int(v)] for v in trainidx])
    cache={};rhs=0
    def ev(z):
        nonlocal rhs
        if cache.get('key')!=z.tobytes():
            a,x=unpack(z);fw=model.forward(a,x,ids,jacobian=True)
            rhs+=fw['work']['rhs_solves_total']
            cache.update(key=z.tobytes(),r=realify(fw['total'][idx]-ytrain,sigma),
                         j=realify(np.c_[fw['A'],fw['B']/POSE_SCALE][idx],sigma))
        return cache['r'],cache['j']
    bb=np.array(bounds(9))
    opt=least_squares(lambda z:ev(z)[0],z0,jac=lambda z:ev(z)[1],bounds=(bb[:,0],bb[:,1]),
                      max_nfev=max_nfev,ftol=1e-9,xtol=1e-9,gtol=1e-7)
    a,x=unpack(opt.x)
    fw=model.forward(a,x,[2],jacobian=False)
    rhs+=fw['work']['rhs_solves_total']
    return dict(alpha=a,x=x,pred=fw['total'],nfev=int(opt.nfev),status=int(opt.status),
                loss=float(opt.fun@opt.fun/2),rhs=rhs)


def run_scene(first):
    seed=first['seed'];a=np.array(first['true_alpha']);x=np.array(first['true_pose']);sigma=first['sigma']
    truth=A2.Model(A2.Config(N=32,kmax=12.)).forward(a,x)['total']
    # Reproduce the original training noise, including its preceding scene draws.
    rng=np.random.default_rng(seed)
    rng.choice(9,3,False);rng.uniform(.35,.9,3)
    oldidx=rows_for(3,TRAIN_RX)
    yold=truth[oldidx]+sigma/np.sqrt(2)*(rng.normal(size=len(oldidx))+1j*rng.normal(size=len(oldidx)))
    # This is a newly drawn acquisition batch, independent of the first report.
    arng=np.random.default_rng(seed+90000)
    acquisition_noise=sigma/np.sqrt(2)*(arng.normal(size=len(truth))+1j*arng.normal(size=len(truth)))
    vrng=np.random.default_rng(seed+100000)
    valtrue=truth[144:216]
    valnoise=sigma/np.sqrt(2)*(vrng.normal(size=(200,72))+1j*vrng.normal(size=(200,72)))
    yy=np.array([realify(valtrue+v,sigma) for v in valnoise])
    policies=[]
    for policy in first['policies']:
        t0=time.perf_counter()
        newidx=rows_for(policy['fi'],policy['rx'])
        idx=np.r_[oldidx,newidx];y=np.r_[yold,truth[newidx]+acquisition_noise[newidx]]
        bank=[]
        for b in first['bank']:
            rec=refit(A2.Model(A2.Config(N=16,kmax=12.)),pack(np.array(b['alpha']),np.array(b['x'])),idx,y,sigma)
            rec['pose_error']=lever_metric_error(rec['x'],x)
            rec['material_error']=float(np.linalg.norm(rec['alpha']-a)/np.linalg.norm(a))
            bank.append(rec)
        mu=np.array([realify(b['pred'],sigma) for b in bank])
        costs=np.sum((yy[:,None,:]-mu[None,:,:])**2,axis=-1)
        choice=costs.argmin(axis=1)
        covers=np.array([b['pose_error']<.05 and b['material_error']<.25 for b in bank])
        rejected=costs[np.arange(len(yy)),choice]>chi2.ppf(.99,144)
        good=covers[choice]
        row=dict(policy=policy['policy'],fi=policy['fi'],rx=policy['rx'],
                 bank_covers_tolerance=bool(covers.any()),accurate_selection_rate=float(good.mean()),
                 rejection_rate=float(rejected.mean()),false_accept_rate=float(np.mean(~good & ~rejected)),
                 selected_pose_error=float(np.mean([bank[i]['pose_error'] for i in choice])),
                 selected_material_error=float(np.mean([bank[i]['material_error'] for i in choice])),
                 rhs_columns=sum(b['rhs'] for b in bank),wall_seconds=time.perf_counter()-t0,
                 candidate_results=[{k:v for k,v in b.items() if k not in ['alpha','x','pred']} for b in bank],
                 no_global_coverage_guarantee=True,certified_beta_available=False)
        policies.append(row)
        print(json.dumps({k:row[k] for k in ['policy','bank_covers_tolerance','accurate_selection_rate','rejection_rate','selected_pose_error','wall_seconds']},default=float),flush=True)
    return dict(seed=seed,policies=policies,
                scope='development rebuild after failed-bank validation; k9 observations kept independent of all fitting and policy selection')


if __name__=='__main__':
    initial=json.loads((OUT/'results'/'branch_acquisition_development.json').read_text())
    dest=OUT/'results'/'branch_repair_development.json'
    rows=json.loads(dest.read_text()) if dest.exists() else []
    for rec in initial:
        if any(r['seed']==rec['seed'] for r in rows):continue
        rows.append(run_scene(rec))
        dest.write_text(json.dumps(rows,indent=2)+'\n')
