"""Physical 2D full-wave development bank/heldout-acquisition experiment.

Grid32 generates truth, grid16 fits nine material coefficients and three rigid
pose coordinates. Candidate construction uses training rows only. Selection
uses fresh noise, with a fixed candidate bank; repetition is CONDITIONAL noise
Monte Carlo, not independent scenes. No certified discrepancy bound is known.
"""
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2

OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT.parents[1]/'delegated'/'a3_rom_reuse'))
from controller import A2, POSE_SCALE, bounds, pack, unpack, realify, lever_metric_error
from gain_graph import branch_bound, span

TRAIN_RX=[0,6]

def rows_for(fi, rxs):
    return np.array([fi*72+p*24+t*12+r for p in range(3) for t in range(2) for r in rxs])


def candidate_fit(model, y, sigma, x0, max_nfev):
    indices=rows_for(0,TRAIN_RX) # relative to single-frequency output
    box=np.array(bounds(9)); cache={}
    rhs=0; t0=time.perf_counter()
    def ev(z):
        nonlocal rhs
        key=z.tobytes()
        if cache.get('key')!=key:
            a,x=unpack(z)
            fw=model.forward(a,x,[3],jacobian=True)
            rhs+=fw['work']['rhs_solves_total']
            jc=np.c_[fw['A'],fw['B']/POSE_SCALE]
            cache.update(key=key,r=realify(fw['total'][indices]-y,sigma),
                         j=realify(jc[indices],sigma))
        return cache['r'],cache['j']
    fit=least_squares(lambda z:ev(z)[0],pack(np.full(9,.08),x0),jac=lambda z:ev(z)[1],
                      bounds=(box[:,0],box[:,1]),max_nfev=max_nfev,
                      ftol=1e-9,xtol=1e-9,gtol=1e-7)
    a,x=unpack(fit.x)
    prediction=model.forward(a,x,[0,1,2,3],jacobian=True)
    rhs+=prediction['work']['rhs_solves_total']
    return dict(alpha=a,x=x,training_loss=float(fit.fun@fit.fun/2),
                status=int(fit.status),nfev=int(fit.nfev),rhs_columns=rhs,
                wall_seconds=time.perf_counter()-t0,
                prediction=prediction['total'],A=prediction['A'],B=prediction['B'])


def scene(seed, candidates=8, repetitions=200, max_nfev=50):
    rng=np.random.default_rng(seed)
    a=np.full(9,.08)
    a[rng.choice(9,3,False)]+=rng.uniform(.35,.9,3)
    x=np.array([.025,-.018,.012])
    true_model=A2.Model(A2.Config(N=32,kmax=12.))
    truth=true_model.forward(a,x,jacobian=False)['total']
    sigma=float(np.sqrt(np.mean(abs(truth)**2)/1000))
    trainidx=rows_for(3,TRAIN_RX)
    ytrain=truth[trainidx]+sigma/np.sqrt(2)*(rng.normal(size=len(trainidx))+1j*rng.normal(size=len(trainidx)))
    # Deterministic multistart does NOT use truth to position initial guesses.
    starts=[np.zeros(3)]
    for j in range(candidates-1):
        theta=2*np.pi*j/(candidates-1)
        starts.append(np.array([.45*np.cos(theta),.45*np.sin(theta),.2*(-1)**j]))
    bank=[]
    for j,start in enumerate(starts):
        rec=candidate_fit(A2.Model(A2.Config(N=16,kmax=12.)),ytrain,sigma,start,max_nfev)
        rec['pose_error']=lever_metric_error(rec['x'],x)
        rec['material_relative_error']=float(np.linalg.norm(rec['alpha']-a)/np.linalg.norm(a))
        bank.append(rec)
        print(json.dumps(dict(seed=seed,candidate=j,training_loss=rec['training_loss'],pose_error=rec['pose_error'],nfev=rec['nfev'])),flush=True)
    # No oracle or validation-based bank pruning. Deduplicate only parameter-
    # close solutions to avoid a zero separation score for repeated minima.
    training_seconds_total=sum(b['wall_seconds'] for b in bank)
    training_rhs_total=sum(b['rhs_columns'] for b in bank)
    unique=[]
    for rec in sorted(bank,key=lambda r:r['training_loss']):
        if not any(np.linalg.norm(pack(rec['alpha'],rec['x'])-pack(s['alpha'],s['x']))<1e-3 for s in unique):
            unique.append(rec)
    bank=unique
    pool=[]
    for fi in [0,1,3]:
        for rx in [[1,2,7,8],[3,4,9,10],[2,5,8,11]]:
            pool.append(dict(fi=fi,rx=rx,indices=rows_for(fi,rx)))
    design_start=time.perf_counter()
    # Both scores depend only on the training-frozen bank/model.
    best=bank[0]
    for acq in pool:
        idx=acq['indices']
        mu=np.array([realify(b['prediction'][idx],sigma) for b in bank])
        cert=branch_bound(mu,0.) if len(bank)>1 else dict(minimum_distance=0.,worst_conditional_bound=1.)
        aa=realify(best['A'][idx],sigma)
        bb=realify(best['B'][idx]/POSE_SCALE,sigma)
        nuisance=span(aa)
        sv=np.linalg.svd(bb-nuisance@(nuisance.T@bb),compute_uv=False)
        acq.update(mu=mu,separation_score=cert['minimum_distance'],
                   fisher_score=float(sv[-1]**2),nominal_bank_risk=cert['worst_conditional_bound'])
    design_seconds=time.perf_counter()-design_start
    chosen={'fixed':0,'random':int(rng.integers(len(pool))),
            'local_information':int(np.argmax([v['fisher_score'] for v in pool])),
            'branch_separation':int(np.argmax([v['separation_score'] for v in pool]))}
    covered=np.array([b['pose_error']<.05 and b['material_relative_error']<.25 for b in bank])
    policies=[]
    # Same underlying noise array across policies (paired; overlap correlated).
    valrng=np.random.default_rng(seed+80000)
    noise=(valrng.normal(size=(repetitions,len(truth)))+1j*valrng.normal(size=(repetitions,len(truth))))*sigma/np.sqrt(2)
    for policy,which in chosen.items():
        acq=pool[which];idx=acq['indices'];mu=acq['mu']
        vv=np.array([realify(truth[idx]+e[idx],sigma) for e in noise])
        costs=np.sum((vv[:,None,:]-mu[None,:,:])**2,axis=-1)
        selected=costs.argmin(axis=1)
        # Noise-only residual envelope; mismatch bound beta is NOT available.
        # Union adjustment controls rejection of a correct *exact* predictor,
        # not false acceptance of an arbitrary wrong physical solution.
        threshold=float(chi2.ppf(.99,len(vv[0])))
        rejected=costs[np.arange(repetitions),selected]>threshold
        accurate=covered[selected]
        mue=realify(truth[idx],sigma)
        mismatch=np.linalg.norm(mu-mue[None,:],axis=1)
        policies.append(dict(policy=policy,pool_index=which,fi=acq['fi'],rx=acq['rx'],
             nominal_separation=acq['separation_score'],nominal_bank_risk_beta0=acq['nominal_bank_risk'],
             selection_accurate_rate=float(accurate.mean()),rejection_rate=float(rejected.mean()),
             false_accept_rate=float(np.mean(~accurate & ~rejected)),
             selected_pose_error_mean=float(np.mean([bank[j]['pose_error'] for j in selected])),
             selected_material_error_mean=float(np.mean([bank[j]['material_relative_error'] for j in selected])),
             true_prediction_discrepancy_whitened=mismatch.tolist(),
             selected_counts=np.bincount(selected,minlength=len(bank)).tolist(),
             coverage_guarantee=False,certified_beta_available=False,
             conditional_bound_applicable=False,validation_refitting=False))
    compact=[{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in b.items() if k not in ['prediction','A','B']} for b in bank]
    return dict(seed=seed,true_alpha=a.tolist(),true_pose=x.tolist(),sigma=sigma,
                bank=compact,bank_covers_tolerance=bool(covered.any()),coverage_tolerance=dict(pose=.05,material_relative=.25),
                policies=policies,design_seconds=design_seconds,training_seconds=training_seconds_total,
                training_rhs_total=training_rhs_total,
                initial_candidates=candidates,distinct_candidates=len(bank),validation_noise_repetitions=repetitions,
                scope='development 2D scalar full-wave independent grids; known electronics; conditional repeated-noise statistics, not independent scene sample size')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--seeds',type=int,nargs='+',default=[7101,7102,7103,7104])
    parser.add_argument('--candidates',type=int,default=8)
    parser.add_argument('--max-nfev',type=int,default=50)
    args=parser.parse_args()
    dest=OUT/'results'/'branch_acquisition_development.json'
    records=json.loads(dest.read_text()) if dest.exists() else []
    for seed in args.seeds:
        if any(r['seed']==seed for r in records): continue
        records.append(scene(seed,args.candidates,max_nfev=args.max_nfev))
        dest.write_text(json.dumps(records,indent=2)+'\n')
