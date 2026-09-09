"""Serial development comparison of standard coarse/fine calibration controls."""
import argparse
import gc
import json
import time
import numpy as np
from scipy.optimize import least_squares
from nonspherical_calibration import OUT, Model, KS, LOSS, reference_field, metrics
from nonspherical3d import reference
from maxwell3d import receivers

SCALE=np.array([1.]+[.1]*12)
LO=np.array([1.2]+[-.25]*3+[-.2]+[-.5]*4+[-1.]*4)
HI=np.array([5.]+[.25]*3+[.2]+[.5]*4+[1.]*4)


def data_for(seed):
    rng=np.random.default_rng(seed);true=np.zeros(13);true[0]=2.5
    true[1:4]=rng.normal(size=3);true[1:4]*=.09/np.linalg.norm(true[1:4])
    true[4]=.06;true[5:9]=rng.normal(0,.06,4);true[9:13]=rng.normal(0,.12,4)
    sources=[reference('ellipsoid',64,k,true[0]+1j*LOSS)[:2] for k in KS]
    mean,_=reference_field(sources,true,receivers())
    sigma=np.linalg.norm(mean)/np.sqrt(mean.size)*10**(-30/20)
    y=mean+sigma/np.sqrt(2)*(rng.normal(size=mean.shape)+1j*rng.normal(size=mean.shape))
    return true,sources,y,sigma


def real(a):return np.r_[a.real.ravel(),a.imag.ravel()]


def residual_jac(mu,j,y,sigma):
    a=j.reshape(-1,13)
    return np.sqrt(2)/sigma*real(mu-y),np.sqrt(2)/sigma*np.r_[a.real,a.imag]


def ls(evaluator,z,y,sigma,bounds=(LO,HI),max_nfev=35):
    cache={}
    def ev(v):
        if cache.get('key')!=v.tobytes():
            mu,j=evaluator(v);r,jr=residual_jac(mu,j,y,sigma)
            cache.update(key=v.tobytes(),r=r,j=jr)
        return cache['r'],cache['j']
    return least_squares(lambda v:ev(v)[0],z,jac=lambda v:ev(v)[1],bounds=bounds,
                         x_scale=SCALE,max_nfev=max_nfev,ftol=1e-9,xtol=1e-9,gtol=1e-7)


def corrected(coarse,anchor,mu_f,j_f,tangent):
    mu_c,j_c=coarse.field_jac(anchor)
    delta=mu_f-mu_c;dj=j_f-j_c if tangent else np.zeros_like(j_f)
    def ev(z):
        mu,j=coarse.field_jac(z)
        return mu+delta+np.tensordot(dj,z-anchor,axes=(-1,0)),j+dj
    return ev


def controls():
    z=np.array([2.4,.02,-.01,.03,.04]+[.02]*8)
    coarse=Model(8);fine=Model(12);mf,jf=fine.field_jac(z)
    ev=corrected(coarse,z.copy(),mf,jf,True);m,j=ev(z)
    value=float(np.linalg.norm(m-mf)/np.linalg.norm(mf))
    jac=float(np.linalg.norm(j-jf)/np.linalg.norm(jf))
    d=np.random.default_rng(82000).normal(size=13)*SCALE;h=1e-5
    fd=(ev(z+h*d)[0]-ev(z-h*d)[0])/(2*h)
    expected=np.tensordot(j,d,axes=(-1,0))
    err=float(np.linalg.norm(fd-expected)/np.linalg.norm(expected))
    assert value<1e-13 and jac<1e-13 and err<1e-7
    result=dict(anchor_value_error=value,anchor_jacobian_error=jac,surrogate_directional_error=err,passed=True)
    (OUT/'results/multifidelity_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


def run_one(seed,method,data):
    true,sources,y,sigma=data
    start=time.perf_counter();models={};events=[];warm_status=None
    base=np.zeros(13);base[0]=2.;z=base.copy()
    if method!='fine_direct':
        models['coarse']=Model(16)
        warm=ls(models['coarse'].field_jac,z,y,sigma)
        z=warm.x;warm_status=int(warm.status)
    if method=='coarse_only':
        status=f'coarse_optimizer_{warm_status}';active=models['coarse']
    else:
        models['fine']=Model(32);fine=models['fine'];active=fine
        if method in ['fine_direct','coarse_warm_fine']:
            opt=ls(fine.field_jac,z,y,sigma);z=opt.x;status=f'fine_optimizer_{opt.status}'
        else:
            radius=.25;status='outer_trial_limit'
            mu,j=fine.field_jac(z);r,jr=residual_jac(mu,j,y,sigma);loss=float(r@r/2)
            for outer in range(8):
                anchor=z.copy();ev=corrected(models['coarse'],anchor,mu,j,method=='tangent_corrected')
                opt=ls(ev,anchor,y,sigma,bounds=(np.maximum(LO,anchor-radius*SCALE),np.minimum(HI,anchor+radius*SCALE)),max_nfev=8)
                candidate=opt.x;sur_mu,_=ev(candidate)
                predicted=loss-float(opt.fun@opt.fun/2)
                step=float(np.linalg.norm((candidate-anchor)/SCALE))
                if predicted<=1e-12 or step<1e-9:
                    status='surrogate_stalled';break
                tf,tj=fine.field_jac(candidate);tr,tjr=residual_jac(tf,tj,y,sigma)
                trial_loss=float(tr@tr/2);actual=loss-trial_loss;ratio=actual/predicted
                accepted=actual>0 and ratio>.1
                events.append(dict(outer=outer,anchor_loss=loss,trial_loss=trial_loss,
                    predicted_reduction=predicted,ratio=ratio,radius=radius,accepted=accepted,
                    surrogate_field_error=float(np.linalg.norm(sur_mu-tf)/np.linalg.norm(tf)),
                    inner_nfev=int(opt.nfev),scaled_step=step))
                if accepted:
                    z,mu,j,r,jr,loss=candidate,tf,tj,tr,tjr,trial_loss
                    if actual<1e-9*(1+loss):status='fine_objective_stalled';break
                if ratio<.25:radius*=.5
                elif ratio>.75 and step>.8*radius:radius=min(.5,2*radius)
    online_seconds=time.perf_counter()-start
    work={k:[m.work.copy() for m in model.models] for k,model in models.items()}
    audit_start=time.perf_counter()
    # Evaluate all endpoints on the same N32 objective; coarse-only gets an
    # explicitly separate offline audit, not an uncharged online fine oracle.
    auditor=active if method!='coarse_only' else Model(32)
    pred,_=auditor.field_jac(z);fine_loss=float(np.linalg.norm(pred-y)**2/sigma**2)
    held=receivers(17,1.6);active.rx=held;deployable,_=active.field_jac(z)
    ht,hs=reference_field(sources,true,held)
    fac=np.exp(z[5:9]+1j*z[9:13]+1j*KS[:,None,None,None]*z[4])
    row=dict(seed=seed,method=method,status=status,warm_status=warm_status,estimated=z.tolist(),
        true=true.tolist(),online_seconds=online_seconds,offline_evaluation_seconds=time.perf_counter()-audit_start,
        work=work,fine_rhs=sum(w['rhs_columns'] for w in work.get('fine',[])),
        coarse_rhs=sum(w['rhs_columns'] for w in work.get('coarse',[])),events=events,
        final_fine_objective=fine_loss,pose_error_m=float(np.linalg.norm(z[1:4]-true[1:4])),
        material_relative_error=float(abs(z[0]-true[0])/true[0]),
        deployable_prediction=metrics(deployable,ht),structural_field=metrics(deployable/fac,hs),
        scope='serial DEVELOPMENT; includes setup and warm start; two previously inspected scenes; no global novelty/performance certificate')
    return row


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--test',action='store_true')
    p.add_argument('--seeds',type=int,nargs='+',default=[8101,8102]);args=p.parse_args()
    if args.test:controls()
    else:
        dest=OUT/'results/multifidelity_calibration.json';rows=json.loads(dest.read_text()) if dest.exists() else []
        for seed in args.seeds:
            data=data_for(seed)
            for method in ['fine_direct','coarse_only','coarse_warm_fine','value_corrected','tangent_corrected']:
                if any(r['seed']==seed and r['method']==method for r in rows):continue
                result=run_one(seed,method,data);rows.append(result)
                dest.write_text(json.dumps(rows,indent=2)+'\n')
                print(json.dumps({k:result[k] for k in ['seed','method','status','online_seconds','fine_rhs','coarse_rhs','final_fine_objective','pose_error_m','material_relative_error']}),flush=True)
                gc.collect()
