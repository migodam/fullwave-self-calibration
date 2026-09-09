"""A3 Fresnel prediction/model-adequacy extension, NOT labelled array calibration.

Frozen A2 files are imported read-only. Compare corrected plane-wave and
line-source cylindrical illumination, with optional bounded effective antenna
radii. Fitted radii lack hardware ground truth and can absorb horn-model error.
Gains and multistart selection use training views only; no test profiling.
"""
import argparse
import hashlib
import importlib.util
import json
import time
from pathlib import Path
import numpy as np
from scipy.special import hankel1,jv
from scipy.optimize import minimize

OUT=Path(__file__).resolve().parent
SOURCE=OUT.parents[1]/'delegated/a2_measured/model.py'
spec=importlib.util.spec_from_file_location('fresnel_a2_readonly',SOURCE)
mie=importlib.util.module_from_spec(spec);spec.loader.exec_module(mie)


def forward(z,frequency,kind,extra_order=0):
    center=np.array(z[:2])*.01;eps=z[2]
    de=.720+.001*(z[3] if len(z)>3 else 0.)
    dr=.760+.001*(z[4] if len(z)>4 else 0.)
    sa=np.deg2rad(np.arange(36)*10);ra=np.deg2rad(np.arange(72)*5)
    src=de*np.c_[np.cos(sa),np.sin(sa)]
    rx=dr*np.c_[np.cos(ra),np.sin(ra)]
    k=2*np.pi*frequency*1e9/mie.C0
    nmax=int(np.ceil(k*mie.A_RAD))+15+extra_order;n=np.arange(-nmax,nmax+1)
    d=rx-center;dist=np.linalg.norm(d,axis=1);angle=np.arctan2(d[:,1],d[:,0])
    receive=hankel1(n[None,:],k*dist[:,None])*np.exp(1j*angle[:,None]*n)
    s=src-center;sd=np.linalg.norm(s,axis=1);theta=np.arctan2(s[:,1],s[:,0])
    if kind=='plane':
        # Propagation direction is from source towards the cylinder.
        weights=(-1j)**n*np.exp(-1j*theta[:,None]*n)*np.exp(1j*k*sd[:,None])
    else:
        # Graf expansion of H0(k |r-source|) regular about the cylinder centre.
        weights=hankel1(n[None,:],k*sd[:,None])*np.exp(-1j*theta[:,None]*n)
    coeff=mie.mie_coefficients(k,complex(eps),mie.A_RAD,nmax)
    # H1 convention is exp(-i omega t); primary data use exp(+i omega t).
    return ((weights*coeff)@receive.T).conj()


def test():
    k=100.;s=np.array([.72,.13]);sd=np.linalg.norm(s);theta=np.arctan2(s[1],s[0])
    angles=np.arange(37)*2*np.pi/37;r=.025
    points=r*np.c_[np.cos(angles),np.sin(angles)];n=np.arange(-30,31)
    expanded=np.sum(hankel1(n,k*sd)*np.exp(-1j*n*theta)*jv(n,k*r)*np.exp(1j*angles[:,None]*n),axis=1)
    exact=hankel1(0,k*np.linalg.norm(points-s,axis=1))
    err=float(np.linalg.norm(expanded-exact)/np.linalg.norm(exact))
    z=[.2,2.6,3.3]
    a=forward(z,8,'line');b=forward(z,8,'line',extra_order=8)
    trunc=float(np.linalg.norm(a-b)/np.linalg.norm(b))
    assert err<1e-12 and trunc<1e-10
    row=dict(graf_addition_relative_error=err,cylinder_order_relative_error=trunc,
             time_convention='conjugated complete H1 model for exp(+i omega t)',passed=True)
    (OUT/'results/measured_extension_checks.json').write_text(json.dumps(row,indent=2)+'\n')
    print(json.dumps(row),flush=True)


def run():
    view,rec,freq,total,inc=mie.load_data();observed=total-inc
    rawfile=SOURCE.parent.parent/'a2_literature/data/2001_iop_17_6_301/dielTM_dec8f.exp'
    digest=hashlib.sha256(rawfile.read_bytes()).hexdigest()
    assert digest=='476cc9d1cfc98797545ab4adf69302dc5aeb45848a24cf8d7b3d222940cc79eb'
    dest=OUT/'results/measured_extension.json'
    rows=json.loads(dest.read_text()) if dest.exists() else []
    for fold in range(3):
        train=(view-1)%3!=fold;testmask=~train
        for method in ['plane','line','line_effective_radii']:
            if any(r['fold']==fold and r['method']==method for r in rows):continue
            t0=time.perf_counter();bounds=[(-6,6),(-6,6),(1.2,6)]
            if method=='line_effective_radii':bounds+=[(-3,3),(-3,3)]
            def evaluate(z,mask,gains=None):
                numerator=denominator=0.;gg={};per=[];phase_num=0.
                for f in range(1,9):
                    ix=mask&(freq==f);pred=forward(z,f,method)[view[ix]-1,rec[ix]-1];obs=observed[ix]
                    gain=np.vdot(pred,obs)/np.vdot(pred,pred) if gains is None else gains[f]
                    error=float(np.linalg.norm(gain*pred-obs)**2);energy=float(np.linalg.norm(obs)**2)
                    numerator+=error;denominator+=energy;gg[f]=gain
                    ph=float(np.sum(abs(obs)**2*np.angle(gain*pred*obs.conj())**2));phase_num+=ph
                    per.append(dict(frequency_ghz=f,normalized_squared_residual=error/energy,
                                    weighted_phase_rmse_rad=float(np.sqrt(ph/energy))))
                return numerator/denominator,gg,per,float(np.sqrt(phase_num/denominator))
            fits=[]
            for start in [[3,0,3],[0,3,3],[-3,0,3],[0,-3,3]]:
                if len(bounds)==5:start=start+[0,0]
                opt=minimize(lambda z:evaluate(z,train)[0],start,method='L-BFGS-B',bounds=bounds,
                             options=dict(maxiter=90,ftol=1e-11,gtol=1e-7))
                fits.append(dict(start=start,estimated=opt.x.tolist(),objective=float(opt.fun),
                                 status=int(opt.status),nit=int(opt.nit),nfev=int(opt.nfev)))
            best=min(fits,key=lambda r:r['objective']);z=best['estimated']
            trainerr,gains,trainper,trainphase=evaluate(z,train)
            testerr,_,testper,testphase=evaluate(z,testmask,gains)
            row=dict(fold=fold,method=method,fits=fits,best_estimated=z,
                training_normalized_squared_residual=trainerr,test_normalized_squared_residual=testerr,
                training_phase_rmse_rad=trainphase,test_phase_rmse_rad=testphase,
                train_per_frequency=trainper,test_per_frequency=testper,
                training_gains={f:[v.real,v.imag] for f,v in gains.items()},
                radii_bound_active=(len(z)==5 and any(abs(v)>2.999 for v in z[3:])),
                source_radius_m=.720+(z[3]*.001 if len(z)==5 else 0.),
                receiver_radius_m=.760+(z[4]*.001 if len(z)==5 else 0.),
                sha256=digest,wall_seconds=time.perf_counter()-t0,
                split='source view modulo 3; 24 training and 12 held-out views per fold; overlapping folds are not independent scenes',
                known_antenna_offset_truth=False,preregistered=False,
                scope='exploratory real-data model adequacy; effective radii may absorb unmodeled horn/illumination error; no calibrated hardware-position claim')
            rows.append(row);dest.write_text(json.dumps(rows,indent=2)+'\n')
            print(json.dumps({k:v for k,v in row.items() if k not in ['fits','train_per_frequency','test_per_frequency','training_gains']}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--test',action='store_true');args=parser.parse_args()
    if args.test:test()
    else:run()
