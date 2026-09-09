"""Matched Rice amplitude likelihood for the existing 3D complex-noise data.

No independent Gaussian amplitude-noise approximation. Electronic phase/delay
are omitted because |g exp(ik delay) E| is exactly insensitive to them here.
This is a strong direct likelihood baseline, NOT a reproduction of PDSOM.
"""
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import i0e,i1e
from maxwell3d import receivers,treams_field
import calibrate3d as cal

OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT.parents[1]/'delegated'/'a3_maxwell_refine'))
from calibration_fft import FFTJointModel


def rice_value_gradient(mu,j,amplitude,sigma):
    nu=np.abs(mu).reshape(-1);r=np.asarray(amplitude).reshape(-1)
    v=2*r*nu/sigma**2
    # Omit terms independent of nu; i0e avoids overflow at high SNR.
    cost=np.sum(nu**2/sigma**2 - np.log(i0e(v)) - v)
    weight=2/sigma**2*(nu-r*i1e(v)/i0e(v))
    jj=j.reshape(len(nu),-1)
    dnu=np.real(mu.reshape(-1).conj()[:,None]*jj)/np.maximum(nu[:,None],1e-300)
    return float(cost),dnu.T@weight


def test():
    rng=np.random.default_rng(61001)
    mu=rng.normal(size=20)+1j*rng.normal(size=20)
    j=rng.normal(size=(20,3))+1j*rng.normal(size=(20,3))
    amp=np.abs(mu+.1*(rng.normal(size=20)+1j*rng.normal(size=20)))
    _,g=rice_value_gradient(mu,j,amp,.15)
    h=1e-5
    fd=np.array([(rice_value_gradient(mu+h*j[:,i],j,amp,.15)[0]-rice_value_gradient(mu-h*j[:,i],j,amp,.15)[0])/(2*h) for i in range(3)])
    error=float(np.linalg.norm(fd-g)/np.linalg.norm(g))
    assert error<1e-7
    phase_j=1j*mu[:,None]
    _,gg=rice_value_gradient(mu,phase_j,amp,.15)
    assert abs(gg[0])<1e-9
    row=dict(rice_gradient_relative_error=error,phase_gradient=float(gg[0]),passed=True)
    (OUT/'results'/'intensity3d_checks.json').write_text(json.dumps(row,indent=2)+'\n')
    return row


def fit(seed,spacing=.03,maxiter=80):
    ztrue,y,_,sigma=cal.truth(seed)
    model=FFTJointModel(spacing,cal.KS)
    idx=np.array([0,1,2,3,4,6,7,8,9])
    base=np.zeros(14);base[:2]=2.
    lo=np.array([1.2,1.2]+[-.3]*3+[-.5]*4)
    hi=np.array([5.,5.]+[.3]*3+[.5]*4)
    scales=np.array([1.,1.]+[.1]*3+[.1]*4)
    cache={};t0=time.perf_counter()
    def ev(v):
        if cache.get('key')!=v.tobytes():
            z=base.copy();z[idx]=v*scales
            pred,j=model.field_jac(z)
            value,g=rice_value_gradient(pred,j[...,idx]*scales,np.abs(y),sigma)
            cache.update(key=v.tobytes(),value=value,g=g)
        return cache['value'],cache['g']
    opt=minimize(ev,base[idx]/scales,jac=True,method='L-BFGS-B',bounds=list(zip(lo/scales,hi/scales)),
                 options=dict(maxiter=maxiter,maxls=30,ftol=1e-12,gtol=1e-6))
    z=base.copy();z[idx]=opt.x*scales
    pred,j=model.field_jac(z)
    held=receivers(17,1.6)
    true=np.array([treams_field(cal.CENTERS,cal.RADII,ztrue[:2]+1j*cal.LOSS,k,held+ztrue[2:5],5)[0] for k in cal.KS])
    est=np.array([treams_field(cal.CENTERS,cal.RADII,z[:2]+1j*cal.LOSS,k,held+z[2:5],5)[0] for k in cal.KS])
    phase=np.angle(est*true.conj())
    return dict(seed=seed,spacing=spacing,method='direct_rice_intensity',estimated=z.tolist(),
                true=ztrue.tolist(),status=str(opt.message),nfev=int(opt.nfev),nit=int(opt.nit),
                pose_error_m=float(np.linalg.norm(z[2:5]-ztrue[2:5])),
                material_relative_error=float(np.linalg.norm(z[:2]-ztrue[:2])/np.linalg.norm(ztrue[:2])),
                heldout_phase_rmse_rad=float(np.sqrt(np.sum(abs(true)**2*phase**2)/np.sum(abs(true)**2))),
                heldout_field_relative_error=float(np.linalg.norm(est-true)/np.linalg.norm(true)),
                magnitude_relative_residual=float(np.linalg.norm(abs(pred)-abs(y))/np.linalg.norm(abs(y))),
                elapsed_seconds=time.perf_counter()-t0,work=[m.work for m in model.models],
                unidentifiable_and_not_fitted=['common_delay','illumination_phases'],
                scope='development same complex samples transformed to magnitudes, exact Rice likelihood, known two supports; not classical phaseless SOM; one common nominal initialization')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--spacing',type=float,default=.03)
    p.add_argument('--seeds',type=int,nargs='+',default=[6101,6102]);args=p.parse_args()
    print(json.dumps(test()),flush=True)
    dest=OUT/'results'/f'intensity3d_{args.spacing}.json'
    rows=json.loads(dest.read_text()) if dest.exists() else []
    for seed in args.seeds:
        if any(r['seed']==seed for r in rows):continue
        row=fit(seed,args.spacing);rows.append(row)
        dest.write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)
