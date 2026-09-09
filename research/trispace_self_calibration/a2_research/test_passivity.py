"""Independent dense evaluation checks of the analytic bound; no estimator oracle."""
from pathlib import Path
import json
import sys
import time
import numpy as np
from scipy.special import hankel1, j0
from scipy.linalg import svdvals, solve, eigvalsh
from certificates import passive_bound, state_bound, derivative_bound

def operators(N,k):
    h=1/N
    z=-.5+(np.arange(N)+.5)*h
    X,Y=np.meshgrid(z,z,indexing='ij')
    pts=np.column_stack([X.ravel(),Y.ravel()])
    d=np.linalg.norm(pts[:,None,:]-pts[None,:,:],axis=-1)
    safe=d.copy(); np.fill_diagonal(safe,1)
    D=k*k*h*h*1j/4*hankel1(0,k*safe)
    a=h/np.sqrt(np.pi)
    np.fill_diagonal(D,k*k*(1j*np.pi*a/(2*k)*hankel1(1,k*a)-1/k**2))
    return pts,h,D,k*k*h*h/4*j0(k*d)

def main():
    rng=np.random.default_rng(50905)
    rows=[]; start=time.perf_counter()
    for N in [8,16]:
      for k in [np.pi,2*np.pi,4*np.pi]:
       pts,h,D,Q=operators(N,k)
       rx=np.column_stack([1.5*np.cos(np.arange(12)*2*np.pi/12),1.5*np.sin(np.arange(12)*2*np.pi/12)])
       S=k*k*h*h*1j/4*hankel1(0,k*np.linalg.norm(rx[:,None,:]-pts[None,:,:],axis=-1))
       u=.12+1.4*np.exp(-np.sum(pts**2,axis=1)/.1)
       for tau in [.15,.6,0.0]:
        cert=passive_bound(u,tau,k,h)
        M=np.eye(len(u))-(1+1j*tau)*u[:,None]*D
        actual_inverse=1/svdvals(M)[-1]
        algebra_error=np.linalg.norm(D.imag-(Q-cert.delta*np.eye(len(u)))) if cert.available else None
        row=dict(N=N,k=k,tau=tau,available=cert.available,reason=cert.reason,
                 alpha=cert.alpha if np.isfinite(cert.alpha) else None,
                 actual_inverse=float(actual_inverse), bound=None, matrix_identity=algebra_error)
        if cert.available:
          assert algebra_error<1e-10
          assert eigvalsh(Q)[0]>-1e-10
          assert cert.euclidean_inverse_bound>=actual_inverse*(1-1e-10)
          b=rng.normal(size=len(u))+1j*rng.normal(size=len(u))
          exact=solve(M,b)
          approx=exact+.03*(rng.normal(size=len(u))+1j*rng.normal(size=len(u)))
          z=M@approx-b
          db=state_bound(cert,z,S)
          err=np.linalg.norm(S@(approx-exact))
          Mv=-.1*(1+1j*tau)*D
          bv=.2*b
          dexact=solve(M,bv-Mv@exact)
          dapprox=dexact+.02*(rng.normal(size=len(u))+1j*rng.normal(size=len(u)))
          zv=M@dapprox-bv+Mv@approx
          d_bound=derivative_bound(cert,z,zv,Mv,S)
          d_error=np.linalg.norm(S@(dapprox-dexact))
          assert db>=err*(1-1e-10) and d_bound>=d_error*(1-1e-10)
          row.update(bound=float(cert.euclidean_inverse_bound),state_ratio=float(db/err),
                     derivative_ratio=float(d_bound/d_error))
        rows.append(row)
    bad=passive_bound(np.ones(64)*1e5,.001,4*np.pi,1/8)
    assert not bad.available
    result=dict(status='pass',cases=rows,high_contrast_refusal=bad.reason,
                seconds=time.perf_counter()-start,
                boundary='analytic model certificate; floating point diagnostics, not interval proof')
    dest=Path(__file__).parent/'results'; dest.mkdir(exist_ok=True)
    (dest/'passivity_checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'}))

if __name__=='__main__': main()
