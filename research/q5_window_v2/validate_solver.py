"""Numerical checks only. These are NOT rigorous continuum error bounds."""
import json,time
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn,spherical_yn
from sphere_cluster import *

HERE=Path(__file__).resolve().parent

def relative(a,b):return float(np.linalg.norm(a-b)/np.linalg.norm(b))

def validate():
    start=time.perf_counter();rows={}
    dirs,w=sphere_quadrature(18);r=.04;points=dirs*r
    basis=vector_waves(points,K,9)
    ww=np.repeat(w,3)
    norms=np.einsum('ij,i,ij->j',basis.conj(),ww,basis).real
    inc=np.exp(1j*K*points[:,2])[:,None]*np.array([1,0,0])
    coeff=(basis.conj().T@(ww*inc.ravel()))/norms
    rows['plane_wave_reconstruction_lmax9']=relative(basis@coeff,inc.ravel())
    rows['basis_gram_error_lmax9']=float(np.linalg.norm((basis.conj().T*ww)@basis/np.sqrt(norms[:,None]*norms[None,:])-np.eye(len(norms))))
    # Every mode satisfies both E/H tangential interface equations.
    errs=[]
    for e in [1.5+.03j,2.4+.05j,4+.03j,5+.05j]:
        for x in [.2,.45,.63]:
            md=mie_diagonal(e,x,5);m=np.sqrt(e);offset=0
            for l in range(1,6):
                j=spherical_jn(l,x);d=spherical_jn(l,x,True)+j/x
                h=j+1j*spherical_yn(l,x);dh=d+1j*(spherical_yn(l,x,True)+spherical_yn(l,x)/x)
                jm=spherical_jn(l,m*x);dm=spherical_jn(l,m*x,True)+jm/(m*x)
                tm,te=md[offset:offset+2]
                cm=(j+tm*h)/jm;ce=(d+te*dh)/dm
                errs += [abs(d+tm*dh-m*cm*dm)/max(abs(d),1e-30),abs(j+te*h-m*ce*jm)/max(abs(j),1e-30)]
                offset+=2*(2*l+1)
    rows['max_relative_tangential_boundary_residual']=float(max(errs))
    # Optical theorem for lossless T coefficient in this convention.
    ts=mie_diagonal(3.0,.2,4)
    rows['lossless_optical_theorem_absolute_residual']=float(np.max(np.abs(ts.real+abs(ts)**2)))
    rx=sphere_quadrature(5,7)[0]*.6
    theta=np.array([2.4,3.7,.3]);f={}
    for l,n in [(3,10),(3,14),(4,14),(5,18)]:
        model=SphereCluster(order=l,ntheta=n);f[(l,n)]=model.forward(theta,rx)
    rows['translation_quadrature_l3_n10_vs14']=relative(f[(3,10)],f[(3,14)])
    rows['truncation_l3_vs4']=relative(f[(3,14)],f[(4,14)])
    rows['truncation_l4_vs5']=relative(f[(4,14)],f[(5,18)])
    rows['wall_seconds']=time.perf_counter()-start
    rows['status']='passed_numerical_checks' if (rows['plane_wave_reconstruction_lmax9']<1e-8 and rows['max_relative_tangential_boundary_residual']<1e-10 and rows['translation_quadrature_l3_n10_vs14']<1e-6) else 'failed'
    return rows

if __name__=='__main__':
    result=validate();path=HERE/'results/solver_checks.json'
    if path.exists():raise FileExistsError(path)
    path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
