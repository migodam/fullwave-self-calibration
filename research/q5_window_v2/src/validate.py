"""Executed numerical checks. None of the floating checks is a continuum bound."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,time
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn as j, spherical_yn as y
from sphere_solver import *
from certify_modal import I,F,boundary,q_enclosure,series,decimal_bound


def run():
    out={};err=[];passive=[]
    for e in (1.5,4,2,5,1.5+.03j,4+.03j,2+.05j,5+.05j):
        for x in (.15,.2,.45,.63):
            ts=mie(e,x,5);z=np.sqrt(complex(e))*x;m=np.sqrt(complex(e))
            for t,(l,_,p) in zip(ts,modes(5)):
                jo=j(l,x);ho=jo+1j*y(l,x);do=j(l,x,True)+jo/x
                dh=do+1j*(y(l,x,True)+y(l,x)/x)
                ji=j(l,z);di=j(l,z,True)+ji/z
                if p==0:
                    c=(jo+t*ho)/ji;res=do+t*dh-m*c*di
                else:
                    c=(do+t*dh)/di;res=jo+t*ho-m*c*ji
                err.append(abs(res)/(abs(do)+abs(jo)+abs(t*dh)+abs(t*ho)))
                passive.append(-t.real-abs(t)**2)
    out['single_sphere_interface_max_relative']=float(max(err))
    out['passivity_min_absorption_channel']=float(min(passive))
    assert max(err)<1e-12 and min(passive)>-1e-14
    rec=[]
    for L in (3,5,7):
        pts,X=extraction(18,.035,L,24)
        inc=np.exp(1j*18*pts[:,2])[:,None]*[1,0,0]
        recon=(vsw(pts,18,L).reshape(-1,len(modes(L)))@(X@inc.ravel())).reshape(-1,3)
        rec.append(float(np.linalg.norm(recon-inc)/np.linalg.norm(inc)))
    out['plane_wave_truncation_relative_L3_L5_L7']=rec
    assert rec[-1]<1e-7
    sols=[SphereSolver(L,n) for L,n in [(3,18),(3,24),(4,24),(5,28)]]
    conv=[]
    for a in ([1.5,2,-2],[4,5,2],[2,3,0],[3,4,-1]):
        fs=[s.forward(a,receivers()) for s in sols]
        conv.append({'theta':a,'quadrature_18_24':float(np.linalg.norm(fs[0]-fs[1])/np.linalg.norm(fs[1])),
                     'L3_L5':float(np.linalg.norm(fs[0]-fs[3])/np.linalg.norm(fs[3])),
                     'L4_L5':float(np.linalg.norm(fs[2]-fs[3])/np.linalg.norm(fs[3]))})
    out['cluster_convergence_diagnostic']=conv
    assert max(c['quadrature_18_24'] for c in conv)<1e-10
    assert max(c['L4_L5'] for c in conv)<1e-6
    # Compare exact rational POINT intervals with a separate double Bessel path.
    dif=[]
    for x,lo,hi in [(F(1,5),F(3,2),F(4)),(F(3,20),F(2),F(5))]:
        bd=boundary(x)
        for n in range(21):
            e=lo+(hi-lo)*F(n,20);qi,qp,_=q_enclosure(I(e),x,bd)
            t=mie(float(e),float(x),1)[1];qs=(-1j*t/(1+t)).real
            dif.append(abs(qs-float((qi.lo+qi.hi)/2)))
            h=1e-4
            def amp(v): return abs(mie(v,float(x),1)[1])
            slope=(amp(float(e)+h)-amp(float(e)-h))/(2*h)
            expected=float((qp.lo+qp.hi)/2)/(1+qs*qs)**1.5
            assert abs(slope-expected)<1e-9
    out['rational_q_vs_scipy_max_absolute']=max(dif)
    assert max(dif)<1e-13
    # Proof of ratio bound, checked algebraically using exact rational numbers.
    for n in range(4,101):
        for r in (1,2):
            for d in range(3):
                ratio=F(2*n+4,2*n+2)**r*F(n+1,n+1-d)/((2*n+5)*(2*n+4))
                assert ratio<F(1,2)
    for v in [F(1,3),F(-1,3),F(2),F(-2),F(1,10**30)]:
        assert F(decimal_bound(v))<=v<=F(decimal_bound(v,upper=True))
    out['status']='passed floating diagnostics and exact arithmetic unit checks; not a C continuum certificate'
    return out

if __name__=='__main__':
    out=run();p=Path(__file__).resolve().parents[1]/'results/validation.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
