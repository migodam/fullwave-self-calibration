"""Small adjoint, finite-difference and local GN checks; no inverse acceptance."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.matrix_free_tangent import TangentState,MatrixFreeTangent

def main():
    rng=np.random.default_rng(2026091103);states=[]
    centers=rng.uniform(-.12,.12,(15,2));theta=rng.uniform(.02,.08,15)
    for f in (1.5e9,2.75e9):
        v=VIE(Geometry(n=8,n_tx=2,n_rx=10),f,cell_integrated=True)
        basis=np.exp(-np.sum((v.points[None,:,:]-centers[:,None,:])**2,axis=2)/(2*.035**2))*(1+.03j)
        c=theta@basis;fw=v.forward(c,rtol=1e-11)
        states.append(TangentState(v,c,basis,fw,np.arange(0,10,2),.01))
    tangent=MatrixFreeTangent(states,rtol=1e-11)
    q=rng.normal(size=15);r=rng.normal(size=2*tangent.m)
    Jq=tangent.matvec(q);Jtr=tangent.rmatvec(r)
    dot=float(abs(r@Jq-q@Jtr)/max(abs(r@Jq),1e-30))
    eps=1e-6;diff=[]
    for s in states:
        dc=q@s.partials
        plus=s.vie.forward(s.chi+eps*dc,rtol=1e-11)['scattered'][s.receivers]
        minus=s.vie.forward(s.chi-eps*dc,rtol=1e-11)['scattered'][s.receivers]
        diff.append(((plus-minus)/(2*eps*s.sigma)).ravel())
    z=np.concatenate(diff);fd=float(np.linalg.norm(np.r_[z.real,z.imag]-Jq)/np.linalg.norm(Jq))
    J=np.column_stack([tangent.matvec(q) for q in np.eye(15)])
    damping=.1;dense=np.linalg.solve(J.T@J+damping*np.eye(15),J.T@r)
    sparse,meta=tangent.damped_step(r,damping,rtol=1e-8,maxiter=100)
    error=float(np.linalg.norm(sparse-dense)/np.linalg.norm(dense))
    assert dot<1e-7 and fd<1e-7 and error<1e-6
    out=ROOT/'runs/a2/matrix_free';out.mkdir(parents=True,exist_ok=True)
    result=dict(scope='small discrete derivative/GN correctness checks, not speed or imaging acceptance',adjoint_relative=dot,finite_difference_relative=fd,GN_step_relative=error,solve=meta)
    (out/'checks.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
