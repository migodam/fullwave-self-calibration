"""Second, bounded fixed-amplitude matrix-free benchmark.

This deliberately tests whether standard randomized data-range preconditioning
ever beats a strong dense/dual GN baseline.  It is a local shared-state
linearization, not nonlinear material recovery or a new SOM theorem.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,hashlib
import numpy as np
from scipy.sparse.linalg import LinearOperator,cg
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.matrix_free_tangent import TangentState,MatrixFreeTangent
OUT=ROOT/'runs/a2/matrix_free_v2';OUT.mkdir(parents=True,exist_ok=True)

def dump(x): (OUT/'benchmark.json').write_text(json.dumps(x,indent=2))
def states_for(p,rng):
    side=int(round(np.sqrt(p)));axis=np.linspace(-.12,.12,side);xx,yy=np.meshgrid(axis,axis,indexing='ij');centers=np.c_[xx.ravel(),yy.ravel()];out=[];tic=time.perf_counter()
    for f in (1.5e9,2.75e9):
        v=VIE(Geometry(n=64,n_tx=4,n_rx=24,aperture='half'),f)
        basis=np.exp(-np.sum((v.points[None]-centers[:,None])**2,axis=2)/(2*.025**2))*(1+.03j)*.1
        chi=np.ones(p)@basis*(49/p);fw=v.forward(chi,rtol=2e-8)
        sigma=.01*np.linalg.norm(fw['scattered'])/np.sqrt(fw['scattered'].size);out.append(TangentState(v,chi,basis,fw,np.arange(0,24,2),sigma))
    return out,time.perf_counter()-tic
def dense_j(t):
    """Use cheaper of p J-columns and 2m J*-rows; never solve p^3."""
    tic=time.perf_counter(); route='columns' if t.p<=2*t.m else 'adjoint_rows'
    if route=='columns': J=np.column_stack([t.matvec(e) for e in np.eye(t.p)])
    else: J=np.vstack([t.rmatvec(e) for e in np.eye(2*t.m)])
    return J,time.perf_counter()-tic,route
def preconditioner(t,rank,rng):
    tic=time.perf_counter();ell=rank+4;omega=rng.normal(size=(t.p,ell));Y=np.column_stack([t.matvec(omega[:,i]) for i in range(ell)]);Q=np.linalg.qr(Y,mode='reduced')[0];Bt=np.column_stack([t.rmatvec(Q[:,i]) for i in range(Q.shape[1])]);U,s,_=np.linalg.svd(Bt,full_matrices=False);V=U[:,:rank];s=s[:rank]
    return V,s,time.perf_counter()-tic,t.cost.copy()
def dense_steps(J,damp,residuals):
    tic=time.perf_counter();R=np.column_stack(residuals)
    if J.shape[1]<=J.shape[0]: return np.linalg.solve(J.T@J+damp*np.eye(J.shape[1]),J.T@R),time.perf_counter()-tic,'primal_dense'
    return J.T@np.linalg.solve(J@J.T+damp*np.eye(J.shape[0]),R),time.perf_counter()-tic,'dual_woodbury'
def run_cg(t,J,damp,residuals,oracle,V=None,s=None,setup=0.,setup_cost=None):
    pre=None
    if V is not None: pre=LinearOperator((t.p,t.p),matvec=lambda z:z/damp+V@((1/(s*s+damp)-1/damp)*(V.T@z)),dtype=float)
    H=LinearOperator((t.p,t.p),matvec=lambda z:t.rmatvec(t.matvec(z))+damp*z,dtype=float)
    tic=time.perf_counter();rows=[]
    for i,r in enumerate(residuals):
        g=t.rmatvec(r);hist=[];q,info=cg(H,g,M=pre,rtol=1e-5,atol=0,maxiter=100,callback=lambda _:hist.append(1))
        nr=np.linalg.norm(J.T@(J@q-r)+damp*q)/max(np.linalg.norm(J.T@r),1e-30)
        rows.append({'info':int(info),'iterations':len(hist),'iteration_cap_reached':bool(info==100),'relative_step_error':float(np.linalg.norm(q-oracle[:,i])/max(np.linalg.norm(oracle[:,i]),1e-30)),'explicit_J_normal_residual':float(nr)})
    return {'wall_seconds_including_range':setup+time.perf_counter()-tic,'cost':t.cost.copy(),'range_setup_seconds':setup,'range_setup_cost':setup_cost,'rows':rows}
def main():
    rng=np.random.default_rng(2026091105);data={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Fixed Gaussian centers/scales and four shared-state residuals. Standard randomized data-range preconditioning; no nonlinear recovery, no parameter birth/death, no novel theorem. A CG iteration cap is reported as a cap, never convergence.','results':[]};dump(data)
    for p in (784,3136):
        states,setup=states_for(p,rng);base=MatrixFreeTangent(states,rtol=2e-8);probes=rng.choice((-1.,1.),size=(p,4));tic=time.perf_counter();trace=np.mean([np.linalg.norm(base.matvec(probes[:,i]))**2 for i in range(4)]);damp=.01*trace/p;probe_seconds=time.perf_counter()-tic;probe_cost=base.cost.copy();resids=[rng.normal(size=2*base.m) for _ in range(4)]
        dense=MatrixFreeTangent(states,rtol=2e-8);J,dense_seconds,route=dense_j(dense);oracle,solve_seconds,solve_route=dense_steps(J,damp,resids)
        record={'parameters':p,'state_cells':4096,'frequencies':2,'illuminations_per_frequency':4,'receivers_per_frequency':12,'residuals':4,'common_state_setup_seconds':setup,'damping':float(damp),'damping_probe_seconds':probe_seconds,'damping_probe_cost':probe_cost,'dense_GN':{'wall_seconds':dense_seconds+solve_seconds,'assembly_seconds':dense_seconds,'solve_seconds':solve_seconds,'assembly_route':route,'solution_route':solve_route,'cost':dense.cost.copy(),'jacobian_bytes':J.nbytes},'methods':{}}
        plain=MatrixFreeTangent(states,rtol=2e-8);record['methods']['plain_cg']=run_cg(plain,J,damp,resids,oracle);data['results'].append(record);dump(data)
        for rank in (12,32,64):
            spec=MatrixFreeTangent(states,rtol=2e-8);V,s,st,sc=preconditioner(spec,rank,rng);record['methods'][f'data_som_rank_{rank}']=run_cg(spec,J,damp,resids,oracle,V,s,st,sc);dump(data)
        print(json.dumps(record),flush=True)
    fig,ax=plt.subplots(figsize=(6,4));
    for rec in data['results']:
        names=['dense_GN','plain_cg']+[f'data_som_rank_{r}' for r in (12,32,64)];times=[rec['dense_GN']['wall_seconds']]+[rec['methods'][x]['wall_seconds_including_range'] for x in names[1:]];ax.plot(names,times,'o-',label=f"p={rec['parameters']}")
    ax.set(yscale='log',ylabel='wall seconds including setup',title='Fixed-amplitude local GN cost');ax.tick_params(axis='x',rotation=25);ax.legend();fig.tight_layout();fig.savefig(OUT/'time_vs_parameter_rank.png',dpi=170)
if __name__=='__main__':main()
