"""Large fixed-Gaussian amplitude local-step benchmark, not full inversion.

Compares full Jacobian GN, matrix-free CG, and data-SOM spectral preconditioned
CG. The latter uses standard randomized range approximation, not novel theory.
All setup/refresh costs are charged. Physical states and parameters stay fixed
while four residuals probe amortization; this is explicitly a local benchmark.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time,hashlib
import numpy as np
from scipy.sparse.linalg import LinearOperator,cg
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.physics import Geometry,VIE
from a2.matrix_free_tangent import TangentState,MatrixFreeTangent

def main():
    rng=np.random.default_rng(2026091104);out=ROOT/'runs/a2/matrix_free';out.mkdir(parents=True,exist_ok=True)
    result=[]
    for sidecount in (7,14):
        axis=np.linspace(-.12,.12,sidecount);xx,yy=np.meshgrid(axis,axis,indexing='ij');centers=np.c_[xx.ravel(),yy.ravel()];p=len(centers)
        states=[];tic=time.perf_counter()
        for frequency in (1.5e9,2.75e9):
            v=VIE(Geometry(n=64,n_tx=4,n_rx=24,aperture='half'),frequency)
            basis=np.exp(-np.sum((v.points[None,:,:]-centers[:,None,:])**2,axis=2)/(2*.025**2))*(1+.03j)*.1
            chi=np.ones(p)@basis*(49/p);fw=v.forward(chi,rtol=2e-8)
            sigma=.01*np.linalg.norm(fw['scattered'])/np.sqrt(fw['scattered'].size)
            states.append(TangentState(v,chi,basis,fw,np.arange(0,24,2),sigma))
        setup=time.perf_counter()-tic;t=MatrixFreeTangent(states,rtol=2e-8)
        # Shared damping protocol: four stochastic trace probes; no dense J used.
        tic=time.perf_counter();probes=rng.choice((-1.,1.),size=(p,4));trace=np.mean([np.linalg.norm(t.matvec(probes[:,i]))**2 for i in range(4)]);damping=.01*trace/p
        damping_seconds=time.perf_counter()-tic;damping_cost=t.cost.copy()
        residuals=[rng.normal(size=2*t.m) for _ in range(4)]
        # Dense exact-J comparator.
        td=MatrixFreeTangent(states,rtol=2e-8);tic=time.perf_counter();J=np.column_stack([td.matvec(e) for e in np.eye(p)]);H=J.T@J+damping*np.eye(p)
        solutions=np.linalg.solve(H,J.T@np.column_stack(residuals));dense_seconds=time.perf_counter()-tic
        # First-fold data-SOM range, preconditioner only. No B treated as data.
        ts=MatrixFreeTangent(states,rtol=2e-8);tic=time.perf_counter();rank=min(12,p);ell=min(rank+4,p)
        omega=rng.normal(size=(p,ell));Y=np.column_stack([ts.matvec(omega[:,i]) for i in range(ell)])
        Q=np.linalg.qr(Y,mode='reduced')[0];Bt=np.column_stack([ts.rmatvec(Q[:,i]) for i in range(Q.shape[1])]);U,s,_=np.linalg.svd(Bt,full_matrices=False);V=U[:,:rank];s=s[:rank]
        pre=LinearOperator((p,p),matvec=lambda z:z/damping+V@((1/(s*s+damping)-1/damping)*(V.T@z)),dtype=float)
        pre_setup=time.perf_counter()-tic;pre_cost=ts.cost.copy()
        methods={}
        for name,operator,M,extra,extra_cost in [('matrix_free_cg',MatrixFreeTangent(states,rtol=2e-8),None,0.,None),('data_som_preconditioned_cg',ts,pre,pre_setup,pre_cost)]:
            tic=time.perf_counter();rows=[]
            Lop=LinearOperator((p,p),matvec=lambda z:operator.rmatvec(operator.matvec(z))+damping*z,dtype=float)
            for i,res in enumerate(residuals):
                g=operator.rmatvec(res);history=[]
                q,info=cg(Lop,g,M=M,rtol=1e-5,atol=0,maxiter=100,callback=lambda q:history.append(1))
                # Dense oracle ONLY used for after-the-fact evaluation, not solve/stopping.
                rel=np.linalg.norm(q-solutions[:,i])/max(np.linalg.norm(solutions[:,i]),1e-30)
                normal=np.linalg.norm(H@q-J.T@res)/max(np.linalg.norm(J.T@res),1e-30)
                rows.append(dict(info=int(info),iterations=len(history),relative_step_error=float(rel),relative_normal_residual=float(normal)))
            methods[name]=dict(wall_seconds_including_range=extra+time.perf_counter()-tic,cost=operator.cost.copy(),range_setup_seconds=extra,range_setup_cost=extra_cost,rows=rows)
        record=dict(parameters=p,state_cells=64**2,frequencies=2,illuminations_per_frequency=4,residuals=4,
                    scope='Fixed Gaussian centers/scales; shared-state amplitude linearization, not move/split/merge or nonlinear image recovery',
                    common_state_setup_seconds=setup,damping=float(damping),common_damping_probe_seconds=damping_seconds,common_damping_probe_cost=damping_cost,
                    dense_GN=dict(wall_seconds=dense_seconds,cost=td.cost.copy(),jacobian_bytes=J.nbytes),methods=methods)
        result.append(record);(out/'benchmark.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),results=result),indent=2));print(json.dumps(record),flush=True)

if __name__=='__main__':main()
