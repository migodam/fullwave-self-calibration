"""Exact discrete Jv/J* v for real material coordinates and complex fields.

This is a strong matrix-free GN baseline and a reusable SOM building block,
not a new inverse-scattering method. It never assembles all state derivatives.
"""
from dataclasses import dataclass
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres, cg

@dataclass
class TangentState:
    vie: object
    chi: np.ndarray
    partials: np.ndarray  # [p,N], already physically scaled; complex permitted
    forward: dict
    receivers: np.ndarray
    sigma: float = 1.

class MatrixFreeTangent:
    def __init__(self, states, rtol=1e-8):
        self.states=list(states);self.rtol=rtol
        self.p=self.states[0].partials.shape[0]
        self.sizes=[s.vie.S[s.receivers].shape[0]*s.vie.E.shape[1] for s in self.states]
        self.offsets=np.r_[0,np.cumsum(self.sizes)];self.m=int(self.offsets[-1])
        self.total=[s.vie.E+s.vie.D.matmat(s.forward['current']) for s in self.states]
        self.cost={'tangent_rhs':0,'adjoint_rhs':0,'Jv_calls':0,'Jstar_calls':0}
        self.operator=LinearOperator((2*self.m,self.p),matvec=self.matvec,rmatvec=self.rmatvec,dtype=float)

    def matvec(self,q):
        q=np.asarray(q).reshape(self.p);result=[];self.cost['Jv_calls']+=1
        for s,u in zip(self.states,self.total):
            dc=q@s.partials
            dj,_=s.vie.solve_linear(s.chi,dc[:,None]*u,rtol=self.rtol)
            result.append((s.vie.S[s.receivers]@dj/s.sigma).ravel())
            self.cost['tangent_rhs']+=u.shape[1]
        z=np.concatenate(result)
        return np.r_[z.real,z.imag]

    def rmatvec(self,r):
        r=np.asarray(r).reshape(2*self.m);z=r[:self.m]+1j*r[self.m:]
        gradient=np.zeros(self.p);self.cost['Jstar_calls']+=1
        for index,(s,u) in enumerate(zip(self.states,self.total)):
            residual=z[self.offsets[index]:self.offsets[index+1]].reshape(-1,u.shape[1])
            rhs=s.vie.S[s.receivers].conj().T@residual/s.sigma
            adjoint=np.empty_like(rhs);M=s.vie.system(s.chi).H
            for t in range(rhs.shape[1]):
                adjoint[:,t],info=gmres(M,rhs[:,t],rtol=self.rtol,atol=0,maxiter=300)
                if info:raise RuntimeError(f'adjoint GMRES did not converge: {info}')
            gradient+=np.real(s.partials.conj()@np.sum(u.conj()*adjoint,axis=1))
            self.cost['adjoint_rhs']+=u.shape[1]
        return gradient

    def damped_step(self,residual,damping=.01,rtol=1e-3,maxiter=50):
        """CG baseline with a reported normal-equation residual and gap bound.

        Bounds assume exact operator arithmetic. State-solve errors are reported
        separately by the caller and require A2 bounds if certification is claimed.
        """
        started=time.perf_counter();before=self.cost.copy();g=self.rmatvec(residual);history=[]
        H=LinearOperator((self.p,self.p),matvec=lambda q:self.rmatvec(self.matvec(q))+damping*q,dtype=float)
        q,info=cg(H,g,rtol=rtol,atol=0,maxiter=maxiter,callback=lambda q:history.append(float(np.linalg.norm(q))))
        s=g-H@q
        return q,dict(status=int(info),iterations=len(history),normal_residual=float(np.linalg.norm(s)),
                      ideal_local_gap_bound=float(np.linalg.norm(s)**2/(2*damping)),
                      wall_seconds=time.perf_counter()-started,cost={k:self.cost[k]-before[k] for k in self.cost})
