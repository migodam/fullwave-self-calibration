"""Gaussian-between-components local-T and response-trained n-port probes.

Local component equations are retained in full grid space first.  The optional
port reduction is constructed from incident/output and interaction response
snapshots; it is not the fixed Hermite-Gaussian patch basis in A1.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .physics import VIE

@dataclass(frozen=True)
class GaussianComponent:
    amplitude: complex
    center: tuple[float,float]
    sigma: tuple[float,float]
    angle: float=0.
    def render(self,points):
        d=points-np.asarray(self.center);c,s=np.cos(self.angle),np.sin(self.angle)
        u=c*d[:,0]+s*d[:,1];v=-s*d[:,0]+c*d[:,1]
        return self.amplitude*np.exp(-.5*((u/self.sigma[0])**2+(v/self.sigma[1])**2))

def pack(components):
    """Stable real parameter vector [Re a, Im a, cx, cy, log sx, log sy, angle]."""
    return np.asarray([[q.amplitude.real,q.amplitude.imag,*q.center,np.log(q.sigma[0]),np.log(q.sigma[1]),q.angle] for q in components],float).ravel()

def unpack(theta):
    a=np.asarray(theta,float).reshape(-1,7)
    return [GaussianComponent(r[0]+1j*r[1],(r[2],r[3]),(float(np.exp(r[4])),float(np.exp(r[5])),),r[6]) for r in a]

def parameter_scales(components):
    """Positive scaling appropriate to the pack convention, not a statistical prior."""
    return np.asarray([[max(abs(q.amplitude),1e-6),max(abs(q.amplitude),1e-6),q.sigma[0],q.sigma[1],1.,1.,1.] for q in components],float).ravel()

def render(components,points): return np.sum([q.render(points) for q in components],axis=0) if components else np.zeros(len(points),complex)

def ownership_blocks(components,points,tail_mass=1e-5):
    """Finite-grid ownership for preconditioning; report tails separately."""
    vals=np.abs(np.column_stack([q.render(points) for q in components])); owner=vals.argmax(1)
    return [np.flatnonzero(owner==i) for i in range(len(components))]

def tail_audit(components,points,side):
    # Values on finite outer boundary quantify truncation, not compact support.
    edge=np.isclose(np.abs(points[:,0]),side/2-side/(2*round(np.sqrt(len(points))))) | np.isclose(np.abs(points[:,1]),side/2-side/(2*round(np.sqrt(len(points)))))
    vals=np.column_stack([np.abs(q.render(points)) for q in components])
    return {"max_grid_edge_amplitude":vals[edge].max(axis=0).tolist(),"tail_statement":"Gaussian terms are nonzero beyond the finite grid; this only audits sampled boundary tails."}

class LocalTNetwork:
    def __init__(self,vie:VIE,components):
        self.vie,self.components=vie,list(components); self.chis=[q.render(vie.points) for q in components]; self._dense_T=[None]*len(self.chis)
    def dense_local_T(self,i):
        """Exact finite-grid Ti matrix, deliberately restricted to small probes."""
        if self.vie.D.npix>512: raise ValueError("dense local T is only for bounded small-grid validation")
        if self._dense_T[i] is None:
            xi=self.chis[i]; d=self.vie.D.dense_direct()
            self._dense_T[i]=np.linalg.solve(np.eye(len(xi))-xi[:,None]*d,np.diag(xi))
        return self._dense_T[i]
    def local_response(self,i,rhs,rtol=1e-8):
        """Exact discretized local Ti rhs = (I-Xi D)^-1 Xi rhs."""
        if self.vie.D.npix<=512: return self.dense_local_T(i)@rhs
        return self.vie.solve(self.chis[i],rhs,rtol=rtol)[0]
    def full_component_identity_error(self,component_currents,incident):
        """Check ci=Ti(E+D sum_{k!=i} ck) before any port truncation."""
        component_currents=list(component_currents)
        if len(component_currents)!=len(self.chis): raise ValueError("one current array per component")
        total=sum(component_currents)
        errs=[]
        for i,ci in enumerate(component_currents):
            rhs=incident+self.vie.D.matmat(total-ci)
            ti=self.local_response(i,rhs); errs.append(float(np.linalg.norm(ci-ti)/max(np.linalg.norm(ci),1e-30)))
        return errs
    def ports(self,rank,rtol=1e-8,enrich=1):
        """Build each component's response range from E and cross-interactions."""
        Q=[]
        for i in range(len(self.chis)):
            snaps=self.local_response(i,self.vie.E,rtol)
            q=np.linalg.qr(snaps)[0][:,:min(rank,snaps.shape[1])]
            Q.append(q)
        for _ in range(enrich):
            new=[]
            for i in range(len(Q)):
                snaps=[self.local_response(i,self.vie.E,rtol)]+[self.local_response(i,self.vie.D.matmat(Q[j]),rtol) for j in range(len(Q)) if j!=i]
                u,_,_=np.linalg.svd(np.hstack(snaps),full_matrices=False);new.append(u[:,:rank])
            Q=new
        return Q
    def solve_ports(self,rank,rtol=1e-8,enrich=1):
        Q=self.ports(rank,rtol,enrich); sizes=[q.shape[1] for q in Q]; off=np.r_[0,np.cumsum(sizes)]; n=off[-1]
        A=np.eye(n,dtype=complex); b=np.zeros((n,self.vie.E.shape[1]),complex)
        V=[]
        for i,q in enumerate(Q):
            # Least-squares test coordinates on a response-trained range.
            # Q was built from exact local-T incident/interaction responses;
            # fixed HG patch columns never enter this construction.
            # Ti is retained in full grid space first, then its output and
            # input response maps are compressed.  This is distinct from an
            # HG patch Galerkin projection.
            if self.vie.D.npix<=512:
                v=q.conj().T@self.dense_local_T(i); V.append(v); b[off[i]:off[i+1]]=v@self.vie.E
            else:
                # Matrix-free Petrov block: q* Ti is applied to each RHS,
                # rather than silently replacing the local response by q*.
                V.append(None); b[off[i]:off[i+1]]=q.conj().T@self.local_response(i,self.vie.E,rtol)
        for i in range(len(Q)):
            for j in range(len(Q)):
                if i!=j:
                    coupling=V[i]@self.vie.D.matmat(Q[j]) if V[i] is not None else Q[i].conj().T@self.local_response(i,self.vie.D.matmat(Q[j]),rtol)
                    A[off[i]:off[i+1],off[j]:off[j+1]]=-coupling
        z=np.linalg.solve(A,b); j=sum(Q[i]@z[off[i]:off[i+1]] for i in range(len(Q)))
        return {"current":j,"scattered":self.vie.S@j,"total":self.vie.direct+self.vie.S@j,"Q":Q,"network_dimension":int(n)}
