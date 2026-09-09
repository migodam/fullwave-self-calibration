"""Independent small dense vector DDA implementation, not a continuum certificate.

Clausius-Mossotti plus radiative-reaction polarizability. All pairwise multiple
scattering retained. Grid geometry and polarizability errors remain.
"""
from __future__ import annotations
from dataclasses import dataclass
from time import perf_counter
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from modal import dyad_batch, dyad_jacobian, regular_incident


def grid_target(shape='sphere',h=.025,axes=(.05,.05,.05)):
    axes=np.asarray(axes,float)
    coords=[np.arange(-int(np.ceil(a/h)),int(np.ceil(a/h))+1)*h for a in axes]
    xyz=np.stack(np.meshgrid(*coords,indexing='ij'),axis=-1).reshape(-1,3)
    if shape in ('sphere','ellipsoid','near_sphere'):
        mask=np.sum((xyz/axes)**2,axis=1)<=1+1e-10
    elif shape=='box': mask=np.all(abs(xyz)<=axes+1e-12,axis=1)
    else: raise ValueError(shape)
    return xyz[mask],h**3


@dataclass
class DDAState:
    points: np.ndarray
    moments: np.ndarray
    k: float
    setup_seconds: float
    solve_seconds: float
    relative_linear_residual: float
    material_moments: np.ndarray

    def field(self,r):
        G=self.k**2/(4*np.pi)*dyad_batch(np.asarray(r)-self.points,self.k)
        return np.einsum('nij,njl->il',G,self.moments)

    def geometry_jacobian(self,r):
        out=np.zeros((3,3,3),complex)
        for x,p in zip(self.points,self.moments):
            dj=dyad_jacobian(np.asarray(r)-x,self.k).reshape(3,3,3)
            out+=self.k**2/(4*np.pi)*np.einsum('ija,jl->ila',dj,p)
        return out.reshape(9,3)

    def material_jacobian(self,r):
        G=self.k**2/(4*np.pi)*dyad_batch(np.asarray(r)-self.points,self.k)
        return np.einsum('nij,njl->il',G,self.material_moments).ravel()


def solve(points,volume,k,epsilon):
    t=perf_counter(); points=np.asarray(points,float); n=len(points)
    alpha0=3*volume*(epsilon-1)/(epsilon+2)
    if abs(alpha0)<1e-16: raise ValueError('zero contrast')
    alpha_inv=1/alpha0-1j*k**3/(6*np.pi)
    disp=points[:,None,:]-points[None,:,:]
    # Avoid zero arguments before replacing the self block.
    disp[np.arange(n),np.arange(n)]=[1.,0.,0.]
    G=k**2/(4*np.pi)*dyad_batch(disp,k)
    G[np.arange(n),np.arange(n)]=0
    M=-G.transpose(0,2,1,3).reshape(3*n,3*n)
    M.flat[::3*n+1]+=alpha_inv
    rhs=regular_incident(points,k).reshape(3*n,3)
    setup=perf_counter()-t; t=perf_counter(); lu=lu_factor(M)
    p=lu_solve(lu,rhs)
    dalpha0=9*volume/(epsilon+2)**2
    dminv=-dalpha0/alpha0**2
    dp=lu_solve(lu,-dminv*p)
    solve_time=perf_counter()-t
    residual=np.linalg.norm(M@p-rhs)/np.linalg.norm(rhs)
    return DDAState(points,p.reshape(n,3,3),k,setup,solve_time,float(residual),dp.reshape(n,3,3))
