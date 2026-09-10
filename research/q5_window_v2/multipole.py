"""Vector Maxwell sphere-cluster solver; SciPy only, exp(-i omega t).

Standard Mie coefficients, vector spherical waves, and numerically projected
regular/outgoing translations. This is a transparent reference implementation,
NOT a new scattering or inversion algorithm and NOT a continuum error bound.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from scipy.special import spherical_jn, spherical_yn, sph_harm_y
from scipy.linalg import solve


def modes(order):
    return [(ell,m) for ell in range(1,order+1) for m in range(-ell,ell+1)]


def angular(directions, order):
    d=np.asarray(directions,float)
    theta=np.arccos(np.clip(d[:,2],-1,1)); phi=np.arctan2(d[:,1],d[:,0])
    st=np.sin(theta); ct=np.cos(theta)
    if np.any(st < 1e-12):
        raise ValueError('polar angular nodes unsupported; rotate quadrature')
    et=np.c_[ct*np.cos(phi),ct*np.sin(phi),-st]
    ep=np.c_[-np.sin(phi),np.cos(phi),np.zeros(len(phi))]
    Y=[]; P=[]; C=[]
    for ell,m in modes(order):
        y=sph_harm_y(ell,m,theta,phi)
        yp=sph_harm_y(ell,m+1,theta,phi) if m<ell else np.zeros_like(y)
        dt=m*ct/st*y+np.sqrt((ell-m)*(ell+m+1))*np.exp(-1j*phi)*yp
        p=(dt[:,None]*et+(1j*m*y/st)[:,None]*ep)/np.sqrt(ell*(ell+1))
        Y.append(y); P.append(p); C.append(np.cross(d,p))
    return np.stack(Y,axis=-1),np.stack(P,axis=-1),np.stack(C,axis=-1)


def vector_waves(points, center, k, order, outgoing=False):
    dr=np.asarray(points,float)-np.asarray(center,float)
    radius=np.linalg.norm(dr,axis=1)
    if np.any(radius==0):
        raise ValueError('evaluation at expansion center unsupported')
    d=dr/radius[:,None]; x=k*radius
    Y,P,C=angular(d,order)
    ell=np.array([v[0] for v in modes(order)])
    z=spherical_jn(ell[None,:],x[:,None])
    zp=spherical_jn(ell[None,:],x[:,None],True)
    if outgoing:
        z=z+1j*spherical_yn(ell[None,:],x[:,None])
        zp=zp+1j*spherical_yn(ell[None,:],x[:,None],True)
    M=C*z[:,None,:]
    N=-(np.sqrt(ell*(ell+1))[None,None,:]*(z/x[:,None])[:,None,:]
         *Y[:,None,:]*d[:,:,None]+P*(zp+z/x[:,None])[:,None,:])
    return np.concatenate([M,N],axis=-1)


def sphere_t(epsilon, size, order):
    """Diagonal scattering factors: M -> -b_l, N -> -a_l."""
    ell=np.arange(1,order+1); n=np.sqrt(complex(epsilon)); x=float(size)
    j=spherical_jn(ell,x); d=spherical_jn(ell,x,True)+j/x
    h=j+1j*spherical_yn(ell,x)
    dh=d+1j*(spherical_yn(ell,x,True)+spherical_yn(ell,x)/x)
    jm=spherical_jn(ell,n*x)
    dm=spherical_jn(ell,n*x,True)+jm/(n*x)
    a=(n*jm*d-j*dm)/(n*jm*dh-h*dm)
    b=(jm*d-n*j*dm)/(jm*dh-n*h*dm)
    ii=np.array([v[0]-1 for v in modes(order)])
    return np.r_[-b[ii],-a[ii]]


def quadrature(nt=16, nphi=32):
    z,w=np.polynomial.legendre.leggauss(nt)
    phi=2*np.pi*(np.arange(nphi)+0.37)/nphi
    zz,pp=np.meshgrid(z,phi,indexing='ij')
    rr=np.sqrt(1-zz*zz)
    dirs=np.c_[rr.ravel()*np.cos(pp.ravel()),rr.ravel()*np.sin(pp.ravel()),zz.ravel()]
    return dirs,np.repeat(w,nphi)*2*np.pi/nphi


class Cluster:
    def __init__(self,centers,radii,k,illuminations,order=3,nt=16,nphi=32,
                 fit_radius=0.04,interactions=True):
        self.centers=np.asarray(centers,float); self.radii=np.asarray(radii,float)
        self.k=float(k); self.order=int(order); self.illuminations=illuminations
        self.nm=2*len(modes(order)); nobj=len(self.radii)
        self.translation=np.zeros((nobj*self.nm,nobj*self.nm),complex)
        self.incident=np.zeros((nobj*self.nm,len(illuminations)),complex)
        dirs,w=quadrature(nt,nphi); self.quad=(nt,nphi)
        Y,P,C=angular(dirs,order)
        ell=np.array([v[0] for v in modes(order)])
        x=k*fit_radius
        j=spherical_jn(ell,x); d=spherical_jn(ell,x,True)+j/x
        projector=np.concatenate([
            (C.conj()*w[:,None,None]).reshape(-1,len(ell)).T/j[:,None],
            -(P.conj()*w[:,None,None]).reshape(-1,len(ell)).T/d[:,None]],axis=0)
        for i,center in enumerate(self.centers):
            sli=slice(i*self.nm,(i+1)*self.nm)
            points=center+fit_radius*dirs
            wave=np.stack([np.exp(1j*k*(points@direction))[:,None]*pol
                           for direction,pol in illuminations],axis=-1)
            self.incident[sli]=projector@wave.reshape(-1,len(illuminations))
            if interactions:
                for j,other in enumerate(self.centers):
                    if i==j: continue
                    slj=slice(j*self.nm,(j+1)*self.nm)
                    outgoing=vector_waves(points,other,k,order,True).reshape(-1,self.nm)
                    self.translation[sli,slj]=projector@outgoing

    @lru_cache(maxsize=32)
    def coefficients(self, epsilon_tuple):
        eps=np.asarray(epsilon_tuple,complex)
        if len(eps)!=len(self.radii):raise ValueError('one permittivity per sphere')
        t=np.concatenate([sphere_t(e,self.k*r,self.order) for e,r in zip(eps,self.radii)])
        matrix=np.eye(len(t),dtype=complex)-t[:,None]*self.translation
        return solve(matrix,t[:,None]*self.incident,check_finite=False)

    def observation_matrix(self,points):
        return np.concatenate([vector_waves(points,c,self.k,self.order,True)
                               for c in self.centers],axis=-1).reshape(-1,len(self.radii)*self.nm)

    def field(self,epsilon,points):
        return (self.observation_matrix(points)@self.coefficients(tuple(epsilon))).reshape(
            len(points),3,len(self.illuminations))


def boundary_residual(epsilon,size,order):
    """Tangential E/H continuity for unit regular waves, no cluster approximation."""
    n=np.sqrt(complex(epsilon)); x=size
    ell=np.array([v[0] for v in modes(order)])
    t=sphere_t(epsilon,x,order); nm=len(ell)
    j=spherical_jn(ell,x); dj=spherical_jn(ell,x,True)+j/x
    h=j+1j*spherical_yn(ell,x)
    dh=dj+1j*(spherical_yn(ell,x,True)+spherical_yn(ell,x)/x)
    jm=spherical_jn(ell,n*x); dm=spherical_jn(ell,n*x,True)+jm/(n*x)
    cm=(j+t[:nm]*h)/jm
    cn=(dj+t[nm:]*dh)/dm
    rm=dj+t[:nm]*dh-n*cm*dm
    rn=j+t[nm:]*h-n*cn*jm
    return float(max(np.max(np.abs(rm)),np.max(np.abs(rn))))
