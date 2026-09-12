"""Known homogeneous spheres: vector Maxwell multipoles with numerical translation.

Classical T-matrix algorithm, NOT a new inverse method. exp(-i omega t),
M=z_l X_lm; N=curl M/k. Nonmagnetic relative permittivities, SI lengths.
Angular translation uses quadrature and is a numerical approximation, not a
continuum error certificate. Independent of the inherited voxel DDA engine.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from scipy.special import sph_harm_y, spherical_jn, spherical_yn
from scipy.linalg import solve

CENTERS=np.array([[-.06,0,0],[.055,.02,0]])
RADII=np.array([.035,.025])
LOSS=np.array([.03,.05])
K=18.


def illuminations():
    return [(np.array(d,float),np.array(p,float)) for d,p in
            [([0,0,1],[1,0,0]),([0,0,1],[0,1,0]),([1,0,0],[0,1,0]),([1,0,0],[0,0,1])]]


def sphere_grid(ntheta, nphi):
    z,w=np.polynomial.legendre.leggauss(ntheta)
    phi=2*np.pi*np.arange(nphi)/nphi
    zz,pp=np.meshgrid(z,phi,indexing='ij')
    zz=zz.ravel();pp=pp.ravel();ss=np.sqrt(1-zz*zz)
    points=np.c_[ss*np.cos(pp),ss*np.sin(pp),zz]
    weights=np.repeat(w,nphi)*2*np.pi/nphi
    return points,weights


def modes(order):
    return [(l,m) for l in range(1,order+1) for m in range(-l,l+1)]


def basis(points,k,order,outgoing=True):
    """Rows=(point,Cartesian component); columns=(M modes,N modes)."""
    p=np.asarray(points,float);r=np.linalg.norm(p,axis=1)
    if np.any(r<=0): raise ValueError('basis origin excluded')
    theta=np.arccos(np.clip(p[:,2]/r,-1,1));phi=np.arctan2(p[:,1],p[:,0])
    st=np.sin(theta);ct=np.cos(theta);cp=np.cos(phi);sp=np.sin(phi)
    if np.any(st<1e-12): raise ValueError('polar-axis evaluation: rotate grid away from poles')
    er=p/r[:,None];et=np.c_[ct*cp,ct*sp,-st];ep=np.c_[-sp,cp,np.zeros_like(sp)]
    columnsM=[];columnsN=[]
    for l,m in modes(order):
        # scipy returns [d/dtheta,d/dphi] in diff_n=1.
        yy,grad=sph_harm_y(l,m,theta,phi,diff_n=1)
        P=(grad[:,0,None]*et+(grad[:,1]/st)[:,None]*ep)/np.sqrt(l*(l+1))
        X=np.cross(er,P)
        z=k*r;j=spherical_jn(l,z);jp=spherical_jn(l,z,True)
        if outgoing:
            j=j+1j*spherical_yn(l,z);jp=jp+1j*spherical_yn(l,z,True)
        M=j[:,None]*X
        N=-(np.sqrt(l*(l+1))*(j/z)*yy)[:,None]*er-(jp+j/z)[:,None]*P
        columnsM.append(M.ravel());columnsN.append(N.ravel())
    return np.stack(columnsM+columnsN,axis=1)


def mie(eps,x,order):
    m=np.sqrt(complex(eps));valuesM=[];valuesE=[]
    for l in range(1,order+1):
        j=spherical_jn(l,x);jp=spherical_jn(l,x,True)
        h=j+1j*spherical_yn(l,x);hp=jp+1j*spherical_yn(l,x,True)
        v=spherical_jn(l,m*x);vp=spherical_jn(l,m*x,True)
        dj=jp+j/x;dh=hp+h/x;dv=vp+v/(m*x)
        a=(m*v*dj-j*dv)/(m*v*dh-h*dv)
        b=(v*dj-m*j*dv)/(v*dh-m*h*dv)
        valuesM.extend([-b]*(2*l+1));valuesE.extend([-a]*(2*l+1))
    return np.array(valuesM+valuesE)


class SphereCluster:
    def __init__(self,order=3,ntheta=18,nphi=36,centers=CENTERS,radii=RADII,k=K):
        self.order=order;self.k=k;self.centers=np.array(centers);self.radii=np.array(radii)
        if len(self.centers)!=len(self.radii): raise ValueError('centers/radii mismatch')
        self.nm=2*order*(order+2);self.n=len(self.radii);self.ntheta=ntheta;self.nphi=nphi
        u,w=sphere_grid(ntheta,nphi);sqrtw=np.repeat(np.sqrt(w),3)
        self.U=np.zeros((self.n*self.nm,self.n*self.nm),complex)
        self.inc=np.zeros((self.n*self.nm,4),complex)
        self.projection_diagnostics=[]
        for i,(c,r) in enumerate(zip(self.centers,self.radii)):
            points=c+r*u
            V=basis(r*u,k,order,False)
            Vw=sqrtw[:,None]*V
            # Weighted projection to all retained regular multipoles.
            P=np.linalg.pinv(Vw,rcond=1e-13)*sqrtw[None,:]
            self.projection_diagnostics.append(float(np.linalg.norm(P@V-np.eye(self.nm))))
            inc=np.stack([np.exp(1j*k*(points@d))[:,None]*p for d,p in illuminations()],axis=-1).reshape(-1,4)
            si=slice(i*self.nm,(i+1)*self.nm)
            self.inc[si]=P@inc
            for j,cj in enumerate(self.centers):
                if i!=j:
                    sj=slice(j*self.nm,(j+1)*self.nm)
                    self.U[si,sj]=P@basis(points-cj,k,order,True)

    @lru_cache(maxsize=24)
    def coefficients(self,material1,material2=None):
        eps=[material1] if self.n==1 else [material1,material2]
        T=np.concatenate([mie(e,self.k*r,self.order) for e,r in zip(eps,self.radii)])
        lhs=np.eye(self.n*self.nm)-T[:,None]*self.U
        return solve(lhs,T[:,None]*self.inc,assume_a='gen',check_finite=False)

    def readout(self,points):
        return np.concatenate([basis(np.asarray(points)-c,self.k,self.order,True) for c in self.centers],axis=1)

    def field(self,eps,points):
        p=self.coefficients(*tuple(complex(e) for e in eps))
        return (self.readout(points)@p).reshape(len(points),3,4)

    def forward(self,theta,rx):
        return self.field(np.asarray(theta[:2])+1j*LOSS,np.asarray(rx)+np.array([.001*theta[2],0,0])).ravel()


class ShiftReadout:
    """Chebyshev interpolation of analytic receiver displacement (numerical only)."""
    def __init__(self,model,rx,degree=8):
        self.model=model;self.rx=np.asarray(rx);self.degree=degree
        nodes=np.cos(np.pi*(np.arange(degree+1)+.5)/(degree+1))
        values=np.array([model.readout(self.rx+np.array([.002*t,0,0])).ravel() for t in nodes])
        self.coef=np.polynomial.chebyshev.chebfit(nodes,values,degree)
        self.shape=(3*len(rx),model.n*model.nm)
    def __call__(self,theta):
        if not -2.0000001<=theta[2]<=2.0000001: raise ValueError('shift outside registered domain')
        R=np.polynomial.chebyshev.chebval(theta[2]/2,self.coef).reshape(self.shape)
        p=self.model.coefficients(*(np.asarray(theta[:2])+1j*LOSS))
        return (R@p).ravel()
