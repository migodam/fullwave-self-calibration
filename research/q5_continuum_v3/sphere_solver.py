"""Transparent vector-spherical-wave cluster solver using NumPy/SciPy only.

Standard Mie/T-matrix algebra, not a new inversion method. Translation and plane
wave coefficients are obtained by angular quadrature. Differences with increased
order/quadrature are diagnostics, NOT certified continuum error bounds.
Conventions: exp(-i*w*t); M=z_l X, X=Psi x rhat; N=curl(M)/k.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'
from functools import lru_cache
import numpy as np
from scipy.special import sph_harm_y, spherical_jn, spherical_yn
from scipy.linalg import solve

CENTERS = np.array([[-.06, 0, 0], [.055, .02, 0]])
RADII = np.array([.035, .025])
LOSS = np.array([.03, .05])
K = 18.
ILLUMINATIONS = [(np.array(d), np.array(p)) for d, p in
    [([0.,0,1], [1.,0,0]), ([0.,0,1], [0.,1,0]),
     ([1.,0,0], [0.,1,0]), ([1.,0,0], [0.,0,1])]]


def angular(points, order):
    """Orthonormal vector spherical harmonics in Cartesian components."""
    points = np.asarray(points)
    r = np.linalg.norm(points, axis=1)
    if np.any(r <= 0):
        raise ValueError('Expansion evaluation at origin is not implemented')
    er = points/r[:,None]
    theta = np.arccos(np.clip(er[:,2], -1, 1))
    phi = np.arctan2(er[:,1], er[:,0])
    st = np.sin(theta)
    if np.any(np.abs(st) < 1e-13):
        raise ValueError('Use off-pole points for this coordinate implementation')
    et = np.c_[np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -st]
    ep = np.c_[-np.sin(phi), np.cos(phi), np.zeros(len(r))]
    harmonics, gradients, degrees = [], [], []
    for ell in range(1, order+1):
        for m in range(-ell, ell+1):
            y, grad = sph_harm_y(ell, m, theta, phi, diff_n=1)
            psi = (grad[:,0,None]*et + grad[:,1,None]/st[:,None]*ep)/np.sqrt(ell*(ell+1))
            harmonics.append(y[:,None]*er)
            gradients.append(psi)
            degrees.append(ell)
    yr = np.stack(harmonics, axis=-1)
    psi = np.stack(gradients, axis=-1)
    xx = np.cross(psi.transpose(0,2,1), er[:,None,:]).transpose(0,2,1)
    return r, np.array(degrees), yr, psi, xx


def waves(points, k, order, outgoing=False, ang=None):
    r, ell, yr, psi, xx = angular(points, order) if ang is None else ang
    z = k*r[:,None]
    jl = spherical_jn(ell[None,:], z)
    dl = spherical_jn(ell[None,:], z, derivative=True) + jl/z
    if outgoing:
        yl = spherical_yn(ell[None,:], z)
        dl = dl + 1j*(spherical_yn(ell[None,:], z, derivative=True)+yl/z)
        jl = jl + 1j*yl
    mm = xx*jl[:,None,:]
    nn = yr*(np.sqrt(ell*(ell+1))*jl/z)[:,None,:] + psi*dl[:,None,:]
    return np.concatenate([mm,nn], axis=-1).reshape(3*len(r), -1)


def mie(eps, radius, k, order):
    """Diagonal outgoing coefficients: magnetic -b_l then electric -a_l."""
    n = np.sqrt(complex(eps)); x = k*radius; z = n*x
    ell = np.array([l for l in range(1,order+1) for m in range(-l,l+1)])
    j = spherical_jn(ell,x); d = spherical_jn(ell,x,True)+j/x
    h = j + 1j*spherical_yn(ell,x)
    dh = d + 1j*(spherical_yn(ell,x,True)+spherical_yn(ell,x)/x)
    u = spherical_jn(ell,z); v = spherical_jn(ell,z,True)+u/z
    tm = (n*v*j-d*u)/(dh*u-n*v*h)
    te = (v*j-n*u*d)/(n*u*dh-v*h)
    return np.r_[tm,te]


def quadrature(ntheta=18, nphi=36):
    c,w = np.polynomial.legendre.leggauss(ntheta)
    phi = (np.arange(nphi)+.37)*2*np.pi/nphi
    cc,pp = np.meshgrid(c,phi,indexing='ij')
    ss=np.sqrt(1-cc*cc)
    pts=np.c_[ (ss*np.cos(pp)).ravel(), (ss*np.sin(pp)).ravel(), cc.ravel() ]
    return pts,np.repeat(w,nphi)*2*np.pi/nphi


class SphereCluster:
    def __init__(self, order=3, ntheta=12, nphi=24, centers=CENTERS,
                 radii=RADII, k=K):
        self.order, self.k = order, k
        self.centers, self.radii = np.asarray(centers),np.asarray(radii)
        self.nm=order*(order+2);self.nc=2*self.nm;self.ns=len(self.radii)
        unit,w=quadrature(ntheta,nphi)
        radius=.04
        if self.ns > 1 and any(np.linalg.norm(a-b) <= radius for i,a in enumerate(self.centers) for j,b in enumerate(self.centers) if i!=j):
            raise ValueError("Projection sphere intersects another expansion center")
        local=radius*unit
        r,ell,yr,psi,xx=angular(local,order)
        jl=spherical_jn(ell,k*radius)
        dl=spherical_jn(ell,k*radius,True)+jl/(k*radius)
        # Complex-conjugate angular products, not conjugation of radial factors.
        pm=(xx.conj()*w[:,None,None]/jl[None,None,:]).reshape(-1,self.nm).T
        pn=(psi.conj()*w[:,None,None]/dl[None,None,:]).reshape(-1,self.nm).T
        project=np.r_[pm,pn]
        self.translation=np.zeros((self.ns*self.nc,self.ns*self.nc),complex)
        self.incident=np.zeros((self.ns*self.nc,4),complex)
        for i,ci in enumerate(self.centers):
            pts=local+ci; ii=slice(i*self.nc,(i+1)*self.nc)
            inc=np.stack([np.exp(1j*k*(pts@d))[:,None]*p for d,p in ILLUMINATIONS],axis=-1)
            self.incident[ii]=project@inc.reshape(-1,4)
            for j,cj in enumerate(self.centers):
                if j!=i:
                    self.translation[ii,j*self.nc:(j+1)*self.nc]=project@waves(pts-cj,k,order,True)
        self.projection_diagnostic=float(np.linalg.norm(project@waves(local,k,order)-np.eye(self.nc)))

    @lru_cache(maxsize=32)
    def coefficients(self,e1,e2):
        eps=[e1,e2] if self.ns==2 else [e1]
        diag=np.concatenate([mie(e,r,self.k,self.order) for e,r in zip(eps,self.radii)])
        return solve(np.eye(len(diag))-diag[:,None]*self.translation,
                     diag[:,None]*self.incident, assume_a='gen',check_finite=False)

    @lru_cache(maxsize=32)
    def observation(self, points_tuple):
        pts=np.array(points_tuple).reshape(-1,3)
        return np.concatenate([waves(pts-c,self.k,self.order,True) for c in self.centers],axis=1)

    def forward(self,theta,receivers,loss=LOSS):
        eps=np.array(theta[:self.ns])+1j*np.array(loss[:self.ns])
        coeff=self.coefficients(complex(eps[0]),complex(eps[1]) if self.ns==2 else 0j)
        pts=np.asarray(receivers)+np.array([theta[-1]*.001,0,0])
        # Same flattening as the inherited DipoleVIE: receiver, component, illumination.
        out=self.observation(tuple(pts.ravel()))@coeff
        return out.ravel()


def boundary_check(eps=3+.05j, radius=.035, k=K, order=4):
    """Mie tangential E/H residual, normalized separately by channel scale."""
    n=np.sqrt(eps);x=k*radius;z=n*x
    ell=np.array([l for l in range(1,order+1) for m in range(-l,l+1)])
    j=spherical_jn(ell,x);d=spherical_jn(ell,x,True)+j/x
    h=j+1j*spherical_yn(ell,x)
    dh=d+1j*(spherical_yn(ell,x,True)+spherical_yn(ell,x)/x)
    u=spherical_jn(ell,z);v=spherical_jn(ell,z,True)+u/z
    tm,te=np.split(mie(eps,radius,k,order),2)
    cm=(j+tm*h)/u;ce=(j+te*h)/(n*u)
    rm=d+tm*dh-n*cm*v;re=d+te*dh-ce*v
    return float(max(np.max(abs(rm)/(abs(d)+abs(tm*dh))),np.max(abs(re)/(abs(d)+abs(te*dh)))))
