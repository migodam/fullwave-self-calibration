"""Transparent vector spherical-wave Maxwell solver (diagnostic, not certified).

exp(-i*omega*t), M=-z_l X, N=curl(M)/k, and outgoing T=(-b_l,-a_l).
Angular projection numerically translates spheres. Unknown coefficients are
scaled by tangential radial factors at the sphere boundary. No Treams import.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from scipy.special import lpmv, gammaln, spherical_jn, spherical_yn
from scipy.linalg import solve


def directions(n=12, radius=0.6):
    u = 1 - 2*(np.arange(n)+0.5)/n
    phi = np.arange(n)*np.pi*(3-np.sqrt(5))
    return radius*np.c_[np.sqrt(1-u*u)*np.cos(phi),
                        np.sqrt(1-u*u)*np.sin(phi), u]


def illuminations():
    return [(np.array(d, float), np.array(p, float)) for d,p in
            [([0,0,1],[1,0,0]), ([0,0,1],[0,1,0]),
             ([1,0,0],[0,1,0]), ([1,0,0],[0,0,1])]]


def sphere_grid(n):
    u, w = np.polynomial.legendre.leggauss(n)
    phi = np.arange(2*n)*np.pi/n
    uu, pp = np.meshgrid(u, phi, indexing='ij')
    uu, pp = uu.ravel(), pp.ravel()
    points = np.c_[np.sqrt(1-uu*uu)*np.cos(pp),
                   np.sqrt(1-uu*uu)*np.sin(pp), uu]
    return points, np.repeat(w, 2*n)*(np.pi/n)


def angular(points, order):
    """Return radial Y, Psi, -X, l, m, radii in a Cartesian basis."""
    points = np.asarray(points, float)
    rr = np.linalg.norm(points, axis=1)
    if np.any(rr <= 0):
        raise ValueError('Expansion cannot be evaluated at its origin')
    er = points/rr[:,None]
    u = np.clip(er[:,2], -1, 1)
    s = np.sqrt(np.maximum(0, 1-u*u))
    if np.any(s < 1e-13):
        raise ValueError('Choose off-pole evaluation points (angular chart)')
    phi = np.arctan2(er[:,1], er[:,0])
    ct, st = np.cos(phi), np.sin(phi)
    et = np.c_[u*ct, u*st, -s]
    ep = np.c_[-st, ct, np.zeros_like(u)]
    ys=[]; dts=[]; ls=[]; ms=[]
    for l in range(1, order+1):
        for m in range(-l,l+1):
            a=abs(m)
            norm=np.exp(0.5*(np.log((2*l+1)/(4*np.pi))+
                            gammaln(l-a+1)-gammaln(l+a+1)))
            p=lpmv(a,l,u)
            pm=lpmv(a,l-1,u) if a<=l-1 else np.zeros_like(u)
            phase=np.exp(1j*a*phi)
            y=norm*p*phase
            dt=norm*(l*u*p-(l+a)*pm)/s*phase
            if m<0:
                y=(-1)**a*np.conj(y);dt=(-1)**a*np.conj(dt)
            ys.append(y);dts.append(dt);ls.append(l);ms.append(m)
    y=np.stack(ys,axis=-1);dt=np.stack(dts,axis=-1)
    ls=np.array(ls);ms=np.array(ms);root=np.sqrt(ls*(ls+1))
    pt=dt/root;pp=1j*y*ms/(s[:,None]*root)
    psi=et[:,:,None]*pt[:,None,:]+ep[:,:,None]*pp[:,None,:]
    negx=et[:,:,None]*pp[:,None,:]-ep[:,:,None]*pt[:,None,:]
    radial=er[:,:,None]*y[:,None,:]
    return radial,psi,negx,ls,ms,rr


def radial_values(l, z, outgoing=False):
    j=spherical_jn(l,z)
    d=spherical_jn(l,z,True)+j/z
    if outgoing:
        y=spherical_yn(l,z)
        d=d+1j*(spherical_yn(l,z,True)+y/z)
        j=j+1j*y
    return j,d


def waves(points, k, order, outgoing=True, radius_scale=None):
    radial,psi,negx,ls,ms,rr=angular(points,order)
    z=k*rr[:,None]
    j,d=radial_values(ls[None,:],z,outgoing)
    m=negx*j[:,None,:]
    n=(radial*(np.sqrt(ls*(ls+1))*j/z)[:,None,:]+psi*d[:,None,:])
    if radius_scale is not None:
        js,ds=radial_values(ls,k*radius_scale,outgoing)
        m=m/js[None,None,:];n=n/ds[None,None,:]
    return np.concatenate([m,n],axis=-1)


def mie(order, eps, x):
    """Magnetic then electric outgoing coefficients, each repeated for all m."""
    ls=np.repeat(np.arange(1,order+1),2*np.arange(1,order+1)+1)
    m=np.sqrt(complex(eps))
    j,d=radial_values(ls,x)
    h,dh=radial_values(ls,x,True)
    ji,di=radial_values(ls,m*x)
    tm=-(ji*d-m*j*di)/(ji*dh-m*h*di)
    te=-(m*ji*d-j*di)/(m*ji*dh-h*di)
    return np.r_[tm,te],np.r_[j,d],np.r_[h,dh]


class Cluster:
    def __init__(self, centers, radii, k=18., order=3, quadrature=None):
        self.centers=np.asarray(centers,float);self.radii=np.asarray(radii,float)
        self.k=float(k);self.order=int(order)
        self.q=int(quadrature or max(14,order+8))
        self.nm=2*order*(order+2)
        unit,w=sphere_grid(self.q)
        _,psi,negx,_,_,_=angular(unit,order)
        basis=np.concatenate([negx,psi],axis=-1)
        projector=basis.conj().transpose(2,0,1).reshape(self.nm,-1)*np.repeat(w,3)[None,:]
        self.projector=projector;self.unit=unit;self.weights=w
        ns=len(radii);self.translation=np.zeros((ns*self.nm,ns*self.nm),complex)
        self.incident=np.zeros((ns*self.nm,4),complex)
        for i in range(ns):
            sl=slice(i*self.nm,(i+1)*self.nm)
            pts=self.centers[i]+self.radii[i]*unit
            inc=np.stack([np.exp(1j*k*(pts@d))[:,None]*p for d,p in illuminations()],axis=-1)
            self.incident[sl]=projector@inc.reshape(-1,4)
            for j in range(ns):
                if i==j:continue
                block=waves(pts-self.centers[j],k,order,True,self.radii[j])
                self.translation[sl,j*self.nm:(j+1)*self.nm]=projector@block.reshape(-1,self.nm)
        self._solutions={}

    def coefficients(self, eps, interaction=True):
        eps=tuple(complex(v) for v in eps);key=(eps,bool(interaction))
        if key in self._solutions:return self._solutions[key]
        ratios=[]
        for e,a in zip(eps,self.radii):
            t,nreg,nout=mie(self.order,e,self.k*a)
            ratios.append(t*nout/nreg)
        ratios=np.concatenate(ratios)
        if interaction:
            b=solve(np.eye(len(ratios))-self.translation*ratios[None,:],self.incident)
        else:b=self.incident.copy()
        c=ratios[:,None]*b
        if len(self._solutions)>24:self._solutions.clear()
        self._solutions[key]=(b,c)
        return b,c

    def readout(self, points):
        return np.concatenate([waves(np.asarray(points)-c,self.k,self.order,True,a)
                               for c,a in zip(self.centers,self.radii)],axis=-1).reshape(-1,len(self.radii)*self.nm)

    def field(self, eps, points, interaction=True):
        return (self.readout(points)@self.coefficients(eps,interaction)[1]).reshape(len(points),3,4)

    def forward(self, theta, points, losses=(.03,.05)):
        theta=np.asarray(theta,float)
        return self.field(theta[:2]+1j*np.array(losses),
                          np.asarray(points)+[theta[2]*.001,0,0]).ravel()

    def excitation_residual_samples(self, eps, n=32):
        """Residual of represented incident field. Numerical diagnostic ONLY."""
        b,c=self.coefficients(eps);u,w=sphere_grid(n);rows=[]
        for i,(ctr,a) in enumerate(zip(self.centers,self.radii)):
            pts=ctr+a*u
            inc=np.stack([np.exp(1j*self.k*(pts@d))[:,None]*p for d,p in illuminations()],axis=-1)
            wanted=(waves(a*u,self.k,self.order,False,a).reshape(-1,self.nm)@
                    b[i*self.nm:(i+1)*self.nm]).reshape(len(u),3,4)
            for j in range(len(self.radii)):
                if i!=j:
                    inc+=(waves(pts-self.centers[j],self.k,self.order,True,self.radii[j]).reshape(-1,self.nm)@
                          c[j*self.nm:(j+1)*self.nm]).reshape(len(u),3,4)
            r=wanted-inc
            rows.append({'surface_l2':float(np.sqrt(np.sum(w[:,None,None]*abs(r)**2))),
                         'sample_max':float(np.max(np.linalg.norm(r,axis=1)))})
        return rows
