"""Transparent numerical VSW translation solver, exp(-i wt), outgoing h1.
Not a continuum certificate. The only geometry variable is receiver translation.
"""
from functools import lru_cache
import numpy as np
from scipy.special import sph_harm_y, spherical_jn, spherical_yn
from scipy.linalg import solve

CENTERS=np.array([[-.06,0,0],[.055,.02,0]])
RADII=np.array([.035,.025]); LOSS=np.array([.03,.05]); K=18.

def illuminations():
    return [(np.array(d,float),np.array(p,float)) for d,p in
            [([0,0,1],[1,0,0]),([0,0,1],[0,1,0]),([1,0,0],[0,1,0]),([1,0,0],[0,0,1])]]

def receivers(n=12,r=.6):
    z=1-2*(np.arange(n)+.5)/n; a=np.arange(n)*np.pi*(3-np.sqrt(5))
    return r*np.c_[np.sqrt(1-z*z)*np.cos(a),np.sqrt(1-z*z)*np.sin(a),z]

def modes(L):
    return [(l,m,p) for l in range(1,L+1) for m in range(-l,l+1) for p in (0,1)]

def angular(points,L):
    """Return angular M=-r cross grad_S Y and grad_S Y; pole-free samples."""
    points=np.asarray(points); r=np.linalg.norm(points,axis=1)
    if np.any(r==0): raise ValueError('VSW basis at origin not supported')
    th=np.arccos(np.clip(points[:,2]/r,-1,1)); ph=np.arctan2(points[:,1],points[:,0])
    st=np.sin(th); ct=np.cos(th); cp=np.cos(ph); sp=np.sin(ph)
    er=points/r[:,None]; et=np.c_[ct*cp,ct*sp,-st]; ep=np.c_[-sp,cp,np.zeros(len(r))]
    if np.min(st)<1e-12: raise ValueError('Choose off-pole quadrature/receiver points')
    out=[]
    for l in range(1,L+1):
        for m in range(-l,l+1):
            y,dy=sph_harm_y(l,m,th,ph,diff_n=1)
            psi=(dy[:,0,None]*et+(1j*m*y/st)[:,None]*ep)/np.sqrt(l*(l+1))
            mm=-np.cross(er,psi)
            out.append((l,y,er,psi,mm))
    return r,out

def vsw(points,k,L,outgoing=False):
    r,angs=angular(points,L); z=k*r; cols=[]
    for l,y,er,psi,mm in angs:
        b=spherical_jn(l,z); db=spherical_jn(l,z,derivative=True)
        if outgoing:
            b=b+1j*spherical_yn(l,z); db=db+1j*spherical_yn(l,z,derivative=True)
        cols.extend([b[:,None]*mm, (np.sqrt(l*(l+1))*b/z*y)[:,None]*er+(db+b/z)[:,None]*psi])
    return np.stack(cols,axis=-1)

def quadrature(n):
    z,w=np.polynomial.legendre.leggauss(n); ph=(np.arange(2*n)+.25)*np.pi/n
    zz,pp=np.meshgrid(z,ph,indexing='ij')
    pts=np.c_[np.sqrt(1-zz.ravel()**2)*np.cos(pp.ravel()),np.sqrt(1-zz.ravel()**2)*np.sin(pp.ravel()),zz.ravel()]
    return pts,np.repeat(w,2*n)*np.pi/n

def extraction(k,R,L,n):
    pts,w=quadrature(n); _,angs=angular(pts,L); rows=[]
    for l,y,er,psi,mm in angs:
        j=spherical_jn(l,k*R); d=spherical_jn(l,k*R,derivative=True)+j/(k*R)
        rows.extend([(w[:,None]*mm.conj()/j).ravel(),(w[:,None]*psi.conj()/d).ravel()])
    return R*pts,np.stack(rows)

def mie(eps,x,L):
    """Magnetic T=-b, electric T=-a; analytic sphere interface formula."""
    m=np.sqrt(complex(eps)); out=[]
    for l in range(1,L+1):
        j=spherical_jn(l,x); d=spherical_jn(l,x,derivative=True)+j/x
        h=j+1j*spherical_yn(l,x); dh=d+1j*(spherical_yn(l,x,derivative=True)+spherical_yn(l,x)/x)
        ji=spherical_jn(l,m*x); di=spherical_jn(l,m*x,derivative=True)+ji/(m*x)
        b=(ji*d-m*j*di)/(ji*dh-m*h*di)
        a=(m*ji*d-j*di)/(m*ji*dh-h*di)
        for _ in range(2*l+1): out.extend([-b,-a])
    return np.array(out)

class SphereSolver:
    def __init__(self,L=3,nquad=18,centers=CENTERS,radii=RADII,k=K):
        self.L=L; self.nquad=nquad; self.centers=np.array(centers); self.radii=np.array(radii); self.k=k
        self.q=len(modes(L)); self.n=len(self.radii)
        self.U=np.zeros((self.n*self.q,self.n*self.q),complex)
        self.inc=np.zeros((self.n*self.q,4),complex)
        for i,c in enumerate(self.centers):
            pts,X=extraction(k,self.radii[i],L,nquad); world=pts+c; si=slice(i*self.q,(i+1)*self.q)
            self.inc[si]=X@np.stack([np.exp(1j*k*(world@d))[:,None]*p for d,p in illuminations()],axis=-1).reshape(len(world)*3,4)
            for j,cj in enumerate(self.centers):
                if i!=j:
                    self.U[si,j*self.q:(j+1)*self.q]=X@vsw(world-cj,k,L,True).reshape(len(world)*3,self.q)
        self._cache={}
    def coefficients(self,eps):
        key=tuple(np.asarray(eps,complex))
        if key in self._cache:return self._cache[key]
        t=np.concatenate([mie(e,self.k*a,self.L) for e,a in zip(eps,self.radii)])
        coeff=solve(np.eye(len(t))-t[:,None]*self.U,t[:,None]*self.inc,check_finite=False)
        if len(self._cache)>32:self._cache.clear()
        self._cache[key]=coeff;return coeff
    def observation(self,points):
        return np.concatenate([vsw(np.asarray(points)-c,self.k,self.L,True) for c in self.centers],axis=-1).reshape(len(points)*3,-1)
    def field(self,eps,points):
        return (self.observation(points)@self.coefficients(eps)).reshape(len(points),3,4)
    def forward(self,theta,rx):
        return self.field(np.asarray(theta[:2])+1j*LOSS,np.asarray(rx)+[theta[2]*.001,0,0]).ravel()
