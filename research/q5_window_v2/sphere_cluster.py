"""Transparent multiple-sphere Maxwell solver for Q5 diagnostics.

Exact analytic single-sphere Mie coefficients and numerically projected vector
spherical-wave translations. Truncation and quadrature are NOT continuum error
certificates. Time dependence exp(-i omega t), outgoing h_l^(1), T_E=-a_l,
T_M=-b_l; relative permeability is one. Lengths in meters.
"""
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve
from scipy.special import spherical_jn, spherical_yn, sph_harm_y

CENTERS=np.array([[-.06,0,0],[.055,.02,0]])
RADII=np.array([.035,.025]); LOSS=np.array([.03,.05]); K=18.
ILLUMINATIONS=[(np.array(d,float),np.array(p,float)) for d,p in [
    ([0,0,1],[1,0,0]),([0,0,1],[0,1,0]),([1,0,0],[0,1,0]),([1,0,0],[0,0,1])]]

def modes(order):
    return [(l,m,kind) for l in range(1,order+1) for m in range(-l,l+1) for kind in [0,1]]

def sphere_quadrature(ntheta=12,nphi=None):
    nphi=nphi or 2*ntheta
    u,w=leggauss(ntheta);phi=np.arange(nphi)*(2*np.pi/nphi)
    uu,pp=np.meshgrid(u,phi,indexing='ij')
    points=np.stack([np.sqrt(1-uu**2)*np.cos(pp),np.sqrt(1-uu**2)*np.sin(pp),uu],axis=-1).reshape(-1,3)
    return points,np.repeat(w,nphi)*2*np.pi/nphi

def vector_waves(points,k,order,outgoing=False):
    points=np.asarray(points,float);r=np.linalg.norm(points,axis=1)
    if np.any(r<=0):
        raise ValueError('Spherical basis evaluation at origin is unsupported')
    theta=np.arccos(np.clip(points[:,2]/r,-1,1));phi=np.arctan2(points[:,1],points[:,0])
    st=np.sin(theta);ct=np.cos(theta);sp=np.sin(phi);cp=np.cos(phi)
    if np.min(np.abs(st))<1e-10:
        raise ValueError('Use off-pole evaluation points; pole limiting formula not implemented')
    er=points/r[:,None];et=np.c_[ct*cp,ct*sp,-st];ep=np.c_[-sp,cp,np.zeros_like(st)]
    rho=k*r;out=[]
    for l in range(1,order+1):
        z=spherical_jn(l,rho).astype(complex);zp=spherical_jn(l,rho,True).astype(complex)
        if outgoing:
            z+=1j*spherical_yn(l,rho);zp+=1j*spherical_yn(l,rho,True)
        dz=zp+z/rho;norm=np.sqrt(l*(l+1))
        ys={m:sph_harm_y(l,m,theta,phi) for m in range(-l,l+1)}
        for m in range(-l,l+1):
            y=ys[m];dtheta=np.zeros_like(y)
            if m<l:dtheta+=.5*np.sqrt((l-m)*(l+m+1))*np.exp(-1j*phi)*ys[m+1]
            if m>-l:dtheta-=.5*np.sqrt((l+m)*(l-m+1))*np.exp(1j*phi)*ys[m-1]
            pt=dtheta/norm;pp=1j*m*y/(st*norm)
            mw=z[:,None]*(pp[:,None]*et-pt[:,None]*ep)
            nw=(norm*z/rho*y)[:,None]*er+dz[:,None]*(pt[:,None]*et+pp[:,None]*ep)
            out.extend([mw,nw])
    return np.stack(out,axis=-1).reshape(3*len(points),-1)

def mie_diagonal(epsilon,x,order):
    m=np.sqrt(complex(epsilon));z=m*x;diagonal=[]
    for l in range(1,order+1):
        j=spherical_jn(l,x);d=spherical_jn(l,x,True)+j/x
        h=j+1j*spherical_yn(l,x);dh=d+1j*(spherical_yn(l,x,True)+spherical_yn(l,x)/x)
        jm=spherical_jn(l,z);dm=spherical_jn(l,z,True)+jm/z
        a=(m*jm*d-j*dm)/(m*jm*dh-h*dm)
        b=(jm*d-m*j*dm)/(jm*dh-m*h*dm)
        for unused_m in range(-l,l+1):diagonal.extend([-b,-a])
    return np.array(diagonal)

class SphereCluster:
    def __init__(self,centers=CENTERS,radii=RADII,k=K,order=3,ntheta=12,projection_radius=.04):
        self.centers=np.asarray(centers,float);self.radii=np.asarray(radii,float);self.k=k;self.order=order
        self.ntheta=ntheta;self.projection_radius=projection_radius
        self.mode_list=modes(order);self.nm=len(self.mode_list);self.ns=len(self.radii)
        direction,w=sphere_quadrature(ntheta)
        local=direction*projection_radius
        regular=vector_waves(local,k,order)
        weights=np.repeat(w,3)
        norms=np.einsum('ij,i,ij->j',regular.conj(),weights,regular).real
        projector=regular.conj().T*weights[None,:]/norms[:,None]
        self.projection_gram_error=float(np.linalg.norm(projector@regular-np.eye(self.nm)))
        self.translation=np.zeros((self.ns*self.nm,self.ns*self.nm),complex)
        self.incident=np.zeros((self.ns*self.nm,4),complex)
        for i,ci in enumerate(self.centers):
            sl=slice(i*self.nm,(i+1)*self.nm)
            p=local+ci
            inc=np.stack([np.exp(1j*k*(p@d))[:,None]*pol for d,pol in ILLUMINATIONS],axis=-1).reshape(3*len(p),4)
            self.incident[sl]=projector@inc
            for j,cj in enumerate(self.centers):
                if i==j:continue
                sj=slice(j*self.nm,(j+1)*self.nm)
                if np.linalg.norm(ci-cj)<=projection_radius:
                    raise ValueError('Projection sphere includes outgoing-wave singularity')
                self.translation[sl,sj]=projector@vector_waves(p-cj,k,order,True)

    @lru_cache(maxsize=64)
    def coefficients(self,e1,e2):
        eps=np.asarray([e1,e2],complex)
        if self.ns!=2:raise ValueError('Cached Q5 adapter expects two spheres')
        diag=np.concatenate([mie_diagonal(e,self.k*r,self.order) for e,r in zip(eps,self.radii)])
        system=np.eye(len(diag))-diag[:,None]*self.translation
        return solve(system,diag[:,None]*self.incident,assume_a='gen',check_finite=False)

    def measurement_matrix(self,points):
        return np.concatenate([vector_waves(np.asarray(points)-c,self.k,self.order,True) for c in self.centers],axis=1)

    def field(self,epsilon,points):
        b=self.coefficients(complex(epsilon[0]),complex(epsilon[1]))
        return (self.measurement_matrix(points)@b).reshape(len(points),3,4)

    def forward(self,theta,points):
        return self.field(np.array(theta[:2])+1j*LOSS,np.asarray(points)+np.array([theta[2]*.001,0,0])).ravel()
