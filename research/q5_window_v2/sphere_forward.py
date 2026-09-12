"""SciPy vector-spherical-wave two-sphere solver; numerical, NOT certified.

An offline fallback for unavailable Treams. Time convention exp(-i omega t).
M=-z_l X, N=curl(M)/k; X=rhat cross grad_S(Y)/sqrt(l(l+1)).
T_N=-a_l and T_M=-b_l. Interactions solve (I-TU)s=Ta, with U
computed by spherical quadrature. All multiple rescattering within the
truncated space is included. No new solver/algorithm priority is claimed.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from scipy.special import sph_harm_y, spherical_jn, spherical_yn
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve

CENTERS=np.array([[-.06,0.,0.],[.055,.02,0.]])
RADII=np.array([.035,.025]); LOSS=np.array([.03,.05]); K=18.

def waves():
    return [(np.array(d,float),np.array(p,float)) for d,p in
            [([0,0,1],[1,0,0]),([0,0,1],[0,1,0]),([1,0,0],[0,1,0]),([1,0,0],[0,0,1])]]

def receivers(count=12,radius=.6):
    z=1-2*(np.arange(count)+.5)/count
    a=np.arange(count)*np.pi*(3-np.sqrt(5))
    return radius*np.c_[np.sqrt(1-z*z)*np.cos(a),np.sqrt(1-z*z)*np.sin(a),z]

def angular(points,order):
    points=np.asarray(points,float); r=np.linalg.norm(points,axis=1)
    if np.any(r==0): raise ValueError('basis cannot be evaluated at center')
    er=points/r[:,None];theta=np.arccos(np.clip(er[:,2],-1,1));phi=np.arctan2(er[:,1],er[:,0])
    if np.any(np.sin(theta)<1e-12): raise ValueError('avoid polar coordinate singularity')
    et=np.c_[np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)]
    ep=np.c_[-np.sin(phi),np.cos(phi),np.zeros_like(phi)]
    ys=[];ps=[];xs=[];ls=[]
    for l in range(1,order+1):
        for m in range(-l,l+1):
            y=sph_harm_y(l,m,theta,phi)
            dy=m/np.tan(theta)*y
            if m<l: dy+=np.sqrt((l-m)*(l+m+1))*np.exp(-1j*phi)*sph_harm_y(l,m+1,theta,phi)
            p=(dy[:,None]*et+(1j*m*y/np.sin(theta))[:,None]*ep)/np.sqrt(l*(l+1))
            ys.append(y);ps.append(p);xs.append(np.cross(er,p));ls.append(l)
    return r,er,np.array(ls),np.stack(ys,axis=1),np.stack(ps,axis=2),np.stack(xs,axis=2)

def basis(points,order,k=K,outgoing=True):
    r,er,ls,Y,P,X=angular(points,order);z=k*r[:,None]
    zz=spherical_jn(ls[None,:],z); dz=spherical_jn(ls[None,:],z,True)
    if outgoing:
        zz=zz+1j*spherical_yn(ls[None,:],z)
        dz=dz+1j*spherical_yn(ls[None,:],z,True)
    M=-zz[:,None,:]*X
    N=(np.sqrt(ls*(ls+1))*zz/z*Y)[:,None,:]*er[:,:,None]+(dz+zz/z)[:,None,:]*P
    return np.concatenate([M,N],axis=2).reshape(3*len(r),-1)

def projection(order,k=K,radius=.04,nt=16,np_=32):
    u,w=leggauss(nt);p=2*np.pi*np.arange(np_)/np_
    uu,pp=np.meshgrid(u,p,indexing='ij')
    pts=radius*np.c_[(np.sqrt(1-uu*uu)*np.cos(pp)).ravel(),(np.sqrt(1-uu*uu)*np.sin(pp)).ravel(),uu.ravel()]
    weights=np.repeat(w,np_)*2*np.pi/np_
    _,_,ls,_,P,X=angular(pts,order)
    j=spherical_jn(ls,k*radius);d=spherical_jn(ls,k*radius,True)+j/(k*radius)
    mat=np.concatenate([-X.conj()/j,P.conj()/d],axis=2)*weights[:,None,None]
    return pts,mat.reshape(3*len(pts),-1).T

def mie(order,k,radius,epsilon):
    l=np.arange(1,order+1);x=k*radius;m=np.sqrt(complex(epsilon));z=m*x
    j=spherical_jn(l,x);h=j+1j*spherical_yn(l,x)
    dj=spherical_jn(l,x,True)+j/x
    dh=dj+1j*(spherical_yn(l,x,True)+spherical_yn(l,x)/x)
    ji=spherical_jn(l,z);di=spherical_jn(l,z,True)+ji/z
    tm=-(ji*dj-m*j*di)/(ji*dh-m*h*di)
    te=-(m*ji*dj-j*di)/(m*ji*dh-h*di)
    mult=2*l+1
    return np.r_[np.repeat(tm,mult),np.repeat(te,mult)]

class SphereCluster:
    def __init__(self,order=3,nt=16,np_=32,centers=CENTERS,radii=RADII,k=K):
        self.order=order;self.centers=np.array(centers);self.radii=np.array(radii);self.k=k
        pts,P=projection(order,k,nt=nt,np_=np_);self.n=P.shape[0];n=self.n
        self.U=np.zeros((n*len(centers),n*len(centers)),complex)
        inc=[]
        for i,c in enumerate(self.centers):
            p=pts+c
            inc.append(P@np.stack([(np.exp(1j*k*(p@d))[:,None]*pol).ravel() for d,pol in waves()],axis=1))
            for j,c2 in enumerate(self.centers):
                if i!=j:self.U[i*n:(i+1)*n,j*n:(j+1)*n]=P@basis(p-c2,order,k)
        self.inc=np.concatenate(inc,axis=0)
    @lru_cache(maxsize=128)
    def coefficients(self,eps):
        t=np.concatenate([mie(self.order,self.k,r,e) for r,e in zip(self.radii,eps)])
        return solve(np.eye(len(t))-t[:,None]*self.U,t[:,None]*self.inc,check_finite=False)
    def field(self,epsilon,rx):
        s=self.coefficients(tuple(complex(e) for e in epsilon))
        E=np.concatenate([basis(np.asarray(rx)-c,self.order,self.k) for c in self.centers],axis=1)@s
        return E.reshape(len(rx),3,4)
    def forward(self,theta,rx):
        return self.field(np.asarray(theta[:2])+1j*LOSS,np.asarray(rx)+np.array([theta[2]*.001,0,0])).ravel()

def self_checks():
    result={}
    pts,P=projection(4,nt=20,np_=40)
    reg=basis(pts,4,outgoing=False)
    result['projection_identity_error']=float(np.max(np.abs(P@reg-np.eye(P.shape[0]))))
    # A high-order regular plane-wave reconstruction independent of scattering.
    pts,P=projection(9,nt=24,np_=48)
    d,pol=waves()[0];a=P@(np.exp(1j*K*(pts@d))[:,None]*pol).ravel()
    check=receivers(21,.02)
    pred=basis(check,9,outgoing=False)@a
    exact=(np.exp(1j*K*(check@d))[:,None]*pol).ravel()
    result['plane_wave_relative_error']=float(np.linalg.norm(pred-exact)/np.linalg.norm(exact))
    t=mie(4,K,.035,2.5)
    result['lossless_optical_identity']=float(np.max(np.abs(t.real+abs(t)**2)))
    # Independent tangential E/H boundary residual for every electric/magnetic mode.
    errs=[]
    for eps in (1.5+.03j,4+.03j,2+.05j,5+.05j):
        m=np.sqrt(eps);x=.63;l=np.arange(1,5)
        j=spherical_jn(l,x);h=j+1j*spherical_yn(l,x)
        dj=spherical_jn(l,x,True)+j/x;dh=dj+1j*(spherical_yn(l,x,True)+spherical_yn(l,x)/x)
        ji=spherical_jn(l,m*x);di=spherical_jn(l,m*x,True)+ji/(m*x)
        tt=mie(4,1,x,eps);n=24;idx=np.cumsum(np.r_[0,2*l[:-1]+1])
        tm=tt[idx];te=tt[n+idx]
        bm=(j+tm*h)/ji;be=(dj+te*dh)/di
        errs.extend(abs(dj+tm*dh-m*bm*di)/np.maximum(abs(dj),1e-12))
        errs.extend(abs(j+te*h-m*be*ji)/np.maximum(abs(j),1e-12))
    result['boundary_relative_error']=float(max(errs))
    theta=[2.5,3.5,.37];rx=receivers()
    fields={l:SphereCluster(l,20,40).forward(theta,rx) for l in (2,3,4,5)}
    result['order_relative_to_5']={str(l):float(np.linalg.norm(v-fields[5])/np.linalg.norm(fields[5])) for l,v in fields.items()}
    fine=SphereCluster(3,24,48).forward(theta,rx)
    for nt,np_ in [(12,24),(16,32),(20,40)]:
        f=SphereCluster(3,nt,np_).forward(theta,rx)
        result[f'quadrature_{nt}_{np_}_relative_to_24_48']=float(np.linalg.norm(f-fine)/np.linalg.norm(fine))
    result['passed']=result['projection_identity_error']<1e-9 and result['plane_wave_relative_error']<1e-7 and result['boundary_relative_error']<1e-10 and result['order_relative_to_5']['3']<1e-3 and result['quadrature_16_32_relative_to_24_48']<1e-6
    result['scope']='diagnostic convergence and algebraic tests, not continuum error bounds'
    return result
if __name__=='__main__':
    import json
    print(json.dumps(self_checks(),indent=2))
