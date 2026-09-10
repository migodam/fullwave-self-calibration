"""B1 moving head, exact point-source incident fields and vector scattering.

Caution: VSW truncation and DDA voxelization are numerical approximations.
Independent agreement is a diagnostic, not a continuum error enclosure.
"""
from __future__ import annotations
import hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from scipy.special import spherical_jn
from vsw import Cluster, angular, modes, quadrature, vector_waves
from dda import DipoleVIE, dipole_kernel

C0 = 299792458.0
FREQ = np.array([.9, 1.3, 1.7])*1e9
F0 = FREQ[1]
K = 2*np.pi*FREQ/C0
CENTERS = np.array([[-.030,-.018,0],[.037,.023,.006]])
RADII = np.array([.022,.026])
TX = np.array([[0,-.19,.075],[0,.18,.11]])
MOMENT = 1e-3
LOSS = np.array([.025,.04])
T = 8
Y_HEAD = np.linspace(-.095,.095,T)
LOWER = np.array([1.3,1.8])
UPPER = np.array([2.8,3.6])
CAL = np.array([1.8,2.2])


def material(s, frequency):
    """Exact declared narrow-band nondispersive permittivity + conductivity prior."""
    s=np.asarray(s,float)
    return 1+s+1j*s*LOSS*F0/frequency


def points(t, dz_mm=0., dx_mm=0.):
    return np.array([[dx_mm*.001, Y_HEAD[t]-.018, .105+dz_mm*.001],
                     [dx_mm*.001, Y_HEAD[t]+.018, .105+dz_mm*.001]])


def incident(points_, k):
    """E of two z electric dipoles (relative polarization moment in m^3)."""
    mat = dipole_kernel(points_,TX,k).reshape(len(points_),3,2,3)
    return mat[:,:,:,2]*MOMENT


def gain(g, tau_ps):
    return np.asarray(g)[:,None]*np.exp(2j*np.pi*(FREQ-F0)[None,:]*np.asarray(tau_ps)[:,None]*1e-12)


def source_projector(k, order, nt=18, nphi=36, radius=.035):
    dirs,w=quadrature(nt,nphi)
    _,P,C=angular(dirs,order)
    ell=np.array([v[0] for v in modes(order)])
    x=k*radius
    j=spherical_jn(ell,x)
    d=spherical_jn(ell,x,True)+j/x
    projector=np.concatenate([(C.conj()*w[:,None,None]).reshape(-1,len(ell)).T/j[:,None],
                               -(P.conj()*w[:,None,None]).reshape(-1,len(ell)).T/d[:,None]],axis=0)
    return dirs,projector


class Forward:
    """Fast vector spherical-wave output interpolated only over scalar head z."""
    def __init__(self, order=4, nt=18, nphi=36, interactions=True, nodes=11):
        self.order=order
        self.engines=[]
        self.rho=2.5
        self.nodes=np.cos(np.pi*(np.arange(nodes)+.5)/nodes)
        self.obs=[]
        self.known_direct=[]
        self.calls=0
        for f,k in zip(FREQ,K):
            eng=Cluster(CENTERS,RADII,k,[(np.array([0,1,0]),np.array([0,0,1]))],
                        order=order,nt=nt,nphi=nphi,fit_radius=.035,interactions=interactions)
            dirs,proj=source_projector(k,order,nt,nphi)
            eng.incident=np.concatenate([proj@incident(c+.035*dirs,k).reshape(-1,2)
                                        for c in CENTERS],axis=0)
            self.engines.append(eng)
            # Node-major field map, then Chebyshev coefficients on [-2.5,2.5] mm.
            mats=[]; direct=[]
            for z in self.nodes:
                pp=np.concatenate([points(t,self.rho*z) for t in range(T)])
                mm=eng.observation_matrix(pp).reshape(T,2,3,-1)
                mats.append(mm)
                direct.append(incident(pp,k).reshape(T,2,3,2))
            shp=np.shape(mats)[1:]
            cf=np.polynomial.chebyshev.chebfit(self.nodes,np.array(mats).reshape(nodes,-1),nodes-1)
            self.obs.append(cf.reshape((nodes,)+shp))
            dc=np.polynomial.chebyshev.chebfit(self.nodes,np.array(direct).reshape(nodes,-1),nodes-1)
            self.known_direct.append(dc.reshape((nodes,T,2,3,2)))

    def coefficients(self,s):
        self.calls+=1
        return [eng.coefficients(tuple(material(s,f))) for eng,f in zip(self.engines,FREQ)]

    def field(self,s,dz_mm,exact=False,dx_mm=None):
        if np.shape(dz_mm)!=(T,):raise ValueError('one normal displacement per frame')
        coeffs=self.coefficients(s)
        result=np.empty((T,3,2,3,2),complex)
        for fi,(eng,cf,dc,coeff) in enumerate(zip(self.engines,self.obs,self.known_direct,coeffs)):
            if not exact:
                vv=np.polynomial.chebyshev.chebvander(np.asarray(dz_mm)/self.rho,len(self.nodes)-1)
                mm=np.einsum('tj,jtrcn->trcn',vv,cf)
                direct=np.einsum('tj,jtrcq->trcq',vv,dc)
                result[:,fi]=np.einsum('trcn,nq->trcq',mm,coeff)+direct
            else:
                for t in range(T):
                    pp=points(t,dz_mm[t],0 if dx_mm is None else dx_mm[t])
                    mm=eng.observation_matrix(pp).reshape(2,3,-1)
                    result[t,fi]=np.einsum('rcn,nq->rcq',mm,coeff)+incident(pp,K[fi])
        return result

    def cross(self,s,dz_mm):return self.field(s,dz_mm)[:,:,:,0,:]

    def known_target(self,dz_mm):return self.cross(CAL,dz_mm)


class Independent:
    def __init__(self,spacing=.0065):
        self.engines=[];self.spacing=spacing
        for k in K:
            eng=DipoleVIE(CENTERS,RADII,spacing,k,fill_quadrature=4)
            eng.incident=incident(eng.points,k).reshape(-1,2)
            self.engines.append(eng)

    def field(self,s,dz_mm,dx_mm=None,z_slope_mm=0.):
        result=np.empty((T,3,2,3,2),complex)
        for fi,(f,k,eng) in enumerate(zip(FREQ,K,self.engines)):
            dz=np.asarray(dz_mm)+z_slope_mm*(f-F0)/(.4e9)
            pp=np.concatenate([points(t,dz[t],0 if dx_mm is None else dx_mm[t]) for t in range(T)])
            p=eng.currents(material(s,f))
            ee=(dipole_kernel(pp,eng.points,k)@p).reshape(T,2,3,2)
            result[:,fi]=ee+incident(pp,k).reshape(T,2,3,2)
        return result


class Born:
    """Exact linear-in-material discrete quadrature using *incident* E only.

    Uses volume quadrature, NOT a Clausius-Mossotti cell polarizability;
    the latter would introduce nonlinear material information into Born.
    """
    def __init__(self,spacing=.0065,nodes=11):
        from dda import voxelize
        pts,lab,fill=voxelize(CENTERS,RADII,spacing,4)
        self.rho=2.5; self.cf=[];self.dc=[];self.calls=0
        self.nodes=np.cos(np.pi*(np.arange(nodes)+.5)/nodes)
        for f,k in zip(FREQ,K):
            inc=incident(pts,k).reshape(-1,2)
            basis=[];direct=[]
            for zz in self.nodes:
                pp=np.concatenate([points(t,self.rho*zz) for t in range(T)])
                green=dipole_kernel(pp,pts,k)
                bb=[]
                for i in range(2):
                    w=np.repeat(spacing**3*fill*(lab==i)*(1+1j*LOSS[i]*F0/f),3)
                    bb.append((green@(w[:,None]*inc)).reshape(T,2,3,2))
                basis.append(np.stack(bb,axis=-1))
                direct.append(incident(pp,k).reshape(T,2,3,2))
            sh=np.shape(basis)[1:]
            self.cf.append(np.polynomial.chebyshev.chebfit(self.nodes,np.array(basis).reshape(nodes,-1),nodes-1).reshape((nodes,)+sh))
            self.dc.append(np.polynomial.chebyshev.chebfit(self.nodes,np.array(direct).reshape(nodes,-1),nodes-1).reshape((nodes,T,2,3,2)))

    def field(self,s,dz_mm):
        self.calls+=1
        result=np.empty((T,3,2,3,2),complex)
        for fi in range(3):
            for t in range(T):
                z=dz_mm[t]/self.rho
                result[t,fi]=np.polynomial.chebyshev.chebval(z,self.cf[fi][:,t])@s+np.polynomial.chebyshev.chebval(z,self.dc[fi][:,t])
        return result

    def cross(self,s,dz_mm):return self.field(s,dz_mm)[:,:,:,0,:]

    def known_target(self,dz_mm):return self.cross(CAL,dz_mm)


def source_hashes():
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}
