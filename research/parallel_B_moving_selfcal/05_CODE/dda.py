"""Independent vector Maxwell reference (Treams) and voxel dipole VIE.

Nonmagnetic homogeneous-background scattering, exp(-i omega t), SI lengths
with relative electric polarizability in volume units. The inverse discretizer
uses Clausius-Mossotti cells plus radiative correction, NOT Treams coefficients.
Finite-cell dispersion errors require convergence tests, not asserted accuracy.
"""
import time
import numpy as np
from scipy.linalg import lu_factor, lu_solve


def dipole_kernel(points, sources, k):
    """k^2 dyadic free-space Green kernel, zero diagonal for cell self term."""
    dr=np.asarray(points)[:,None,:]-np.asarray(sources)[None,:,:]
    d=np.linalg.norm(dr,axis=-1)
    safe=np.where(d>0,d,1.)
    u=dr/safe[:,:,None]
    uu=u[:,:,:,None]*u[:,:,None,:]
    eye=np.eye(3)
    g=np.exp(1j*k*safe)/(4*np.pi*safe)
    out=g[:,:,None,None]*(k*k*(eye-uu)+(1j*k/safe-1/safe**2)[:,:,None,None]*(eye-3*uu))
    out[d==0]=0
    return out.transpose(0,2,1,3).reshape(3*len(points),3*len(sources))


def voxelize(centers, radii, spacing, fill_quadrature=1):
    centers=np.asarray(centers,float); radii=np.asarray(radii,float)
    lo=np.floor(np.min(centers-radii[:,None],axis=0)/spacing).astype(int)
    hi=np.ceil(np.max(centers+radii[:,None],axis=0)/spacing).astype(int)
    xyz=np.stack(np.meshgrid(*[(np.arange(a,b)+.5)*spacing for a,b in zip(lo,hi)],indexing='ij'),axis=-1).reshape(-1,3)
    offsets=(np.arange(fill_quadrature)+.5)/fill_quadrature-.5
    offsets=np.stack(np.meshgrid(offsets,offsets,offsets,indexing='ij'),axis=-1).reshape(-1,3)*spacing
    weights=np.zeros((len(xyz),len(centers)))
    for offset in offsets:
        distance=np.linalg.norm(xyz[:,None,:]+offset-centers[None,:,:],axis=-1)
        weights+=(distance<radii)/len(offsets)
    if np.any((weights>0).sum(axis=1)>1):
        raise ValueError('overlapping spheres not supported by voxel label adapter')
    active=(weights>0).any(axis=1)
    return xyz[active],np.argmax(weights[active],axis=1),weights[active].max(axis=1)


def illuminations():
    # Two linearly independent transverse polarizations per propagation direction.
    return [(np.array([0.,0.,1.]),np.array([1.,0.,0.])),
            (np.array([0.,0.,1.]),np.array([0.,1.,0.])),
            (np.array([1.,0.,0.]),np.array([0.,1.,0.])),
            (np.array([1.,0.,0.]),np.array([0.,0.,1.]))]


def receivers(count=12, radius=1.3):
    z=1-2*(np.arange(count)+.5)/count
    a=np.arange(count)*np.pi*(3-np.sqrt(5))
    return radius*np.c_[np.sqrt(1-z*z)*np.cos(a),np.sqrt(1-z*z)*np.sin(a),z]


def treams_field(centers,radii,eps,k,rx,lmax=4):
    import treams
    spheres=[treams.TMatrix.sphere(lmax,k,r,[treams.Material(e),treams.Material()]) for r,e in zip(radii,eps)]
    tm=treams.TMatrix.cluster(spheres,np.asarray(centers)).interaction.solve()
    fields=[]; xs=[]
    for direction,pol in illuminations():
        inc=treams.plane_wave(k*direction,pol.tolist(),k0=k,material=tm.material)
        sca=tm@inc.expand(tm.basis)
        fields.append(np.asarray(sca.efield(rx)))
        xs.append([float(x) for x in tm.xs(inc)])
    return np.stack(fields,axis=-1),xs


class DipoleVIE:
    def __init__(self,centers,radii,spacing,k,fill_quadrature=1):
        self.points,self.labels,self.fill=voxelize(centers,radii,spacing,fill_quadrature)
        self.spacing=spacing; self.k=k
        self.kernel=dipole_kernel(self.points,self.points,k)
        self.incident=np.stack([np.exp(1j*k*(self.points@d))[:,None]*p for d,p in illuminations()],axis=-1).reshape(3*len(self.points),-1)
        self.work={'factorizations':0,'rhs_columns':0,'wall_seconds':0.}

    def currents(self,eps,derivatives=False):
        start=time.perf_counter()
        e=np.asarray(eps)[self.labels]
        a0=3*self.spacing**3*self.fill*(e-1)/(e+2)
        alpha=a0/(1-1j*self.k**3*a0/(6*np.pi))
        aa=np.repeat(alpha,3)
        m=np.eye(len(aa),dtype=complex)-aa[:,None]*self.kernel
        lu=lu_factor(m)
        p=lu_solve(lu,aa[:,None]*self.incident)
        self.work['factorizations']+=1
        self.work['rhs_columns']+=p.shape[1]
        if derivatives:
            total=self.incident+self.kernel@p
            dalpha=(9*self.spacing**3*self.fill/(e+2)**2)/(1-1j*self.k**3*a0/(6*np.pi))**2
            rhs=np.stack([np.repeat(dalpha*(self.labels==j),3)[:,None]*total for j in range(len(eps))],axis=-1)
            dp=lu_solve(lu,rhs.reshape(len(aa),-1)).reshape(rhs.shape)
            self.work['rhs_columns']+=rhs.shape[1]*rhs.shape[2]
            self.work['wall_seconds']+=time.perf_counter()-start
            return p,dp
        self.work['wall_seconds']+=time.perf_counter()-start
        return p

    def field(self,eps,rx):
        p=self.currents(eps)
        return (dipole_kernel(rx,self.points,self.k)@p).reshape(len(rx),3,-1)
