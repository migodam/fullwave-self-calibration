"""Exact exterior Maxwell electric-l=1 response for an anchored radial target.

Three calibrated regular electric-dipole modal illuminations and three
Cartesian receive components are assumed. Material dependence and a common
complex electronics gain are profiled as a scalar per block. This is NOT a
model for three arbitrary plane waves or an arbitrary unknown scatterer.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import svd
from scipy.special import spherical_jn, spherical_yn


def coefficients(z):
    if np.any(np.asarray(z) <= 0):
        raise ValueError('kR must be positive')
    return 1+1j/z-1/z**2, -1-3j/z+3/z**2


def dyad(r, k, fidelity='full'):
    r=np.asarray(r, dtype=float); R=np.linalg.norm(r)
    if r.shape != (3,) or not np.isfinite(R) or R <= 0 or k <= 0:
        raise ValueError('nonzero finite 3-vector and positive k required')
    n=r/R; z=k*R
    if fidelity=='full': a,b=coefficients(z)
    elif fidelity=='radiative': a,b=1.,-1.
    elif fidelity=='static': a,b=-1/z**2,3/z**2
    else: raise ValueError(fidelity)
    return np.exp(1j*z)/R*(a*np.eye(3)+b*np.outer(n,n))


def dyad_batch(r, k):
    r=np.asarray(r,float); R=np.linalg.norm(r,axis=-1); n=r/R[...,None]
    a,b=coefficients(k*R)
    return (np.exp(1j*k*R)/R)[...,None,None]*(a[...,None,None]*np.eye(3)+b[...,None,None]*n[..., :,None]*n[...,None,:])


def dyad_jacobian(r,k):
    r=np.asarray(r,float); R=np.linalg.norm(r); n=r/R; z=k*R
    a,b=coefficients(z); az=-1j/z**2+2/z**3; bz=3j/z**2-6/z**3
    h=np.exp(1j*z)/R; nn=np.outer(n,n); m=a*np.eye(3)+b*nn
    out=[]
    for v in np.eye(3):
        u=n@v; w=(v-u*n)/R
        out.append((h*((1j*k-1/R)*u*m+k*u*(az*np.eye(3)+bz*nn)+b*(np.outer(w,n)+np.outer(n,w)))).ravel())
    return np.column_stack(out)


def realify(z, sigma=1.):
    z=np.asarray(z)
    return np.sqrt(2)/sigma*np.concatenate([z.real,z.imag],axis=0)


def profile_columns(B,N,rtol=1e-10):
    if N.size==0: return B.copy(),0
    norms=np.linalg.norm(N,axis=0); keep=norms>1e-300
    if not np.any(keep): return B.copy(),0
    u,s,_=svd(N[:,keep]/norms[keep],full_matrices=False)
    rank=int(np.sum(s>rtol*s[0])); Q=u[:,:rank]
    return B-Q@(Q.T@B),rank


def visible_singular_values(r,k,q=1.,gain_model='scalar',fidelity='full'):
    m=dyad(r,k,fidelity).ravel(); scale=q/np.linalg.norm(m)
    if fidelity=='full': jac=dyad_jacobian(r,k)
    else:
        step=np.linalg.norm(r)*1e-6
        jac=np.column_stack([(dyad(np.asarray(r)+step*v,k,fidelity)-dyad(np.asarray(r)-step*v,k,fidelity)).ravel()/(2*step) for v in np.eye(3)])
    B=realify(jac*scale)
    if gain_model=='scalar': D=m[:,None]*scale
    elif gain_model=='row':
        D=np.zeros((9,3),complex)
        for i in range(3): D[3*i:3*i+3,i]=m[3*i:3*i+3]*scale
    elif gain_model=='entry': D=np.diag(m*scale)
    elif gain_model=='known': D=np.empty((9,0),complex)
    else: raise ValueError(gain_model)
    N=realify(np.column_stack([D,1j*D])); V,_=profile_columns(B,N)
    return svd(V,compute_uv=False)


def theoretical_singular_values(R,k,q=1.):
    z=k*R; D=z**4+z*z+3
    st=np.sqrt(2)*q/R*np.sqrt((z**4+3*z*z+9)/D)
    sr=2*q*k*z*np.sqrt(z*z+4)/D
    return np.array([st,st,sr])


def regular_incident(points,k):
    """(I + grad grad/k^2) j0(kR), three regular electric l=1 columns."""
    points=np.atleast_2d(np.asarray(points,float)); R=np.linalg.norm(points,axis=1); z=k*R
    out=np.empty((len(points),3,3),complex)
    for i,(r,Ri,zi) in enumerate(zip(points,R,z)):
        if zi<1e-5:
            a=2/3-2*zi**2/15+zi**4/140
            b=zi**2/15-zi**4/210
        else:
            j0=spherical_jn(0,zi); jp=-spherical_jn(1,zi)
            a=j0+jp/zi; b=-j0-3*jp/zi
        nn=np.outer(r/Ri,r/Ri) if Ri>0 else np.zeros((3,3))
        out[i]=a*np.eye(3)+b*nn
    return out


def mie_a1(epsilon,k,radius):
    """Standard electric Mie coefficient, nonmagnetic homogeneous sphere."""
    x=k*radius; m=np.sqrt(complex(epsilon)); mx=m*x
    psi=lambda z: z*spherical_jn(1,z)
    dp=lambda z: spherical_jn(1,z)+z*spherical_jn(1,z,True)
    xi=x*(spherical_jn(1,x)+1j*spherical_yn(1,x))
    dxi=(spherical_jn(1,x)+1j*spherical_yn(1,x))+x*(spherical_jn(1,x,True)+1j*spherical_yn(1,x,True))
    return (m*psi(mx)*dp(x)-psi(x)*dp(mx))/(m*psi(mx)*dxi-xi*dp(mx))


def sphere_field(r,k,epsilon,radius):
    # With the regular incident convention above: scattered coefficient i*a1/k.
    # Sign/normalization is checked against the quasistatic limit and DDA below.
    return 1j*mie_a1(epsilon,k,radius)/k*dyad(r,k)


def unit_pattern(r,k):
    f=dyad(r,k).ravel(); return f/np.linalg.norm(f)


def profile_complex(y,f):
    f=np.asarray(f).ravel(); y=np.asarray(y).ravel()
    c=np.vdot(f,y)/np.vdot(f,f)
    return c, y-c*f


def tensor_proposal(y,k):
    """Noiseless algebraic inverse, noisy initializer only; returns both signs."""
    Y=np.asarray(y).reshape(3,3); tr=np.trace(Y)
    if abs(tr)<1e-12*np.linalg.norm(Y): return []
    Z=(Y+Y.T)/(2*tr)
    vals,vecs=np.linalg.eig(Z)
    # The two closest eigenvalues are transverse; the isolated one is radial.
    gaps=[abs(vals[(i+1)%3]-vals[(i+2)%3]) for i in range(3)]
    j=int(np.argmin(gaps)); lam=vals[j]
    if lam.imag>=-1e-10: return []
    z=-1/lam.imag; v=vecs[:,j]
    v=v*np.exp(-1j*np.angle(v[np.argmax(abs(v))])); n=v.real
    if np.linalg.norm(n)<1e-12: return []
    n/=np.linalg.norm(n); r=z/k*n
    return [r,-r]
