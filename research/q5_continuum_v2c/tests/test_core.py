import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from multipole import angular,sphere_grid,waves,radial_values,mie,Cluster,illuminations
from run_v2 import profile


def test_vsh_orthonormality():
    pts,w=sphere_grid(12)
    _,psi,negx,*_=angular(pts,5)
    b=np.concatenate([negx,psi],axis=-1).reshape(-1,70)
    gram=b.conj().T@(np.repeat(w,3)[:,None]*b)
    assert np.linalg.norm(gram-np.eye(70))<5e-12


def test_dielectric_boundary_conditions():
    for e in [1.5+.03j,4+.05j,2+0j]:
        t,jh,hh=mie(6,e,.63);n=len(t)//2
        ls=np.repeat(np.arange(1,7),2*np.arange(1,7)+1)
        m=np.sqrt(e);ji,di=radial_values(ls,m*.63)
        j,d=radial_values(ls,.63);h,dh=radial_values(ls,.63,True)
        cm=(j+t[:n]*h)/ji
        cn=(d+t[n:]*dh)/di
        assert np.max(abs(d+t[:n]*dh-m*cm*di))<2e-13
        assert np.max(abs(j+t[n:]*h-m*cn*ji))<2e-13


def test_plane_wave_reconstruction():
    c=Cluster([[0,0,0]],[.035],18,8,20)
    pts,w=sphere_grid(13);pts=pts*.027
    exact=np.stack([np.exp(1j*18*(pts@d))[:,None]*p for d,p in illuminations()],axis=-1)
    got=(waves(pts,18,8,False,.035).reshape(-1,c.nm)@c.incident).reshape(exact.shape)
    assert np.linalg.norm(got-exact)/np.linalg.norm(exact)<3e-9


def test_annulus_profile_and_gls_identity():
    rng=np.random.default_rng(912)
    f=rng.normal(size=11)+1j*rng.normal(size=11)
    sigma=.1;sr=.01;z=.97+.1j;y=z*f+.1*(rng.normal(size=11)+1j*rng.normal(size=11))
    r,g=profile(y,f,sigma,z,sr)
    cov=sigma**2*np.eye(len(f))+sr**2*np.outer(f,f.conj())
    v=y-z*f
    assert .75<abs(g)<1.25
    assert abs(r@r-2*np.vdot(v,np.linalg.solve(cov,v)).real)<1e-9
    # Deliberately force the upper annulus boundary and compare scalar search.
    rr,gg=profile(3*f,f,sigma)
    assert abs(abs(gg)-1.25)<1e-14
    assert np.linalg.norm(rr)>0


def test_order_quadrature_diagnostics():
    centers=[[-.06,0,0],[.055,.02,0]];radii=[.035,.025]
    from multipole import directions
    rx=directions(12,.6);theta=[2.8,4.1,.35]
    a=Cluster(centers,radii,18,3,14).forward(theta,rx)
    b=Cluster(centers,radii,18,3,22).forward(theta,rx)
    c=Cluster(centers,radii,18,6,20).forward(theta,rx)
    assert np.linalg.norm(a-b)/np.linalg.norm(b)<1e-9
    assert np.linalg.norm(a-c)/np.linalg.norm(c)<2e-5
