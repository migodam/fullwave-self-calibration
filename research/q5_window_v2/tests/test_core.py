"""Proof-kernel and solver regression tests; diagnostics are not continuum proofs."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from fractions import Fraction as F
import numpy as np
from certify_modal import I,q_enclosure,boundary,decimal_bound,series
from certify_ratio import compact,magnetic
from continuum_bound import constants
from finite_geometry import chord_squared,covered_chord_lower_squared
from window_cover import CertificateAssumptions,finite_pair_cover,class_C_status
from sphere_solver import mie
from recovery import profile

def test_outward_decimal_and_compaction():
    for q in [F(1,3),F(-1,3),F(1,10**40),F(-1,10**40),F(7,8)]:
        assert F(decimal_bound(q))<=q<=F(decimal_bound(q,upper=True))
        assert compact(I(q)).lo<=q<=compact(I(q)).hi

def test_interval_arithmetic():
    a,b=I(-2,3),I(F(1,3),4)
    for x in [a.lo,F(0),a.hi]:
        for y in [b.lo,F(1),b.hi]:
            for v,z in [(a+b,x+y),(a*b,x*y),(a/b,x/y)]:assert v.lo<=z<=v.hi
    import pytest
    with pytest.raises(ValueError):a/I(-1,1)

def test_series_tail_ratios():
    for r in (1,2):
        for d in range(3):
            for n in range(16,100):
                ratio=F((2*n+4)**r,(2*n+2)**r)*F(n+1,n+1-d)/((2*n+5)*(2*n+4))
                assert ratio<F(1,2)

def test_mie_and_ratio_diagnostics():
    for x in [F(1,5),F(3,20)]:
        for e in [F(3,2),F(2),F(3),F(4)]:
            a,ap,_=q_enclosure(I(e),x,boundary(x));b,bp=magnetic(I(e),x,boundary(x))
            t=mie(float(e),float(x),1)
            qa=1j*(-t[1])/(1+t[1]);qb=1j*(-t[0])/(1+t[0])
            assert abs(qa.real-float(a.lo))<1e-13 and abs(qa.imag)<1e-13
            assert abs(qb.real-float(b.lo))<1e-13 and abs(qb.imag)<1e-13
            da,db=float(ap.lo),float(bp.lo);aa,bb=float(a.lo),float(b.lo)
            deriv=((db*aa-bb*da)/aa**2)*(1+aa*bb)/(1+bb**2)
            deriv+=(bb/aa)*(da*bb+aa*db)/(1+bb**2)
            deriv-=(bb/aa)*(1+aa*bb)*2*bb*db/(1+bb**2)**2
            h=1e-4
            def ratio(eps):
                mm=mie(eps,float(x),1);return (mm[0]/mm[1]).real
            numerical=(ratio(float(e)+h)-ratio(float(e)-h))/(2*h)
            assert abs(deriv-numerical)<1e-8

def test_finite_dipole_chord_diagnostic():
    def v(z):return np.array([np.sqrt(2)*(z*z+1j*z-1),2*(1-1j*z)])
    for a,b in [(F(1,4),F(3,4)),(F(1),F(11,10)),(F(2),F(3))]:
        va,vb=v(float(a)),v(float(b))
        numerical=abs(np.linalg.det(np.stack([va,vb])))**2/(np.vdot(va,va).real*np.vdot(vb,vb).real)
        assert abs(numerical-float(chord_squared(a,b)))<1e-14
    low=covered_chord_lower_squared(F(3,4),F(5,4),F(1,10))
    assert low<=chord_squared(F(9,10),F(11,10))

def test_gain_annulus_and_proper_whitening():
    f=np.array([1+2j,3-1j]);y=np.array([4+1j,-2+3j]);sigma=.2;ref=.9+.1j
    r,g=profile(y,f,sigma,ref)
    objective=2*np.sum(abs(y-g*f)**2)/sigma**2+2*abs(ref-g)**2/.01**2
    assert abs(r@r-objective)<1e-9
    assert .75-1e-14<=abs(g)<=1.25+1e-14
    for a in np.linspace(.75,1.25,11):
        for p in np.linspace(-np.pi,np.pi,33):
            gg=a*np.exp(1j*p)
            obj=2*np.sum(abs(y-gg*f)**2)/sigma**2+2*abs(ref-gg)**2/.01**2
            assert objective<=obj+1e-9
    # c=0 has arbitrary phase but the optimum radius is a_min, not zero.
    r,g=profile(np.zeros_like(f),f,1.)
    assert abs(abs(g)-.75)<1e-14

def test_uniform_cover_and_fail_closed():
    good=CertificateAssumptions('identity interval proof','exact forward model',F(1,100),F(0))
    ans=finite_pair_cover([I(0,1)],lambda x:x,[0],[F(1,10)],good,10000)
    assert ans['status']=='certified'
    fail=finite_pair_cover([I(0,1)],lambda x:(I(0),),[0],[F(1,10)],good,100)
    assert fail['status']=='unresolved' and fail['unprocessed_boxes']>0
    missing=CertificateAssumptions('','','0','0')
    assert finite_pair_cover([I(0,1)],lambda x:x,[0],[F(1,10)],missing)['status']=='unresolved'
    assert class_C_status()['status']=='unresolved'

def test_coercivity_rational_constants():
    c=constants();assert F(c['coercivity_lower'])==F(43,500)
    assert c['actual_continuum_residual_bound'] is None

def test_certified_readout_implementation():
    from modal_readout import amplitude_readout,ratio_readout
    for i,e in [(0,2.3),(1,3.1)]:
        x=[.2,.15][i];m=mie(e,x,1);g=.9+.2j
        bracket=amplitude_readout(g*m[1],g,i)
        assert abs(float((bracket.lo+bracket.hi)/2)-e)<1e-9
        assert bracket.hi-bracket.lo<=F(1,10**10)
        bracket=ratio_readout(g*m[1],g*m[0],i)
        assert abs(float((bracket.lo+bracket.hi)/2)-e)<1e-8

def test_rational_hardware_budget_and_size_allocation():
    from certify_budgets import certificate
    c=certificate()
    assert all(F(v['absolute_material_error_upper'])<F(1,10) for v in c['conditional_reference_readout'])
    assert c['hardware_attainability'].startswith('NOT')

def test_sqrt_outward_and_secant_evidence():
    from certify_secant import sqrt_bounds
    import json
    for x in [F(0),F(1,3),F(2),F(12345,678)]:
        v=sqrt_bounds(I(x));assert v.lo*v.lo<=x<=v.hi*v.hi
    p=Path(__file__).resolve().parents[1]/'results/secant_certificate_final.json'
    c=json.loads(p.read_text())
    assert all(F(r['q_second_derivative_upper'])<0 for r in c['results'])
    assert all(F(r['additional_amplitude_bias_lower'])>0 for r in c['results'])
