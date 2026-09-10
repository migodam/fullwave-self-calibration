"""Regression checks complement, but do not replace, continuum proofs."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,math
import numpy as np
import pytest
from scipy.special import spherical_jn,spherical_yn
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from interval_certificate import I,SCALE,coefficient,modal_intervals,series
from multipole import Cluster,angular,quadrature,sphere_t,vector_waves,boundary_residual,modes
from readout_certificate import tail_term,tail_bound,odd_df
from modal_estimator import readout_operator,project_modal_field,material_from_amplitude,require_continuum_certificate
from run_experiment import profile,HERE


def encloses(iv,x):return F(iv.lo,SCALE)<=x<=F(iv.hi,SCALE)

@pytest.mark.parametrize('a,b',[(F(1,3),F(2,7)),(-F(1,7),F(3,11)),(-F(4,9),-F(2,3))])
def test_integer_outward_arithmetic(a,b):
    x,y=I.of(a),I.of(b)
    assert encloses(x+y,a+b) and encloses(x*y,a*b) and encloses(x/y,a/b)


def test_interval_input_rejections():
    with pytest.raises(ZeroDivisionError): I.of(-1,1).reciprocal()
    with pytest.raises(ValueError): series('J',I.of(-1,0))
    with pytest.raises(ValueError): I(1,0)

@pytest.mark.parametrize('kind',['J','D','C','S'])
def test_analytic_tail_ratio_regression(kind):
    for r in range(3):
        for n in range(12,80):
            a=abs(coefficient(kind,n))*math.prod(range(n-r+1,n+1))
            b=abs(coefficient(kind,n+1))*math.prod(range(n-r+2,n+2))
            assert b/a<F(1,2)


def test_complete_rational_cover_and_proof_source():
    cert=json.loads((HERE/'results/modal_certificate.json').read_text())
    assert cert['source_sha256']==hashlib.sha256((HERE/'interval_certificate.py').read_bytes()).hexdigest()
    assert len(cert['boxes'])==5500
    for obj,lo,hi,count in [(1,F(3,2),F(4),2500),(2,F(2),F(5),3000)]:
        rows=[r for r in cert['boxes'] if r['object']==obj]
        assert len(rows)==count
        previous=lo
        for r in rows:
            a,b=map(F,r['eps']); assert a==previous and b-a==F(1,1000); previous=b
            assert F(r['q'][0])>0 and F(r['q_prime'][0])>0 and F(r['q_second'][1])<0
        assert previous==hi

@pytest.mark.parametrize('eps,x',[(1.5,F(1,5)),(4,F(1,5)),(2,F(3,20)),(5,F(3,20))])
def test_mie_convention_and_interval_regression(eps,x):
    a=-sphere_t(eps,float(x),1)[4];q=1j*a/(1-a)
    iv=modal_intervals(I.of(F(eps)),x)[0]
    assert abs(q.imag)<2e-17
    assert float(F(iv.lo,SCALE))-2e-16<=q.real<=float(F(iv.hi,SCALE))+2e-16


def test_parent_Dj_second_derivative():
    for z in np.linspace(.1,.5,9):
        j=spherical_jn(1,z);jp=spherical_jn(1,z,True)
        jpp=-2*jp/z-(1-2/z**2)*j
        jppp=-(1-4/z**2)*jp-2*jpp/z-4*j/z**3
        direct=jppp+jpp/z-2*jp/z**2+2*j/z**3
        corrected=-jpp/z-jp+2*jp/z**2-2*j/z**3
        assert abs(direct-corrected)<1e-10


def test_spherical_polynomial_quadrature():
    d,w=quadrature(3,6)
    for a in range(6):
        for b in range(6-a):
            for c in range(6-a-b):
                got=np.sum(w*d[:,0]**a*d[:,1]**b*d[:,2]**c)
                if a%2 or b%2 or c%2: exact=0.
                else: exact=2*math.gamma((a+1)/2)*math.gamma((b+1)/2)*math.gamma((c+1)/2)/math.gamma((a+b+c+3)/2)
                assert abs(got-exact)<2e-14


def test_readout_noise_row_and_low_mode_orthogonality():
    d,row=readout_operator(); assert abs(np.sum(abs(row)**2)-224/1053)<2e-15
    v=vector_waves(2*d,[0,0,0],1,4,True)
    val=np.einsum('pc,pcj->j',row,v)
    target=len(modes(4))+1
    val[target]=0
    assert np.max(abs(val))<1e-12


def test_all_order_bounds_regression_not_proof():
    for x in (.15,.2):
        for e in (1.5,3.,5.):
            t=sphere_t(e,x,12);nm=len(modes(12))
            for ell in range(1,13):
                idx=modes(12).index((ell,0)); scale=x**(2*ell+1)/(odd_df(2*ell-1)*odd_df(2*ell+1))
                assert abs(t[idx])<=4*scale and abs(t[nm+idx])<=10*scale
    for ell in range(2,30):
        assert tail_term(ell+1,F(1,5))/tail_term(ell,F(1,5))<F(1,100)


def test_plane_wave_reconstruction_and_boundary():
    wave=(np.array([1.,0,0]),np.array([0.,0,1.]))
    model=Cluster([[0,0,0]],[.2],1,[wave],order=8)
    dirs,_=quadrature(4,7);pts=.1*dirs
    reconstructed=vector_waves(pts,[0,0,0],1,8)@model.incident[:,0]
    truth=np.exp(1j*pts[:,0])[:,None]*wave[1]
    assert np.linalg.norm(reconstructed-truth)/np.linalg.norm(truth)<1e-12
    assert boundary_residual(3+.05j,.63,5)<1e-13


def test_modal_recovery_from_independent_field_evaluation():
    d,_=readout_operator();illum=[(np.array([1.,0,0]),np.array([0.,0,1.]))]
    for i,(e,x) in enumerate([(2.7,.2),(3.6,.15)]):
        field=Cluster([[0,0,0]],[x],1,illum,order=8).field([e],2*d)[:,:,0]
        g=.82*np.exp(.4j);y=project_modal_field(g*field)
        assert abs(material_from_amplitude(y,abs(g),i)-e)<1e-8
    with pytest.raises(ValueError):material_from_amplitude(1,-1,0)


def test_gain_profile_annulus_and_real_whitening():
    f=np.array([1+.3j,-.2+.4j]);sigma=.1
    for gain in (.1+0j,1+.2j,2j):
        y=gain*f;ref=(.9+.1j,.01)
        r,gh=profile(y,f,sigma,ref)
        assert .75-1e-15<=abs(gh)<=1.25+1e-15
        def obj(g):return 2*(np.linalg.norm((y-g*f)/sigma)**2+abs((ref[0]-g)/ref[1])**2)
        assert abs(r@r-obj(gh))<1e-8
        phase=np.linspace(-np.pi,np.pi,101)
        assert obj(gh)<=min(obj(rad*np.exp(1j*p)) for rad in (.75,1.,1.25) for p in phase)+1e-8


def test_frozen_counts_raw_data_and_source_hashes():
    d=json.loads((HERE/'results/independent_v2.json').read_text());assert d['complete'] and len(d['rows'])==12
    counts={key:0 for key in ['no_reference','noisy_reference_GLS','fixed_EM','random_EM']}
    for row in d['rows']:
        assert np.array(row['observed_all_fields']).shape==(336,2)
        for key,m in row['methods'].items():
            assert len(m['starts'])==3
            assert m['task_success']==all(x<=.1 for x in m['material_errors'])
            counts[key]+=m['task_success']
    assert list(counts.values())==[2,7,2,2]
    root=HERE.parents[1]
    for p,h in d['source_hashes'].items():assert hashlib.sha256((root/p).read_bytes()).hexdigest()==h


def test_finite_competitor_domain_and_scope():
    d=json.loads((HERE/'results/finite_pair.json').read_text())
    for pair in d['pairs']:
        a,b=pair['world_A'],pair['world_B'];gb=complex(*pair['gain_B'])
        assert abs(a[1]-b[1])>.2 and 1.5<=b[0]<=4 and -2<=b[2]<=2 and .75<=abs(gb)<=1.25
        assert 'NOT certified' in pair['status'] and pair['complex_sigma_shared']>0


def test_fail_closed_missing_model_error():
    with pytest.raises(ValueError,match='UNRESOLVED'):require_continuum_certificate({'optimizer_converged':True})


def test_conditional_budget_has_positive_margin():
    d=json.loads((HERE/'results/readout_certificate.json').read_text())
    assert d['source_sha256']==hashlib.sha256((HERE/'readout_certificate.py').read_bytes()).hexdigest()
    assert F(d['uniform_failure_probability_upper'])<F(1,1000)
    for row in d['rows']:
        assert F(row['uniform_material_error_bound'][1])<F(1,10)
        assert F(row['remaining_margin'][0])>0


def test_finite_window_cover_and_noise_event():
    from finite_windows import radial_squared
    d=json.loads((HERE/'results/finite_windows.json').read_text())
    assert F(d['reference_real_gaussian_tail_upper'])<F(1,1000)
    for row in d['rows']:
        lo,hi=map(F,row['radial_endpoint_exact_enclosure'])
        blo,bhi=map(F,row['max_total_normalized_amplitude_error'])
        assert hi-lo<F(1,10**12)
        assert radial_squared(lo)>(F(1,10**6)/blo)**2
        assert radial_squared(hi)<(F(1,10**6)/bhi)**2
