import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pytest
from modal import *
from coverage import Block,Cell,ModalProblem,cover
from dda import grid_target,solve
from risk import weighted_task_risk, Evidence,action

@pytest.mark.parametrize('z',[.02,.1,.5,1.,3.,20.,100.])
def test_singular_formula(z):
    n=np.array([.3,-.4,.5]); n/=np.linalg.norm(n); k=10.;r=z/k*n
    assert np.allclose(visible_singular_values(r,k),theoretical_singular_values(z/k,k),rtol=3e-8,atol=1e-9)

@pytest.mark.parametrize('k',[3.,10.,40.])
def test_derivative(k):
    r=np.array([.12,.07,.11]);h=1e-6
    fd=np.column_stack([(dyad(r+h*v,k)-dyad(r-h*v,k)).ravel()/(2*h) for v in np.eye(3)])
    assert np.linalg.norm(fd-dyad_jacobian(r,k))/np.linalg.norm(fd)<2e-8


def test_sphere_proposal_and_parity():
    r=np.array([.12,.07,.11]); k=15.; y=(2-3j)*dyad(r,k)
    assert np.allclose(y,(2-3j)*dyad(-r,k))
    p=tensor_proposal(y,k)
    assert min(np.linalg.norm(a-r) for a in p)<1e-12


def test_radial_optimum():
    zs=np.logspace(-3,3,6001);sr=np.array([theoretical_singular_values(z/10,10)[2] for z in zs])
    assert abs(zs[sr.argmax()]-1)<1e-10
    assert np.isclose(sr.max(),20/np.sqrt(5))


def test_electronic_degeneracy():
    r=np.array([0.,0.,.1]); k=10
    assert visible_singular_values(r,k,gain_model='row')[-1]<1e-10
    r=np.array([.07,.04,.08])
    assert max(visible_singular_values(r,k,gain_model='entry'))<1e-10

@pytest.mark.parametrize('model',['static','radiative'])
def test_asymptotic_radial_rank_loss(model):
    assert visible_singular_values(np.array([.07,.04,.08]),10,fidelity=model)[-1]<1e-7


def test_dda_derivatives():
    xyz,v=grid_target(h=.025);k=10.;eps=2.5+.05j;r=np.array([.12,.07,.11]);s=solve(xyz,v,k,eps)
    h=1e-6
    fd=np.column_stack([(s.field(r+h*a)-s.field(r-h*a)).ravel()/(2*h) for a in np.eye(3)])
    assert np.linalg.norm(fd-s.geometry_jacobian(r))/np.linalg.norm(fd)<1e-7
    de=(solve(xyz,v,k,eps+h).field(r)-solve(xyz,v,k,eps-h).field(r)).ravel()/(2*h)
    assert np.linalg.norm(de-s.material_jacobian(r))/np.linalg.norm(de)<1e-7
    assert s.relative_linear_residual<1e-12


def test_lipschitz_lower_bound():
    rng=np.random.default_rng(222);r=np.array([.15,.07,.13]); y=100*unit_pattern(r,10)+rng.normal(size=9)+1j*rng.normal(size=9)
    p=ModalProblem([Block(y,10,np.zeros(3))]);c=Cell(r+np.array([.01,-.01,0]),np.array([.02]*3));lb=p.lower_bound(c)
    for q in c.center+rng.uniform(-1,1,(300,3))*c.half:
        assert p.norm(q)+1e-9>=lb


def test_budget_never_forgets_other_branch():
    r=np.array([.15,.07,.13]);y=1000*unit_pattern(r,10)
    p=ModalProblem([Block(y,10,np.zeros(3))]);roots=[Cell(r,np.ones(3)*.02),Cell(-r,np.ones(3)*.02)]
    out=cover(p,roots,r,max_cells=2)
    assert not out['accepted'] and out['unresolved_cells']>0
    assert any(c.contains(-r) for c in out['cells'])


def test_sandwich_not_inverse_hessian():
    rng=np.random.default_rng(3);J=rng.normal(size=(20,4));W=np.diag(np.linspace(.1,1,20));T=np.eye(4)[:2]
    a=weighted_task_risk(J,W,T);L=a['L']
    assert np.isclose(a['variance'],np.trace(T@L@L.T@T.T))
    assert not np.isclose(a['variance'],np.trace(T@np.linalg.inv(J.T@W@J)@T.T))


def test_gate_actions():
    e=Evidence(.1,1.,True,True);assert action(e)=='accept'
    e.electronics_ambiguity=True;assert action(e)=='reject'
    e.reference_available=True;assert action(e)=='add_electronics_reference'
    e.electronics_ambiguity=False;e.coverage_ok=False;e.acquisition_available=True;assert action(e)=='acquire_new_data'
    e.coverage_ok=True;e.fidelity_supported=False;e.refinement_available=True;assert action(e)=='refine_model'
    e.fidelity_supported=True;e.geometry_risk=2.;e.downweight_risk=.5;assert action(e)=='downweight'


def test_approximation_error_mean_sign():
    from experiments import calibration_forward
    theta=np.array([.13,.07,.11,3.,1.,0.,0.])
    yc=calibration_forward(theta,[4.,8.,22.],.045,'coarse')
    yf=calibration_forward(theta,[4.,8.,22.],.045,'fine')
    error=yf-yc
    assert np.linalg.norm(yc-yf+error)<1e-14
    assert np.linalg.norm(yc-yf-error)>1e-5
