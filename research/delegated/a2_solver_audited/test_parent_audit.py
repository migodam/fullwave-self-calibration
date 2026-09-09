"""Parent regression tests: scaling, reference covariance, and charged work."""
import sys
from pathlib import Path
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(Path(__file__).parent),str(ROOT/'research/delegated/a2_physics'),str(ROOT/'research/trispace_self_calibration/a2_research')]
from physics import Config,Model
from common import CostLedger,Calibration,default_settings
from solver import RunState,StageGuard,make_exact_fg,_scaled_from_params

class ParentTests(unittest.TestCase):
    def test_scaled_gradient_both_likelihoods(self):
        m=Model(Config(N=8)); a=np.full(9,.5); x=np.array([.01,-.03,.02])
        y=m.forward(a,np.zeros(3))['total']; sigma=.004
        z=_scaled_from_params(a,x,.5)
        v=np.random.default_rng(54).normal(size=12);v/=np.linalg.norm(v)
        for ph in (False,True):
            r=RunState('phaseless' if ph else 'direct',a,x,CostLedger(),default_settings(),1000000)
            fg=make_exact_fg(m,y,sigma,(0,1),r,StageGuard(r,3),ph,None)
            _,g=fg(z);e=1e-6
            fd=(fg(z+e*v)[0]-fg(z-e*v)[0])/(2*e)
            self.assertLess(abs(fd-g@v)/max(1,abs(fd)),2e-6)
    def test_products_are_not_free(self):
        l=CostLedger(calibration=Calibration());l.charge_model_work({'rhs_solves_total':2,'operator_products':3})
        self.assertGreater(l.units,2)
        u=l.units;l.charge_basis_products(12);self.assertGreater(l.units,u)
    def test_joint_pose_reference(self):
        from reference import reference_metrics
        m=Model(Config(N=8));mr=Model(Config(N=16))
        ref=reference_metrics(m,mr,.002,'full')
        fw=m.forward(np.full(9,.5),np.zeros(3),jacobian=True)
        j=np.sqrt(2)/.002*np.vstack([np.c_[fw['A'],fw['B']].real,np.c_[fw['A'],fw['B']].imag])
        c=np.linalg.inv(j.T@j)[9:,9:];t=np.diag([1,1,1.5])
        self.assertAlmostEqual(ref['tr_C_u_per_3'],np.trace(t@c@t)/3,places=11)
    def test_stack_certificate_quadrature(self):
        from cert import data_error_bound
        from reduced import reduced_evaluate,build_sensing_basis
        m=Model(Config(N=8));a=np.full(9,.5);x=np.zeros(3);l=CostLedger()
        u=build_sensing_basis(m,x,(0,),8,l)['U']
        o=reduced_evaluate(m,a,x,u,(0,),l,True)
        b,detail=data_error_bound(m,a,o,.004)
        self.assertAlmostEqual(b,np.sqrt(2)*np.linalg.norm(detail['per_state_complex']),places=8)

if __name__=='__main__':unittest.main(verbosity=2)
