"""Audit finite claims and preserved raw evidence; no new recovery campaign."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,unittest
import numpy as np
from finite_window import feasible_intervals,require_class_c_certificate
from run_v2 import profile
HERE=Path(__file__).resolve().parent
METRICS={}

class Audit(unittest.TestCase):
    def test_exact_born_pair_and_finite_geometry(self):
        ell=np.array([.03,.05]);g=(30+1j)/(40+1j)
        np.testing.assert_allclose((30+1j)*ell,g*(40+1j)*ell,rtol=2e-15,atol=1e-15)
        self.assertTrue(.75<=abs(g)<=1.25)
        self.assertTrue(np.all(np.abs((1+30*ell)-(1+40*ell))>.2))
        beta=.002;k=18.;ds=2.;X=.002
        phase_bound=2*np.arcsin(beta)
        self.assertLess(k*ds*X+phase_bound,np.pi)
        error_mm=1000*phase_bound/(k*ds)
        self.assertLess(error_mm,.12)
        METRICS['restricted_Born_geometry_error_mm']=error_mm
        METRICS['Born_partner_gain_modulus']=abs(g)

    def test_annulus_profile(self):
        rng=np.random.default_rng(41)
        for gain in [.2*np.exp(.2j),1.0*np.exp(.8j),2.0*np.exp(-.3j)]:
            f=rng.normal(size=8)+1j*rng.normal(size=8)
            y=gain*f+.05*(rng.normal(size=8)+1j*rng.normal(size=8))
            z=gain+.02j
            r,g=profile(y,f,.1,z,.02)
            v=np.r_[f/.1,1/.02];b=np.r_[y/.1,z/.02]
            uncon=np.vdot(v,b)/np.vdot(v,v).real
            expected=np.clip(abs(uncon),.75,1.25)*np.exp(1j*np.angle(uncon))
            self.assertLess(abs(g-expected),1e-13)
            # Test many legal alternatives, not a new optimizer comparison.
            gs=rng.uniform(.75,1.25,2000)*np.exp(1j*rng.uniform(-np.pi,np.pi,2000))
            losses=2*np.sum(abs(b[:,None]-v[:,None]*gs)**2,axis=0)
            self.assertLessEqual(float(r@r),float(losses.min())+1e-10)

    def test_full_rank_not_finite_injectivity(self):
        self.assertEqual((-1.)**2,1.**2)
        self.assertNotEqual(2*(-1.),0)
        self.assertNotEqual(2*(1.),0)
        self.assertEqual(require_class_c_certificate(None,None,1)['status'],'unresolved')
        self.assertEqual(require_class_c_certificate(1,1,1)['status'],'certificate_not_established')

    def test_modal_and_m2_coverage(self):
        cert=json.loads((HERE/'results/modal_certificate.json').read_text())
        path=HERE/'results/modal_cover.jsonl'
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),cert['trace_sha256'])
        rows=[json.loads(s) for s in path.read_text().splitlines()]
        self.assertEqual(len(rows),5500)
        for region in cert['regions']:
            group=[r for r in rows if r['x']==region['size']]
            last=Q(region['material'][0])
            for row in group:
                self.assertEqual(Q(row['lo']),last);last=Q(row['hi'])
                self.assertGreater(Q(row['bounds']['hp'][0]),0)
            self.assertEqual(last,Q(region['material'][1]))
        cert=json.loads((HERE/'results/m2/finite_certificate.json').read_text())
        path=HERE/'results/m2/finite_cover.jsonl'
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),cert['trace_sha256'])
        rows=[json.loads(s) for s in path.read_text().splitlines()]
        self.assertEqual(len(rows),5650)
        for i,(lo,hi) in enumerate([(Q(3,2),Q(4)),(Q(2),Q(5))]):
            for kind,stop in [('finite_increment_derivative',hi-Q(1,5)),('complex_strip_and_size',hi)]:
                last=lo
                for row in [r for r in rows if r['kind']==kind and r['channel']==i]:
                    self.assertEqual(Q(row['box'][0]),last);last=Q(row['box'][1])
                    if kind=='finite_increment_derivative':self.assertLess(Q(row['upper']),0)
                self.assertEqual(last,stop)
        self.assertTrue(all(r['scenarios']['16']['physical_budget_pass'] for r in cert['regions']))
        self.assertFalse(cert['regions'][1]['scenarios']['1']['physical_budget_pass'])

    def test_recovery_raw_integrity(self):
        report=json.loads((HERE/'results/v2/recovery.json').read_text())
        raw=HERE/'results/v2/observations.jsonl'
        self.assertEqual(hashlib.sha256(raw.read_bytes()).hexdigest(),report['observations_sha256'])
        self.assertEqual(len(raw.read_text().splitlines()),12)
        self.assertTrue(report['complete'])
        for name in report['primary_methods']+report['secondary_ablations']:
            count=sum(np.all(np.array(r['methods'][name]['material_errors'])<=.1) for r in report['rows'])
            self.assertEqual(count,report['summary'][name]['successes'])
        for name,sha in report['source_sha256'].items():
            self.assertEqual(hashlib.sha256((HERE.parents[1]/name).read_bytes()).hexdigest(),sha)

if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    (HERE/'results/audit_checks.json').write_text(json.dumps({'tests':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'metrics':METRICS},indent=2)+'\n')
    raise SystemExit(not result.wasSuccessful())
