"""Regression and independent numerical checks; not substitutes for proofs."""
import json,sys,unittest
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import eval_legendre,spherical_jn,spherical_yn
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from dyadic_interval import I,SCALE,modal,spherical1
from certify_geometry_window import orientation_polynomial
from certify_readout import upper_tail
from sphere_cluster import SphereCluster,mie_diagonal,sphere_quadrature,vector_waves
from recovery_v2 import profiled_gain,realwhite

class Checks(unittest.TestCase):
    def contains(self,iv,v):
        v=F(v)
        self.assertLessEqual(F(iv.lo,SCALE),v)
        self.assertGreaterEqual(F(iv.hi,SCALE),v)

    def test_exact_rational_and_rejection(self):
        for v in [F(1,5),F(-7,13),F(0),F(3,20)]:
            self.contains(I.rational(v.numerator,v.denominator),v)
        with self.assertRaises(TypeError): I.rational(.2)
        with self.assertRaises(TypeError): I.rational(1)+.2

    def test_integer_operations(self):
        for a,b in [(F(-2,3),F(-7,9)),(F(1,17),F(4,3)),(F(-1,2),F(7,3))]:
            x=I.rational(a.numerator,a.denominator);y=I.rational(b.numerator,b.denominator)
            for iv,val in [(x+y,a+b),(x-y,a-b),(x*y,a*b),(x/y,a/b),(x**4,a**4)]:self.contains(iv,val)
        self.assertEqual((I.box(-1,1,2,1)**2).lo,0)
        with self.assertRaises(ZeroDivisionError): I.box(-1,1,1,1).reciprocal()

    def test_square_root_integer_enclosure(self):
        x=I.rational(2).sqrt()
        self.assertLessEqual(x.lo*x.lo,2*SCALE*SCALE)
        self.assertGreaterEqual(x.hi*x.hi,2*SCALE*SCALE)

    def test_decimal_export_outward(self):
        for x in [I.rational(-1,7),I.rational(1,9),I(-3,7)]:
            lo,hi=map(lambda s:F(Decimal(s)),x.decimal_bounds())
            self.assertLessEqual(lo,F(x.lo,SCALE));self.assertGreaterEqual(hi,F(x.hi,SCALE))

    def test_trig_against_high_precision(self):
        mp.mp.dps=80
        for p in [-9,-1,2,7,10]:
            x=I.rational(p,10)
            for iv,v in [(x.sin(),mp.sin(mp.mpf(p)/10)),(x.cos(),mp.cos(mp.mpf(p)/10))]:
                self.assertLessEqual(mp.mpf(iv.lo)/SCALE,v)
                self.assertGreaterEqual(mp.mpf(iv.hi)/SCALE,v)

    def test_correct_second_derivative(self):
        mp.mp.dps=80
        def j(z):return mp.sin(z)/z**2-mp.cos(z)/z
        def dj(z):return mp.diff(j,z)+j(z)/z
        for p in [15,20,45,90]:
            iv=spherical1(I.rational(p,100))[5];v=mp.diff(dj,mp.mpf(p)/100,2)
            self.assertLessEqual(mp.mpf(iv.lo)/SCALE,v)
            self.assertGreaterEqual(mp.mpf(iv.hi)/SCALE,v)

    def test_mie_sign_and_reaction_identity(self):
        for e,x in [(2.,.2),(4.7,.15),(3.4,.2)]:
            a=-mie_diagonal(e,x,2)[1];q=1j*a/(1-a)
            iv=modal(I.rational(round(e*10),10),I.rational(round(x*100),100))['q']
            self.assertLess(abs(q.imag),1e-14)
            self.assertLess(abs(q.real-iv.diagnostic_midpoint()),1e-15)

    def test_gain_annulus_and_gls_objective(self):
        rng=np.random.default_rng(21);f=rng.normal(size=7)+1j*rng.normal(size=7)
        y=.3*f+.1*(rng.normal(size=7)+1j*rng.normal(size=7));z=.4+.2j;s=.13;sr=.01
        g=profiled_gain(f,y,s,z,sr)
        self.assertGreaterEqual(abs(g),.75-1e-14);self.assertLessEqual(abs(g),1.25+1e-14)
        def loss(c):return np.linalg.norm(y-c*f)**2/s**2+abs(z-c)**2/sr**2
        phase=np.angle(np.vdot(f,y)/s**2+z/sr**2)
        opt=minimize_scalar(lambda a:loss(a*np.exp(1j*phase)),bounds=(.75,1.25),method='bounded')
        self.assertLessEqual(loss(g),opt.fun+1e-9)
        self.assertAlmostEqual(np.dot(realwhite(y),realwhite(y)),2*np.vdot(y,y).real,places=11)

    def test_finite_chord_and_orientation_polynomial(self):
        def a(z):return z*z+1j*z-1
        def c(z):return 2-2j*z
        for zp,wp in [(6,8),(12,14),(30,40)]:
            z=zp/10;w=wp/10
            v=np.array([np.sqrt(2)*a(z),c(z)]);u=np.array([np.sqrt(2)*a(w),c(w)])
            chord=1-abs(np.vdot(v,u))**2/(np.vdot(v,v).real*np.vdot(u,u).real)
            formula=2*(z-w)**2*((z+w)**2+z*z*w*w)/((z**4+z*z+3)*(w**4+w*w+3))
            self.assertAlmostEqual(chord,formula,places=13)
            s1=2*np.conj(a(z))*a(w)+np.conj(c(z))*c(w)
            s0=np.conj(a(z))*a(w)+np.conj(a(z))*c(w)+np.conj(c(z))*a(w)
            p=orientation_polynomial(I.rational(zp,10),I.rational(wp,10))
            self.assertAlmostEqual((abs(s1)**2-abs(s0)**2)/p.diagnostic_midpoint(),1,places=12)

    def test_coverage_metadata(self):
        d=json.loads((ROOT/'results/modal_certificate.json').read_text())
        self.assertEqual(sum(r['coverage']['count'] for r in d['rows']),11000)
        for r in d['rows']:
            c=r['coverage'];self.assertEqual(c['count'],c['exclusive_last_index']-c['first_index'])
            self.assertFalse(r['nonpositive_boxes'])
        g=json.loads((ROOT/'results/geometry_window_certificate.json').read_text())
        bins=[r['nominal_z_bin'] for r in g['rows'] if r['finite_radius_error_certified']]
        self.assertEqual([[F(*a),F(*b)] for a,b in bins],[[F(i,10),F(i+1,10)] for i in range(6,16)])
        self.assertTrue(all(r['orientation_unproved_boxes']==0 for r in g['rows']))

    def test_radial_quadrature_orthogonality(self):
        n,w=sphere_quadrature(8,4);u=n[:,2];psi=n[:,0]
        den=np.sum(w*psi**2)
        self.assertAlmostEqual(den,4*np.pi/3,places=13)
        for l in range(2,15):
            dp=l*(eval_legendre(l-1,u)-u*eval_legendre(l,u))/(1-u*u)
            self.assertLess(abs(np.sum(w*psi*n[:,0]*dp)/den),2e-13)

    def test_readout_against_vector_spherical_waves(self):
        model=SphereCluster(centers=np.zeros((1,3)),radii=[.2/6],k=6,order=7,ntheta=12)
        n,w=sphere_quadrature(8,4);b=mie_diagonal(3.,.2,7)*model.incident[:,0]
        e=(vector_waves(n/6,6,7,True)@b).reshape(-1,3);er=np.sum(e*n,axis=1);psi=n[:,0]
        h=spherical_jn(1,1)+1j*spherical_yn(1,1)
        estimate=np.sum(w*psi*er)/(-3*h*np.sum(w*psi**2))
        self.assertLess(abs(estimate+mie_diagonal(3.,.2,7)[1]),1e-14)

    def test_tail_majorant_monotonic(self):
        for x in [F(3,20),F(1,5),F(1,4)]:
            self.assertLess(upper_tail(x,15),F(1,10**24))
            self.assertLess(upper_tail(x,3),upper_tail(x,2))

    def test_readout_cross_terms_and_budget(self):
        d=json.loads((ROOT/'results/readout_certificate.json').read_text())
        for r in d['budget_rows']:
            q=r['quantities'];bounds=lambda k:list(map(Decimal,q[k]['outward_decimal']))
            self.assertGreater(bounds('node_calibration')[0],0)
            self.assertEqual(r['pass_delta_0p1'],r['repetitions']==16 or r['sphere']==1)
        source=(ROOT/'certify_readout.py').read_text()
        self.assertIn('F(1,1000)*F(1001,1000)',source)
        self.assertIn('F(5,4)*F(1001,1000)*tail',source)

    def test_finite_modal_competitors(self):
        ta=-(-mie_diagonal(4.70,.15,1)[1]);tb=-(-mie_diagonal(4.91,.15,1)[1])
        ka=2*tb/(ta+tb);kb=2*ta/(ta+tb)
        self.assertLess(abs(ka*ta-kb*tb),1e-17)
        self.assertLess(max(abs(ka-1),abs(kb-1)),.0125)
        self.assertGreater(4.91-4.70,.2)
        self.assertTrue(.75<abs(ta/tb)<1.25)

    def test_shared_nuisance_fusion_identity(self):
        rng=np.random.default_rng(31);h=np.array([.2,-.1]);blocks=[]
        for _ in range(3):
            a=rng.normal(size=(8,2));n=rng.normal(size=(8,3));hh=n.T@n
            v=-np.linalg.solve(hh,n.T@a@h);q=np.linalg.norm(a@h+n@v)**2
            blocks.append((a,n,hh,v,q))
        H=sum(b[2] for b in blocks);vm=np.linalg.solve(H,sum(b[2]@b[3] for b in blocks))
        joint=sum(np.linalg.norm(a@h+n@vm)**2 for a,n,_,_,_ in blocks)
        rhs=sum(q+(v-vm)@hh@(v-vm) for _,_,hh,v,q in blocks)
        self.assertAlmostEqual(joint,rhs,places=12)

    def test_geometry_to_readout_budget(self):
        d=json.loads((ROOT/'results/geometry_readout_budget.json').read_text())
        for e in d['endpoint_transfer_errors']:
            self.assertGreater(Decimal(e['relative_electric_dipole_transfer_error']['outward_decimal'][0]),Decimal('.04'))
        self.assertLess(Decimal(d['fine_halfwidth_dipole_transfer_error_bound']['outward_decimal'][1]),Decimal('.001'))

    def test_raw_preservation_and_scene_count(self):
        import hashlib
        path=ROOT/'results/recovery_v2';s=json.loads((path/'summary.json').read_text())
        self.assertEqual(hashlib.sha256((path/'raw_observations.json').read_bytes()).hexdigest(),s['raw_sha256'])
        self.assertEqual(len(list(path.glob('scene_*.json'))),12)
        self.assertEqual([s['summary'][k]['successes'] for k in ['no_reference','noisy_reference_GLS','fixed_EM','random_EM']],[2,7,2,2])

if __name__=='__main__':unittest.main(verbosity=2)
