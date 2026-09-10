"""Numerical consistency checks. Passing these is not a continuum solver bound."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
from pathlib import Path
import sys,json,unittest,hashlib
import numpy as np
from scipy.special import spherical_jn,spherical_yn
import mpmath as mp
from multipole import basis,mie,SphereCluster,ShiftReadout,illuminations,LOSS
from modal_interval import rat,enc,enclose,bounds,Q,series
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'research/trispace_self_calibration/a3_research'))
import maxwell3d as mx
RESULT={}

def rel(a,b): return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))

class Checks(unittest.TestCase):
    def test_01_boundary_and_passivity(self):
        errors=[];lossless=[];passive=[]
        for eps in (1.5,2.,4.,5.,2+.03j,4+.05j):
            for x in (.15,.2,.45,.63):
                t=mie(eps,x,5); nm=35;offset=0;m=np.sqrt(complex(eps))
                for l in range(1,6):
                    j=spherical_jn(l,x);h=j+1j*spherical_yn(l,x)
                    dj=spherical_jn(l,x,True)+j/x
                    dh=dj+1j*(spherical_yn(l,x,True)+spherical_yn(l,x)/x)
                    v=spherical_jn(l,m*x);dv=spherical_jn(l,m*x,True)+v/(m*x)
                    tm=t[offset];te=t[nm+offset]
                    cm=(j+tm*h)/v;ce=(j+te*h)/(m*v)
                    errors.extend([abs(dj+tm*dh-m*cm*dv)/max(abs(dj),1e-300),abs(dj+te*dh-ce*dv)/max(abs(dj),1e-300)])
                    vals=np.real([tm,te])+np.abs([tm,te])**2
                    if np.imag(eps)==0: lossless.extend(abs(vals))
                    else: passive.extend(vals)
                    offset+=2*l+1
        RESULT['boundary_relative_residual_max']=max(errors)
        RESULT['lossless_energy_residual_max']=float(max(lossless))
        RESULT['lossy_energy_residual_max']=float(max(passive))
        self.assertLess(max(errors),1e-10);self.assertLess(max(lossless),1e-12)
        self.assertLessEqual(max(passive),1e-12)

    def test_02_curl_and_plane_wave(self):
        points=np.array([[.071,.048,.091],[.081,-.09,.065]])
        k=18.;L=2;h=1e-6;B=basis(points,k,L)
        deriv=[]
        for axis in range(3):
            e=np.eye(3)[axis]*h
            deriv.append(((basis(points+e,k,L)-basis(points-e,k,L))/(2*h)).reshape(2,3,-1))
        curl=np.stack([deriv[1][:,2]-deriv[2][:,1],deriv[2][:,0]-deriv[0][:,2],deriv[0][:,1]-deriv[1][:,0]],axis=1).reshape(6,-1)/k
        nm=L*(L+2);target=np.c_[B[:,nm:],B[:,:nm]]
        err=rel(curl,target);RESULT['curl_identity_relative_error']=err;self.assertLess(err,1e-7)
        c=SphereCluster(8,18,36,centers=[[0,0,0]],radii=[.035])
        rx=mx.receivers(17,.023)
        pred=basis(rx,k,8,False)@c.inc
        true=np.stack([np.exp(1j*k*(rx@d))[:,None]*p for d,p in illuminations()],axis=-1).reshape(-1,4)
        err=rel(pred,true);RESULT['plane_wave_L8_relative_error']=err;self.assertLess(err,1e-7)

    def test_03_independent_rayleigh(self):
        radius=.0008;eps=2.3+.03j;k=18.
        model=SphereCluster(3,18,36,centers=[[0,0,0]],radii=[radius])
        rx=mx.receivers(11,.2)
        E=np.stack([p for d,p in illuminations()],axis=-1)
        alpha=4*np.pi*radius**3*(eps-1)/(eps+2)
        ray=(mx.dipole_kernel(rx,np.zeros((1,3)),k)@(alpha*E)).ravel()
        full=model.field([eps],rx).ravel()
        err=rel(full,ray);RESULT['rayleigh_relative_error']=err;self.assertLess(err,.002)

    def test_04_order_quadrature_and_interpolation(self):
        rx=mx.receivers(12,.6);theta=np.array([2.,3.,.725])
        m3=SphereCluster(3);m4=SphereCluster(4);m5=SphereCluster(5);q3=SphereCluster(3,26,52)
        f5=m5.forward(theta,rx);f3=m3.forward(theta,rx)
        RESULT['L3_vs_L5_relative_difference']=rel(f3,f5)
        RESULT['L4_vs_L5_relative_difference']=rel(m4.forward(theta,rx),f5)
        RESULT['quadrature_relative_difference']=rel(q3.forward(theta,rx),f3)
        interp=ShiftReadout(m3,rx)
        errors=[]
        for shift in (-2.,-1.733,-.37,.725,1.917,2.):
            theta[2]=shift;errors.append(rel(interp(theta),m3.forward(theta,rx)))
        RESULT['readout_interpolation_relative_error_max']=max(errors)
        self.assertLess(RESULT['L3_vs_L5_relative_difference'],1e-4)
        self.assertLess(RESULT['quadrature_relative_difference'],1e-8)
        self.assertLess(max(errors),1e-10)

    def test_05_interval_formulas(self):
        mp.iv.dps=60;mp.mp.dps=70
        errors=[];enclosure_pass=0;d2errors=[]
        def qdirect(e,x):
            m=mp.sqrt(e)
            j=lambda z: mp.sin(z)/z**2-mp.cos(z)/z
            y=lambda z: -mp.cos(z)/z**2-mp.sin(z)/z
            D=lambda fun,z:mp.diff(fun,z)+fun(z)/z
            return (m*j(m*x)*D(j,x)-j(x)*D(j,m*x))/(m*j(m*x)*D(y,x)-y(x)*D(j,m*x))
        for xq,eqs in [(Q(1,5),[Q(3,2),Q(2),Q(7,2),Q(4)]),(Q(3,20),[Q(2),Q(3),Q(9,2),Q(5)])]:
            x=mp.mpf(xq.numerator)/xq.denominator
            for eq in eqs:
                e=mp.mpf(eq.numerator)/eq.denominator
                z=enclose(rat(eq),rat(xq))
                direct={'q':qdirect(e,x),'qp':mp.diff(lambda ee:qdirect(ee,x),e)}
                for key,val in direct.items():
                    lo,hi=bounds(z[key]);a=mp.mpf(lo.numerator)/lo.denominator;b=mp.mpf(hi.numerator)/hi.denominator
                    # Independent finite precision reference is tested with 1e-60 slack,
                    # never used to widen or constitute the interval certificate.
                    self.assertTrue(a-mp.mpf('1e-60')<=val<=b+mp.mpf('1e-60'))
                    enclosure_pass+=1;errors.append(float(abs((a+b)/2-val)))
                zz=mp.sqrt(e)*x
                j=lambda t:mp.sin(t)/t**2-mp.cos(t)/t
                D=lambda t:mp.diff(j,t)+j(t)/t
                rhs=-mp.diff(j,zz,2)/zz-mp.diff(j,zz)+2*mp.diff(j,zz)/zz**2-2*j(zz)/zz**3
                d2errors.append(float(abs(mp.diff(D,zz,2)-rhs)))
        RESULT['independent_modal_point_checks']=enclosure_pass
        RESULT['independent_modal_absolute_difference_max']=max(errors)
        RESULT['correct_Dj_second_derivative_difference_max']=max(d2errors)
        self.assertLess(max(errors),1e-40);self.assertLess(max(d2errors),1e-50)

if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(Checks)
    run=unittest.TextTestRunner(verbosity=2).run(suite)
    RESULT.update({'tests':run.testsRun,'failures':len(run.failures),'errors':len(run.errors),'passed':run.wasSuccessful(),
                   'scope':'numerical consistency checks, not certified Maxwell truncation bounds',
                   'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('multipole.py'),Path(__file__).with_name('modal_interval.py')]}})
    out=Path(__file__).with_name('results')/'checks.json';out.write_text(json.dumps(RESULT,indent=2)+'\n')
    print(json.dumps(RESULT,indent=2));sys.exit(0 if run.wasSuccessful() else 1)
