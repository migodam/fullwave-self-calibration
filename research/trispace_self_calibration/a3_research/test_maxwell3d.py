"""Independent derivative and Maxwell symmetry checks before final experiments."""
import json
from pathlib import Path
import numpy as np
from calibrate3d import JointModel
from maxwell3d import treams_field,receivers,dipole_kernel

if __name__=='__main__':
    rows=[];rng=np.random.default_rng(5102)
    z=np.zeros(14);z[:2]=[2.4,3.];z[2:5]=[.03,-.02,.01];z[5]=.04
    z[6:10]=[.03,-.02,.01,.02];z[10:]=[.1,-.1,.05,.03]
    model=JointModel(.06)
    y,j=model.field_jac(z)
    for label,v in [('material',np.r_[rng.normal(size=2),np.zeros(12)]),
                    ('geometry',np.r_[np.zeros(2),rng.normal(size=3),np.zeros(9)]),
                    ('electronics',np.r_[np.zeros(5),rng.normal(size=9)])]:
        v/=np.linalg.norm(v);h=2e-5
        yp,_=model.field_jac(z+h*v);ym,_=model.field_jac(z-h*v)
        fd=(yp-ym)/(2*h);an=j@v
        error=np.linalg.norm(fd-an)/np.linalg.norm(fd)
        rows.append(dict(name=label+'_directional_derivative',error=float(error),tolerance=2e-7,passed=bool(error<2e-7)))
    # Separate adjoint identity for real parameter Jacobian.
    v=rng.normal(size=14);w=rng.normal(size=y.shape)+1j*rng.normal(size=y.shape)
    left=np.real(np.vdot(w,j@v));right=v@np.real(j.reshape(-1,14).conj().T@w.ravel())
    error=abs(left-right)/max(abs(left),1e-30)
    rows.append(dict(name='real_parameter_adjoint_identity',error=float(error),tolerance=1e-12,passed=bool(error<1e-12)))
    # Full vector dyadic reciprocity, ordinary transpose (not Hermitian).
    a=rng.normal(size=(4,3));b=rng.normal(size=(5,3))
    error=np.linalg.norm(dipole_kernel(a,b,6)-dipole_kernel(b,a,6).T)
    rows.append(dict(name='dyadic_reciprocity',error=float(error),tolerance=1e-12,passed=bool(error<1e-12)))
    ref,xs=treams_field([[0,0,0]],[.12],[2.4+.03j],9,receivers(),5)
    ok=all(0<=sca<=ext for sca,ext in xs)
    rows.append(dict(name='passive_reference_extinction_exceeds_scattering',passed=ok))
    report=dict(scope='development derivative/reference consistency, not inverse gate',checks=rows,passed=all(r['passed'] for r in rows))
    (Path(__file__).resolve().parent/'results/maxwell3d_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    assert report['passed']
