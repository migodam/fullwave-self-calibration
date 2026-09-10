"""Independent scalar reconstruction; not interval certification."""
from recovery_cycle1 import *
from scipy.special import spherical_jn as jn, spherical_yn as yn
from scipy.optimize import brentq
from scipy.stats import norm


def q(e,x):
    m=np.sqrt(e);z=m*x
    dj=lambda z:jn(1,z,True)+jn(1,z)/z
    dy=lambda z:yn(1,z,True)+yn(1,z)/z
    return (m*jn(1,z)*dj(x)-jn(1,x)*dj(z))/(m*jn(1,z)*dy(x)-yn(1,x)*dj(z))


def main():
    x=np.array([.2,.15]);e0=np.array([2.,3.]);q0=np.array([q(e,s) for e,s in zip(e0,x)])
    e1=np.array([brentq(lambda e:q(e,s)-1.25*b,lo,hi,xtol=1e-13)
                 for s,b,lo,hi in zip(x,q0,[1.5,2],[4,5])])
    t=lambda qs:1j*qs/(1-1j*qs)
    sig=np.array([2.6664081317897358e-6,1.8047986208244827e-6])
    q1=np.array([q(e,s) for e,s in zip(e1,x)])
    distance=float(np.linalg.norm((t(q0)-.8*t(q1))/sig))
    risk=float(norm.cdf(-distance/np.sqrt(2)))
    assert np.all(abs(e1-e0)>.2) and .443<risk<.444
    m=np.array([.000429,.000129]);h=np.array([.00272,.00131])
    bn=np.array([7.008016195904829e-6,4.74348162019559e-6]);ba=.0032905267314919254
    result={'world0':e0.tolist(),'world1':e1.tolist(),'gains':[1.,.8],
        'complex_whitened_distance':distance,'equal_prior_binary_error':risk,
        'conditional_reference_material_bound':((bn+h*ba)/((.75-ba)*m)).tolist(),
        'conditional_remaining_model_error_budget':(.1*(.75-ba)*m-bn-h*ba).tolist(),
        'status':'scalar numerical reconstruction passed; interval and readout feasibility unresolved',
        'scope':'Separately measured lossless electric modes, known geometry; not simultaneous lossy imaging',
        'source_hash':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    target=HERE/'results/parent_boundary.json'
    if target.exists():raise RuntimeError('Preserve immutable result')
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
