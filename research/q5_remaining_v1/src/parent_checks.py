"""Parent audit: outgoing Mie convention and profiled complex gain."""
from recovery_cycle1 import *
from scipy.special import spherical_jn as jn, spherical_yn as yn


def main():
    errors=[]
    for eps in (1.5,2.,3.,4.,5.):
        for x in (.15,.2):
            m=np.sqrt(eps);z=m*x
            psi=lambda z:z*jn(1,z)
            dp=lambda z:jn(1,z)+z*jn(1,z,True)
            cy=lambda z:z*yn(1,z)
            dy=lambda z:yn(1,z)+z*yn(1,z,True)
            n=m*psi(z)*dp(x)-psi(x)*dp(z)
            d=m*psi(z)*dy(x)-cy(x)*dp(z)
            q=n/d
            a=n/(n+1j*d)
            assert abs(1j*a/(1-a)-q)<1e-15
            t=1j*q/(1-1j*q)
            eigen=np.linalg.eigvals(np.asarray(treams.TMatrix.sphere(1,1.,x,[treams.Material(eps),treams.Material()])))
            errors.append(float(min(abs(eigen-t))))
    assert max(errors)<1e-14
    rng=np.random.default_rng(12)
    f=rng.normal(size=20)+1j*rng.normal(size=20)
    y=(.9+.2j)*f+.01*(rng.normal(size=20)+1j*rng.normal(size=20))
    r,g=profile(y,f,.02,.9+.21j)
    direct=np.sum(abs(y-g*f)**2)/.02**2+abs((.9+.21j)-g)**2/.01**2
    assert np.isclose(r@r,direct)
    # Whitened proper complex Gaussian real coordinates have variance one.
    w=(rng.normal(size=100000)+1j*rng.normal(size=100000))/np.sqrt(2)
    assert abs(np.mean(2*abs(w)**2)-2)<.03
    target=HERE/'results/parent_checks.json'
    if target.exists():raise RuntimeError('Preserve immutable check')
    result={'pass':True,'mie_points':10,'max_treams_eigenmode_error':max(errors),
        'mie_convention':'a=N/(N+iD), T=-a, q=i*a/(1-a)=-i*T/(1+T)',
        'profile_objective_check':True,'proper_complex_whitening_check':True,
        'scope':'Unit checks only; no continuum or finite-domain material certificate',
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
