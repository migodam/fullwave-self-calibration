"""Electronics-reference and polarization controls for the restricted EM model."""
from __future__ import annotations
import numpy as np
from scipy.optimize import least_squares
from modal import *
from experiments import complex_noise


def diversity_controls():
    rows=[]
    for orientation,r in [('generic',np.array([.11,.06,.12])),('axial',np.array([0.,0.,.18]))]:
      for nr in [1,2,3]:
       for ni in [1,2,3]:
        for gain in ['scalar','row','entry']:
         Y=dyad(r,10.);D=dyad_jacobian(r,10.).reshape(3,3,3)
         m=Y[:nr,:ni].ravel();B=D[:nr,:ni,:].reshape(-1,3)
         if gain=='scalar':N=m[:,None]
         elif gain=='row':
          N=np.zeros((nr*ni,nr),complex)
          for i in range(nr):N[i*ni:(i+1)*ni,i]=m[i*ni:(i+1)*ni]
         else:N=np.diag(m)
         V,rank=profile_columns(realify(B)/np.linalg.norm(m),realify(np.column_stack([N,1j*N])))
         sv=np.linalg.svd(V,compute_uv=False);sv=np.r_[sv,np.zeros(3-len(sv))]
         rows.append({'orientation':orientation,'receiver_components':nr,'illuminations':ni,'gain_model':gain,
                      'singular_values':sv,'visible_rank':int(np.sum(sv>1e-7)),'nuisance_rank':rank})
    return rows


def reference_controls(seed=94032,n_scenes=16):
    rng=np.random.default_rng(seed);rows=[]
    # Geometry is externally known in this material-only diagnostic. Do not
    # count this as independent success of joint geometry/material calibration.
    r=np.array([.13,.07,.11]);ks=[6.,12.,22.];a=.045
    for i in range(n_scenes):
        eps=rng.uniform(1.5,6.);g=(.8+.4*rng.random(3))*np.exp(1j*rng.uniform(-.5,.5,3))
        patterns=[sphere_field(r,k,eps+.05j,a).ravel() for k in ks]
        sigma=np.array([np.linalg.norm(p)/180 for p in patterns]);rsigma=.015
        y=[g[j]*patterns[j]+sigma[j]*complex_noise(rng,9) for j in range(3)]
        refs=g+rsigma*complex_noise(rng,3)
        # Profiled no-reference loss is exactly epsilon-independent as long as
        # the electric Mie coefficient is nonzero.
        profile_losses=[]
        for e in [1.6,2.7,4.8,6.8]:
            profile_losses.append(sum(np.linalg.norm(profile_complex(y[j],sphere_field(r,k,e+.05j,a))[1]/sigma[j])**2 for j,k in enumerate(ks)))
        def residual(t):
            gg=t[1:4]+1j*t[4:7]
            vals=[(gg[j]*sphere_field(r,k,t[0]+.05j,a).ravel()-y[j])/sigma[j] for j,k in enumerate(ks)]
            return realify(np.r_[np.concatenate(vals),(gg-refs)/rsigma])
        t0=np.r_[2.7,refs.real,refs.imag]
        fit=least_squares(residual,t0,bounds=(np.r_[1.15,np.full(6,-2.)],np.r_[7.,np.full(6,2.)]),max_nfev=150,x_scale='jac')
        J=fit.jac;V,_=profile_columns(J[:,:1],J[:,1:]);data_rows=np.r_[np.arange(27),np.arange(30,57)]
        before,_=profile_columns(J[data_rows,:1],J[data_rows,1:])
        rows.append({'scene':i,'true_epsilon':eps,'estimated_epsilon':fit.x[0],
                     'material_relative_error':float(abs(fit.x[0]-eps)/eps),'no_reference_profile_loss_range':float(np.ptp(profile_losses)),
                     'material_visible_norm_with_reference':float(np.linalg.norm(V)),
                     'material_visible_norm_without_reference':float(np.linalg.norm(before)),
                     'action_before':'add_electronics_reference','reference_complex_cost':3,
                     'reference_noise_sigma':rsigma,'no_reference_status':'material_unidentifiable_reject',
                     'scope':'known geometry; free per-frequency complex gains; one real material parameter',
                     'success':bool(fit.success)})
    return rows
