"""A2 bounded discrete dissipative certificates; no continuous-model certification."""
from dataclasses import dataclass
import numpy as np
from scipy.special import j1, j0
import time, hashlib
from pathlib import Path

def verify_scalar_kernel(kernel):
    """Verify against the trusted scalar constructor, including FFT buffers.

    This guards accidental kernel/buffer/metadata changes. It is not a security
    boundary against monkey-patching Python or a formal floating-point proof.
    Regeneration cost is intentionally exposed instead of called free.
    """
    from .physics import FFTGreen
    started=time.perf_counter()
    if type(kernel) is not FFTGreen:
        return False, {'reason':'requires exact trusted FFTGreen constructor type'}
    reference=FFTGreen(kernel.g,kernel.k,cell_integrated=kernel.cell_integrated,
                       quadrature_order=kernel.quadrature_order)
    same=(kernel.n==reference.n and kernel.h==reference.h and kernel.npix==reference.npix
          and kernel.fft_shape==reference.fft_shape
          and np.array_equal(kernel.kernel,reference.kernel)
          and np.array_equal(kernel._fk,reference._fk)
          and np.array_equal(kernel._fk_adj,reference._fk_adj))
    return bool(same), {'kernel_verified':bool(same),'verification_seconds':time.perf_counter()-started,
                       'time_convention':'exp(-i omega t)','outgoing_kernel':'i H0^(1)/4',
                       'passive_contrast_sign':'Im chi >= 0',
                       'discretization':'uniform equal-square source cells, center collocation',
                       'inner_product':'Euclidean uniform-cell mass-white coordinates',
                       'self_term':'square polar integral' if kernel.cell_integrated else 'equal-area disk integral',
                       'constructor':'a2.physics.FFTGreen',
                       'constructor_source_sha256':hashlib.sha256(Path(__file__).with_name('physics.py').read_bytes()).hexdigest(),
                       'reason':'matched trusted construction' if same else 'kernel or FFT buffers differ from trusted construction'}

@dataclass
class Certificate:
    applicable: bool
    reason: str
    gamma_min: float
    corrected_field: np.ndarray | None = None
    corrected_bound: float | None = None
    raw_bound: float | None = None

def gamma_diagonal(chi,k,h):
    chi=np.asarray(chi,complex)
    if np.any(chi==0):return None,'zero contrast: eliminate known inactive cells before L formulation'
    a=h/np.sqrt(np.pi);delta=max(0.,k*k*h*h/4-np.pi*k*a*j1(k*a)/2)
    gamma=chi.imag/np.abs(chi)**2-delta
    return gamma,('applicable' if np.min(gamma)>0 else 'dissipation lower bound not positive')

def gamma_for_fft_kernel(chi, kernel):
    """Sufficient dissipativity bound for this implementation's scalar kernel.

    Tensor GL source quadrature has a nonnegative angular Fourier weight when
    k*h <= pi. Exact square self replacement contributes only a diagonal shift.
    Algebraic bound for the specified ideal kernel; not an interval-arithmetic
    bound on special functions, FFT roundoff, or continuum modeling error.
    """
    chi = np.asarray(chi, complex)
    verified, verification=verify_scalar_kernel(kernel)
    if not verified:return None, {'applicable':False,**verification}
    if np.any(chi == 0):
        return None, {**verification,'applicable': False, 'reason': 'eliminate known inactive zero-contrast cells'}
    if kernel.cell_integrated:
        if kernel.k*kernel.h > np.pi:
            return None, {**verification,'applicable': False, 'reason': 'sufficient quadrature positivity condition k*h<=pi not met'}
        nodes, weights = np.polynomial.legendre.leggauss(kernel.quadrature_order)
        x, y = np.meshgrid(nodes*kernel.h/2, nodes*kernel.h/2, indexing='ij')
        w = np.outer(weights, weights)*(kernel.h/2)**2
        reference_self = kernel.k**2/4*np.sum(w*j0(kernel.k*np.hypot(x,y)))
        label = 'tensor GL positive angular Gram plus actual square-self diagonal shift'
    else:
        reference_self = kernel.k**2*kernel.h**2/4
        label = 'point J0 angular Gram plus actual disk-self diagonal shift'
    actual_self = kernel.kernel[kernel.n-1,kernel.n-1].imag
    shift = float(actual_self-reference_self)
    gamma = chi.imag/np.abs(chi)**2 + shift
    return gamma, {**verification,'applicable': bool(np.min(gamma)>0), 'reason': label,
                   'self_shift': shift, 'reference_self_imag': float(reference_self),
                   'kh': float(kernel.k*kernel.h), 'gamma_min': float(np.min(gamma)),
                   'scope': 'specified discrete scalar kernel; floating-point diagnostics, not interval or continuum certification'}

def certify(L,E,C,J,Z,gamma):
    """Output-centric primal/adjoint residual bound in Euclidean mass-white coordinates.

    L may be a dense array or scipy LinearOperator. gamma must be established
    independently for this exact discrete kernel. E,J [N,Ls]; Z [N,M].
    """
    if gamma is None or np.min(gamma)<=0:return Certificate(False,'missing positive stability bound',float('nan') if gamma is None else float(np.min(gamma)))
    rp=E-L@J; rd=C.conj().T-L.conj().T@Z if isinstance(L,np.ndarray) else C.conj().T-L.H@Z
    corr=Z.conj().T@rp; fc=C@J+corr
    gp=rp/np.sqrt(gamma[:,None]);gd=rd/np.sqrt(gamma[:,None]);bound=float(np.linalg.norm(gd,2)*np.linalg.norm(gp,'fro'))
    return Certificate(True,'discrete corrected-output certificate',float(np.min(gamma)),fc,bound,float(bound+np.linalg.norm(corr,'fro')))

def galerkin(L,rhs,Q,adjoint=False):
    mat=L.conj().T if adjoint else L
    return Q@np.linalg.solve(Q.conj().T@mat@Q,Q.conj().T@rhs)

def structure_decision(old_lower,new_upper,budget,old_witness=None,new_lower=None):
    """Only certified class bounds/witnesses are eligible; local minima are not lower bounds."""
    if old_lower>budget and new_upper<=budget:return 'supported_new_structure'
    if new_lower is not None and old_lower>budget and new_lower>budget:return 'model_incompatible'
    if old_witness is not None and old_witness<=budget and new_upper<=budget:return 'both_compatible'
    return 'undetermined'
