"""Exact class-C coercivity and residual-to-observation constants.
This does not supply the missing continuum residual of any discrete solver.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from certify_modal import decimal_bound

def constants():
    k=F(18);a,b=F(7,200),F(1,40);d=F(1167,10000)
    assert d*d < F(115,1000)**2+F(1,50)**2
    assert (2*k*a)**2<2 and (2*k*b)**2<2
    diagonal=[F(3)/(9+F(3,100)**2)-k*k*a*a/2,
              F(4)/(16+F(1,20)**2)-k*k*b*b/2]
    cross2=k**4*(a*b)**3/(4*(d-a)*(d-b))
    cross=F(97,2000);alpha=F(43,500)
    assert cross2<cross*cross and min(diagonal)-cross>alpha
    # All receivers have nominal radius .6, shift <= .002.
    # Sphere 1: .6-.002-.06-.035=.503.
    # Sphere 2: norm(.055,.02)<.06, so .513>.503.
    assert F(55,1000)**2+F(2,100)**2 < F(6,100)**2
    R=F(503,1000);pi_lower=F(314159,100000);N=12
    # S^2 <= N*(a^3+b^3)/(4*pi) * kernel_factor^2.
    C=k*k/R+2*k/(R*R)+2/(R**3)
    S2=F(N)*(a**3+b**3)/(4*pi_lower)*C*C
    assert S2<36
    gain=F(5,4);factor=gain*6/alpha
    assert factor<88
    return {'scope':'entire declared class C; distributional VIE on L2(D)^3',
        'coercivity_lower':str(alpha),'inverse_norm_upper':str(1/alpha),
        'diagonal_lower':[decimal_bound(v) for v in diagonal],
        'cross_block_norm_upper':str(cross),'minimum_source_receiver_distance':'503/1000',
        'receiver_operator_norm_upper':'6','gain_times_receiver_over_alpha_upper':str(factor),
        'simplified_total_field_error_bound':'88 * ||E_inc - (chi^-1-G_k) p_tilde||_L2',
        'actual_continuum_residual_bound':None,
        'finite_material_separation_bound':None,
        'status':'operator constant certified; class-C recovery certificate unresolved'}

if __name__=='__main__':
    p=Path(__file__).resolve().parents[1]/'results/continuum_bound.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(constants(),indent=2)+'\n');print(p.read_text())
