"""Finite spherical modal readout: rational tail/budget certificate.

The physical instrument-error caps are declared acceptance assumptions, not
measured accuracies. This applies ONLY to separately read lossless spheres.
See MANUSCRIPT_ZH.md for the all-order tail and quadrature proof.
"""
from fractions import Fraction as F
from pathlib import Path
import json, math, hashlib
from interval_certificate import I, SCALE, series, decimal
HERE=Path(__file__).resolve().parent


def odd_df(n):
    return math.prod(range(1,n+1,2))


def tail_term(ell,x):
    # Upper bound for the entire tangential l-th outgoing field at kR=2.
    return 2*x**(2*ell+1)/((2**ell)*odd_df(2*ell-1))*(8+10*ell+F(40,2*ell-1))


def tail_bound(first_ell,x):
    assert first_ell>=2 and 0<x<=F(1,5)
    # T_(l+1)/T_l <= x^2/(2(2l+1))*(l+1)/l < 1/100.
    return F(100,99)*tail_term(first_ell,x)


def magnetic_cap(a,b,x):
    xi=I.of(x);t=xi*xi; Jt,Dt=series('J',t),series('D',t)
    C,S=series('C',t),series('S',t)
    Y=-C-t*S;Z=(1-t)*C+t*S
    cap=0;step=F(1,1000);count=int((b-a)/step)
    for j in range(count):
        eps=I.of(a+j*step,a+(j+1)*step);s=eps*t
        J,D=series('J',s),series('D',s)
        den=J*Z-Y*D
        assert den.lo>0
        qb=xi**3*(J*Dt-Jt*D)/den
        cap=max(cap,abs(qb.lo),abs(qb.hi))
    return F(cap,SCALE),count


def exp_lower(t,n=70):
    return sum((t**j/F(math.factorial(j)) for j in range(n)),F(0))


def main(output=None):
    target=Path(output) if output else HERE/'results/readout_certificate.json'
    if target.exists():raise FileExistsError('Immutable output exists')
    cert=json.loads((HERE/'results/modal_certificate.json').read_text())
    # Simple exact sufficient inequalities behind the Gaussian union bound.
    assert exp_lower(F(9))>8000 and exp_lower(F(49,8))>450
    failure_bound=F(2,8000)+F(8,35*450)
    assert failure_bound<F(1,1000)
    eta=F('0.0005');refbias=F('0.00025');other=F('0.000002')
    ba=F('0.001')*F(7,2)/4+refbias
    rows=[]
    for index,(a,b,x) in enumerate([(F(3,2),F(4),F(1,5)),(F(2),F(5),F(3,20))]):
        summary=cert['summaries'][index]
        H=F(summary['amplitude_cap_upper']);m=F(summary['endpoint_slope_lower'])
        magnetic,count=magnetic_cap(a,b,x)
        alias=F(20,9)*tail_bound(5,x)
        Cread=3*H+4*magnetic+F(20,9)*tail_bound(2,x)
        sigma=[F('2.6664081e-6'),F('1.8047986e-6')][index]
        bn=3*sigma/4
        bmodel=F(5,4)*(alias+eta*Cread)+other
        err=(bmodel+bn+H*ba)/((F(3,4)-ba)*m)
        margin=F(1,10)*(F(3,4)-ba)*m-(bmodel+bn+H*ba)
        assert err<F(1,10) and margin>0
        export=lambda v:I.of(v).export()
        rows.append({'object':index+1,'size':str(x),'magnetic_boxes':count,
            'magnetic_dipole_cap_upper':export(magnetic),
            'quadrature_alias_bound':export(alias),'relative_readout_sensitivity_bound':export(Cread),
            'instrument_model_bias_bound':export(bmodel),
            'noise_bound_after_16_repeats':export(bn),
            'reference_total_error_bound':export(ba),
            'uniform_material_error_bound':export(err),'remaining_margin':export(margin),
            'per_tangential_component_complex_sigma_sufficient':export(2*sigma)})
    result={'status':'conditional_finite_Maxwell_readout_certificate',
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'separately read lossless known-radius spheres, NOT simultaneous class C',
        'directions':18,'tangential_complex_components_per_direction':2,'repetitions':16,
        'total_complex_field_samples':1152,'total_real_reference_samples':16,
        'squared_noise_row_norm_exact':'224/1053','uniform_failure_probability_upper':str(failure_bound),
        'assumed_relative_channel_operator_error':str(eta),'assumed_reference_absolute_bias':str(refbias),
        'assumed_other_projected_model_error':str(other),
        'hardware_status':'specifications not experimentally verified; no lab attainability claim',
        'rows':rows}
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--output');main(p.parse_args().output)
