"""Do not transfer a coarse geometry certificate into a fine modal readout."""
from pathlib import Path
import json,hashlib
from dyadic_interval import I, SCALE, modal
HERE=Path(__file__).resolve().parent

def transfer(z):
    # [h_1(z)/z]/h_1(1) = exp(i(z-1))*(z+i)/(1+i)/z^3.
    d=z-1;c=d.cos();s=d.sin()
    re=((z+1)*c-(1-z)*s)/(2*z**3)
    im=((z+1)*s+(1-z)*c)/(2*z**3)
    return ((re-1)**2+im**2).sqrt()

def run():
    out=HERE/'results/geometry_readout_budget.json'
    if out.exists():raise FileExistsError(out)
    endpoints=[{'kR':[p,100],'relative_electric_dipole_transfer_error':transfer(I.rational(p,100)).record()} for p in [98,102]]
    z=I.box(9997,10000,10003,10000)
    amplitude=((z*z+1)/2).sqrt()/z**3
    logderiv=((3/z-z/(z*z+1))**2+(z*z/(z*z+1))**2).sqrt()
    lipschitz=amplitude*logderiv
    bound=lipschitz*I.rational(3,10000)
    ha=modal(I.rational(47,10),I.rational(3,20))['h']
    hb=modal(I.rational(491,100),I.rational(3,20))['h']
    gain=ha/hb
    assert gain.lo>I.rational(3,4).hi and gain.hi<I.rational(5,4).lo
    result={'kind':'geometry_to_readout_gap_audit','status':'not_a_joint_recovery_certificate',
      'coarse_geometry_error_target_kR':[2,100],'endpoint_transfer_errors':endpoints,
      'fine_geometry_halfwidth_kR':[3,10000],'dipole_transfer_derivative_enclosure':lipschitz.record(),
      'fine_halfwidth_dipole_transfer_error_bound':bound.record(),
      'material_ambiguity_in_class_G':{'epsilon_A':[47,10],'epsilon_B':[491,100],
        'x':[3,20],'gain_A':1,'gain_B_modulus':gain.record(),
        'gain_B_complex_definition':'a_1(epsilon_A)/a_1(epsilon_B)',
        'allowed_gain_annulus':[.75,1.25]},
      'limitations':['Only electric dipole radial transfer is bounded here',
        'The coarse G guarantee cannot be inserted as a 0.1-percent M readout calibration',
        'No claim that observed geometry errors attain the worst-case bound',
        'No joint G-to-M or simultaneous class C recovery theorem'],
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print([(r['kR'],r['relative_electric_dipole_transfer_error']['outward_decimal']) for r in endpoints])
    print('fine bound',bound.decimal_bounds(),'gain',gain.decimal_bounds())
if __name__=='__main__':run()
