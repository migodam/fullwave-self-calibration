"""Independent numerical cross-check of reported Q5 M2 endpoint arithmetic.

This is NOT interval coverage, a Maxwell model-error bound, or a new test stream.
Run: python check_reported_endpoints.py
"""
from pathlib import Path
import json
import mpmath as mp

mp.mp.dps = 70

def h(e, x):
    m = mp.sqrt(e)
    def j(z): return mp.sin(z)/z**2-mp.cos(z)/z
    def dj(z): return mp.sin(z)/z+mp.cos(z)/z**2-mp.sin(z)/z**3
    def y(z): return -mp.cos(z)/z**2-mp.sin(z)/z
    def dy(z): return -mp.cos(z)/z+mp.sin(z)/z**2+mp.cos(z)/z**3
    numerator = m*j(m*x)*dj(x)-j(x)*dj(m*x)
    denominator = m*j(m*x)*dy(x)-y(x)*dj(m*x)
    q = numerator/denominator
    return q/mp.sqrt(1+q*q)

def main():
    output = {'status': 'numerical arithmetic cross-check only',
              'precision_decimal_digits': mp.mp.dps, 'channels': []}
    aa=mp.mpf(3)/4
    delta=mp.mpf(1)/10
    Ba=mp.mpf('0.003291')
    b=mp.mpf('1e-6')
    sigmas=[mp.mpf('2.6664081e-6'),mp.mpf('1.8047986e-6')]
    for i,(x,U) in enumerate([(mp.mpf(1)/5,mp.mpf(4)),(mp.mpf(3)/20,mp.mpf(5))]):
        high=h(U,x); low=h(U-2*delta,x)
        C=aa*(high-low)-2*Ba*low
        def prop(z): return mp.sqrt(z**(-2)-z**(-4)+z**(-6))
        left,right=mp.mpf(1),mp.mpf(100)
        for _ in range(240):
            mid=(left+right)/2
            if prop(mid)*C >= 2*b: left=mid
            else: right=mid
        allowances={}
        for repeats,drift in [(1,mp.mpf(0)),(16,mp.mpf('0.0005'))]:
            bn=sigmas[i]*mp.sqrt(mp.log(1000))/mp.sqrt(repeats)
            br=mp.mpf('3.290526731492')*mp.mpf('0.001')/mp.sqrt(repeats)+drift
            allowances[str(repeats)]=str((aa*(high-low)-2*br*low)/2-bn)
        output['channels'].append({'size':str(x),'h_upper_endpoint':str(high),
            'h_endpoint_minus_0_2':str(low),'finite_C_at_Ba_0_003291':str(C),
            'numerical_radial_boundary':str((left+right)/2),
            'model_bias_allowance':allowances})
    path=Path(__file__).with_name('endpoint_crosscheck.json')
    path.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))

if __name__=='__main__':main()
