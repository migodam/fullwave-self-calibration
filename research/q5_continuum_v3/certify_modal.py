"""Exact-integer inclusion certificate for the A5 lossless electric-dipole class.

No binary floating-point arithmetic is used in the certificate. A number is a
closed interval with endpoints integers / 2**BITS. All operations round outward.
Special functions are bounded Taylor polynomials with an explicit absolute tail.
The reported endpoint minima follow from a proved, covered concavity statement.
Run: python research/q5_window_v2/certify_modal.py --output PATH
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from math import factorial, isqrt
from pathlib import Path
import argparse
import hashlib
import json
import time

BITS = 192
SCALE = 1 << BITS
TERMS = 16


def ceildiv(a: int, b: int) -> int:
    return -((-a) // b)


@dataclass(frozen=True)
class I:
    lo: int
    hi: int

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('Inverted interval')

    @staticmethod
    def q(n=0, d=1):
        f = n if isinstance(n, Fraction) else Fraction(n, d)
        return I((f.numerator*SCALE)//f.denominator,
                 ceildiv(f.numerator*SCALE, f.denominator))

    @staticmethod
    def box(a, b):
        return I(I.q(Fraction(a)).lo, I.q(Fraction(b)).hi)

    def __add__(self, other):
        b = other if isinstance(other, I) else I.q(other)
        return I(self.lo+b.lo, self.hi+b.hi)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -(other if isinstance(other, I) else I.q(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        b = other if isinstance(other, I) else I.q(other)
        p = [self.lo*b.lo, self.lo*b.hi, self.hi*b.lo, self.hi*b.hi]
        return I(min(p)//SCALE, ceildiv(max(p), SCALE))
    __rmul__ = __mul__

    def __truediv__(self, other):
        b = other if isinstance(other, I) else I.q(other)
        if b.lo <= 0 <= b.hi:
            raise ZeroDivisionError('Interval divisor contains zero')
        pairs = [(a*SCALE, d) for a in (self.lo, self.hi) for d in (b.lo, b.hi)]
        return I(min(a//d for a,d in pairs), max(ceildiv(a,d) for a,d in pairs))

    def __rtruediv__(self, other):
        return I.q(other)/self

    def __pow__(self, n):
        if not isinstance(n, int) or n < 0:
            raise ValueError('Only nonnegative integer powers')
        if n == 0:
            return I.q(1)
        if n % 2 == 0 and self.lo < 0 < self.hi:
            return I(0, max(abs(self.lo), abs(self.hi)))**n
        result = I.q(1)
        for _ in range(n):
            result = result*self
        return result

    def sqrt(self):
        if self.lo < 0:
            raise ValueError('Negative square-root endpoint')
        a, b = isqrt(self.lo*SCALE), isqrt(self.hi*SCALE)
        return I(a, b if b*b == self.hi*SCALE else b+1)

    def contains(self, x):
        f = Fraction(x)
        return self.lo*f.denominator <= f.numerator*SCALE <= self.hi*f.denominator

    def export(self, digits=35):
        # Decimal strings are rounded OUTWARD, not nearest-rounded float exports.
        p = 10**digits
        def dec(v):
            sign = '-' if v < 0 else ''
            v = abs(v)
            return f'{sign}{v//p}.{v%p:0{digits}d}'
        return {'lower':dec(self.lo*p//SCALE), 'upper':dec(ceildiv(self.hi*p,SCALE)),
                'lower_integer':str(self.lo), 'upper_integer':str(self.hi),
                'denominator_power_of_two':BITS}


def coef(kind, n):
    if kind == 'J':
        return Fraction((-1)**n*(2*n+2), factorial(2*n+3))
    if kind == 'D':
        return Fraction((-1)**n*(2*n+2)**2, factorial(2*n+3))
    if kind == 'C':
        return Fraction((-1)**n, factorial(2*n))
    if kind == 'S':
        return Fraction((-1)**n, factorial(2*n+1))
    raise ValueError(kind)


def series(kind, t, derivative=0):
    """Entire functions J, D, cos(sqrt(t)), sinc(sqrt(t)).

    Domain 0<=t<=1, derivative<=2. For n>=16, successive absolute Taylor
    terms of any implemented derivative have ratio <1/2. Therefore twice
    the first omitted term bounds the complete tail, including at t=0.
    """
    if not (0 <= t.lo <= t.hi <= SCALE) or derivative not in (0,1,2):
        raise ValueError('Taylor proof domain violated')
    def dc(n):
        return coef(kind,n)*Fraction(factorial(n), factorial(n-derivative))
    value = I.q(0)
    for n in range(TERMS-1, derivative-1, -1):
        value = value*t + I.q(dc(n))
    first = I.q(abs(dc(TERMS))) * I(t.hi,t.hi)**(TERMS-derivative)
    tail = 2*first.hi
    return value + I(-tail,tail)


def modal(e, x):
    """q and its first two real-permittivity derivatives, without sqrt(e)."""
    s = x*x
    t = e*s
    J, D = series('J',t), series('D',t)
    D1, D2 = series('D',t,1), series('D',t,2)
    J0, D0 = series('J',s), series('D',s)
    c, sinc = series('C',s), series('S',s)
    Y, Z = c+s*sinc, (1-s)*c+s*sinc
    U = e*J
    A, B = U*D0-J0*D, U*Z+Y*D
    # J(s) Z(s)+D(s) Y(s)=1, and d[e J(e s)]/de=D(e s)/2.
    W = D*D/2-t*J*D1
    W1 = s*(D*D1/2-t*J*D2)
    B1 = Z*D/2+Y*s*D1
    q = x**3*A/B
    q1 = x**3*W/(B*B)
    q2 = x**3*(W1*B-2*W*B1)/(B**3)
    h = q/(1+q*q).sqrt()
    h1 = q1/((1+q*q)*(1+q*q).sqrt())
    return dict(J=J,D=D,D1=D1,D2=D2,Y=Y,Z=Z,A=A,B=B,W=W,W1=W1,B1=B1,
                q=q,q1=q1,q2=q2,h=h,h1=h1)


def certify():
    start = time.perf_counter()
    domains = [('1/5','3/2','4'), ('3/20','2','5')]
    rows=[]
    for xstr, low, high in domains:
        x=I.q(Fraction(xstr)); full=modal(I.box(low,high),x)
        left=modal(I.q(Fraction(low)),x); right=modal(I.q(Fraction(high)),x)
        required_positive=('J','D','D2','Y','Z','B','W','B1')
        required_negative=('D1','W1','q2')
        for key in required_positive:
            if full[key].lo <= 0:
                raise ArithmeticError(f'Failed strict positive sign: {key}')
        for key in required_negative:
            if full[key].hi >= 0:
                raise ArithmeticError(f'Failed strict negative sign: {key}')
        if left['q'].lo <= 0 or right['h1'].lo <= 0:
            raise ArithmeticError('Endpoint sign failed')
        # With q>0, q'>0 and q''<0, h''<0 on the ENTIRE interval.
        # Hence max h=h(high), min h'=h'(high); no midpoint or grid inference.
        rows.append({'size_parameter_exact':xstr, 'material_interval_exact':[low,high],
            'covered_boxes':1, 'coverage':'One rational box equals the full declared interval',
            'global_sign_enclosures':{k:full[k].export() for k in required_positive+required_negative},
            'h_min':left['h'].export(), 'h_max':right['h'].export(),
            'minimum_derivative':right['h1'].export(),
            'proof':'q>0, q1>0, q2<0 => h1>0 and h2<0; h1 minimum at upper endpoint'})
    return {'complete':True, 'backend':'exact Python integer fixed-point interval arithmetic',
        'bits':BITS, 'terms':TERMS, 'floating_point_used_for_enclosures':False,
        'size_parameters_are_exact_rationals':True,
        'truncation':'Absolute tail <=2 times first omitted absolute term; ratio<1/2 for n>=16,0<=t<=1,d<=2',
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'wall_seconds_diagnostic_only':time.perf_counter()-start,'rows':rows}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,
        default=Path(__file__).parent/'results/modal_certificate.json')
    a=p.parse_args(); data=certify()
    if a.output.exists():
        raise FileExistsError('Refusing to overwrite evidence; choose a new output path')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(data,indent=2)+'\n')
    for row in data['rows']:
        print(row['size_parameter_exact'], 'm >=',row['minimum_derivative']['lower'],
              'H <=',row['h_max']['upper'])

if __name__=='__main__':
    main()
