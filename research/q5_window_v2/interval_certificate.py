"""Exact-integer interval certificate for the inherited lossless modal class.

Only Python's integer/Fraction arithmetic is trusted; no floating-point special
function or decimal-to-binary size-parameter conversion enters the certificate.
All intervals live on the outward-rounded 10**(-36) lattice.  The final slope
uses (1+q*q)**2 >= (1+q*q)**(3/2), avoiding uncertified roots.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import time

DIGITS = 36
SCALE = 10 ** DIGITS
TERMS = 12


def ceildiv(a: int, b: int) -> int:
    if b <= 0:
        raise ValueError('positive denominator required')
    return -((-a) // b)


@dataclass(frozen=True)
class I:
    lo: int
    hi: int

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('empty interval')

    @classmethod
    def of(cls, lo, hi=None):
        lo = F(lo)
        hi = lo if hi is None else F(hi)
        return cls((lo.numerator * SCALE) // lo.denominator,
                   ceildiv(hi.numerator * SCALE, hi.denominator))

    @staticmethod
    def coerce(x):
        return x if isinstance(x, I) else I.of(x)

    def __add__(self, other):
        other = self.coerce(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        v = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return I(min(v) // SCALE, ceildiv(max(v), SCALE))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError('interval denominator contains zero')
        if self.hi < 0:
            return -((-self).reciprocal())
        return I(SCALE * SCALE // self.hi, ceildiv(SCALE * SCALE, self.lo))

    def __truediv__(self, other):
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def __pow__(self, n):
        if not isinstance(n, int) or n < 0:
            raise ValueError('only nonnegative integer powers supported')
        out = I.of(1)
        for _ in range(n):
            out = out * self
        return out

    def export(self):
        return [decimal(self.lo), decimal(self.hi)]


def decimal(x):
    sign = '-' if x < 0 else ''
    x = abs(x)
    return f'{sign}{x // SCALE}.{x % SCALE:0{DIGITS}d}'


def coefficient(kind, n):
    if kind == 'J':
        return F((-1)**n * (2*n+2), math.factorial(2*n+3))
    if kind == 'D':
        return F((-1)**n * (2*n+2)**2, math.factorial(2*n+3))
    if kind == 'C':
        return F((-1)**n, math.factorial(2*n))
    if kind == 'S':
        return F((-1)**n, math.factorial(2*n+1))
    raise ValueError(kind)


def series(kind, s, derivative=0):
    """Polynomial plus a rigorous symmetric tail; requires 0<=s<=1, r<=2.

    For n>=12, every successive absolute derivative-series term is at most
    half its predecessor.  The sum of the omitted tail is therefore bounded
    by twice its first absolute term, evaluated at the interval upper end.
    """
    if not 0 <= s.lo <= s.hi <= SCALE or derivative not in (0, 1, 2):
        raise ValueError('outside proved series domain')
    coeffs = [coefficient(kind, n) * math.prod(range(n-derivative+1, n+1))
              for n in range(derivative, TERMS)]
    out = I.of(0)
    for a in reversed(coeffs):
        out = out * s + a
    first = abs(coefficient(kind, TERMS)) * math.prod(
        range(TERMS-derivative+1, TERMS+1))
    tail = 2 * first * F(s.hi, SCALE) ** (TERMS-derivative)
    return out + I.of(-tail, tail)


def modal_intervals(eps, x):
    """Returns q,q',q'',h' lower enclosure and original Mie denominator."""
    x = I.of(x)
    t = x*x
    s = eps*t
    J, D = series('J', s), series('D', s)
    J1, D1 = series('J', s, 1), series('D', s, 1)
    J2, D2 = series('J', s, 2), series('D', s, 2)
    Jt, Dt = series('J', t), series('D', t)
    C, S = series('C', t), series('S', t)
    Y = -C-t*S
    Z = (1-t)*C+t*S
    n = eps*J*Dt-Jt*D
    d = eps*J*Z-Y*D
    np = (J+eps*t*J1)*Dt-Jt*t*D1
    dp = (J+eps*t*J1)*Z-Y*t*D1
    npp = (2*t*J1+eps*t*t*J2)*Dt-Jt*t*t*D2
    dpp = (2*t*J1+eps*t*t*J2)*Z-Y*t*t*D2
    q = x**3*n/d
    qp = x**3*(np*d-n*dp)/(d*d)
    qpp = x**3*((npp*d-n*dpp)*d-2*(np*d-n*dp)*dp)/(d**3)
    # This is a rigorous lower bound interval for h', NOT h' itself.
    slope_lower_expression = qp/(1+q*q)**2
    return q, qp, qpp, slope_lower_expression, d/(x*x)


def certificate(output, step=F(1,1000)):
    start = time.perf_counter()
    domains = [(F(3,2), F(4), F(1,5)), (F(2), F(5), F(3,20))]
    summaries, rows = [], []
    for index, (a,b,x) in enumerate(domains):
        count_fraction = (b-a)/step
        if count_fraction.denominator != 1:
            raise ValueError('step must exactly tile domain')
        count = int(count_fraction)
        min_q = min_qp = min_slope = min_den = None
        max_q = max_qpp = None
        previous = a
        for k in range(count):
            left, right = a+k*step, a+(k+1)*step
            assert left == previous
            previous = right
            q,qp,qpp,slope,den = modal_intervals(I.of(left,right),x)
            assert q.lo > 0 and qp.lo > 0 and den.lo > 0
            # Negative q'' and positive q,q' imply h''<0.
            assert qpp.hi < 0
            min_q = q.lo if min_q is None else min(min_q,q.lo)
            min_qp = qp.lo if min_qp is None else min(min_qp,qp.lo)
            min_slope = slope.lo if min_slope is None else min(min_slope,slope.lo)
            min_den = den.lo if min_den is None else min(min_den,den.lo)
            max_q = q.hi if max_q is None else max(max_q,q.hi)
            max_qpp = qpp.hi if max_qpp is None else max(max_qpp,qpp.hi)
            rows.append({'object':index+1,'eps':[str(left),str(right)],
                         'q':q.export(),'q_prime':qp.export(),
                         'q_second':qpp.export(),'h_prime_lower_expr':slope.export(),
                         'denominator':den.export()})
        assert previous == b
        # Justified endpoint use AFTER the full q'>0,q''<0 coverage above.
        q,qp,qpp,slope,den = modal_intervals(I.of(b),x)
        safe_slope = F(slope.lo,SCALE)
        safe_cap = F(q.hi,SCALE)  # h<=q; no uncertified square root.
        summaries.append({'object':index+1,'domain':[str(a),str(b)],
          'size_parameter_exact':str(x),'boxes':count,'tile_step_exact':str(step),
          'min_q_lower':decimal(min_q),'min_q_prime_lower':decimal(min_qp),
          'max_q_second_upper':decimal(max_qpp),'min_denominator_lower':decimal(min_den),
          'direct_box_slope_lower':decimal(min_slope),
          'endpoint_slope_lower':decimal(slope.lo),'amplitude_cap_upper':decimal(q.hi),
          'published_slopes_pass':safe_slope >= [F('0.000429'),F('0.000129')][index],
          'published_caps_pass':safe_cap <= [F('0.00272'),F('0.00131')][index],
          'h_strictly_increasing':True,'h_strictly_concave':True})
    body={'status':'certified_under_exact_integer_implementation',
          'arithmetic':'integer fixed-point interval, outward floor/ceil at every operation',
          'decimal_places':DIGITS,'series_terms':TERMS,
          'series_tail':'2 times first omitted absolute derivative-series term; 0<=s<=1; r<=2',
          'coverage':'closed rational boxes exactly tile both intervals, adjacent endpoints equal',
          'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'summaries':summaries,'box_count':len(rows),'boxes':rows,
          'wall_seconds':time.perf_counter()-start,
          'limits':['not a hardware certificate','not class C','not a proof-assistant verification']}
    output=Path(output)
    if output.exists():
        raise FileExistsError('certificate output is immutable; choose a new filename')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(body,indent=2)+'\n')
    print(json.dumps({'summaries':summaries,'boxes':len(rows),
                      'wall_seconds':body['wall_seconds']},indent=2))
    return body


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',default=str(Path(__file__).parent/'results/modal_certificate.json'))
    args=p.parse_args()
    certificate(args.output)
