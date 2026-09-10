"""Small, auditable fixed-point interval arithmetic. No floating-point certificate.

An interval [lo, hi] represents exact integer endpoints / 2**BITS. Every
operation rounds outwards using Python integer division and integer sqrt.
Only sin/cos on [-1,1] are required. Taylor polynomials have explicit Lagrange
remainders; no libm or third-party interval backend enters a proof.
"""
from dataclasses import dataclass
from math import factorial, isqrt

BITS = 160
SCALE = 1 << BITS

def ceildiv(a: int, b: int) -> int:
    return -((-a) // b)

@dataclass(frozen=True)
class I:
    lo: int
    hi: int

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('Reversed interval')

    @classmethod
    def rational(cls, p: int, q: int = 1):
        if not isinstance(p, int) or not isinstance(q, int) or q <= 0:
            raise TypeError('Exact integer numerator and positive denominator required')
        return cls(p * SCALE // q, ceildiv(p * SCALE, q))

    @classmethod
    def box(cls, lp: int, lq: int, hp: int, hq: int):
        a, b = cls.rational(lp, lq), cls.rational(hp, hq)
        return cls(a.lo, b.hi)

    @staticmethod
    def cast(a):
        if isinstance(a, I):
            return a
        if isinstance(a, int):
            return I.rational(a)
        raise TypeError('No implicit float inputs')

    def __add__(self, other):
        b = I.cast(other)
        return I(self.lo + b.lo, self.hi + b.hi)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-I.cast(other))

    def __rsub__(self, other):
        return I.cast(other) + (-self)

    def __mul__(self, other):
        b = I.cast(other)
        v = [self.lo*b.lo, self.lo*b.hi, self.hi*b.lo, self.hi*b.hi]
        return I(min(v)//SCALE, ceildiv(max(v), SCALE))
    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError('Interval contains zero')
        return I(SCALE*SCALE//self.hi, ceildiv(SCALE*SCALE, self.lo))

    def __truediv__(self, other):
        return self * I.cast(other).reciprocal()

    def __rtruediv__(self, other):
        return I.cast(other) * self.reciprocal()

    def __pow__(self, n):
        if not isinstance(n, int) or n < 0:
            raise ValueError('Only nonnegative integer powers')
        if n == 0:
            return I.rational(1)
        if n % 2 == 0:
            lo = 0 if self.lo <= 0 <= self.hi else min(self.lo**n,self.hi**n)
            hi = max(self.lo**n,self.hi**n)
        else:
            lo,hi = self.lo**n,self.hi**n
        den = SCALE**(n-1)
        return I(lo//den, ceildiv(hi,den))

    def sqrt(self):
        if self.lo < 0:
            raise ValueError('Negative square root')
        l, h = isqrt(self.lo*SCALE), isqrt(self.hi*SCALE)
        return I(l, h + (h*h < self.hi*SCALE))

    def abs_upper(self):
        return max(abs(self.lo), abs(self.hi))

    def trig(self, cosine=False, terms=20):
        # Polynomial degree 2*terms-1 (sin) or 2*terms-2 (cos).
        # The first omitted Taylor coefficient may vanish; the stated generic
        # degree+1 Lagrange remainder is used, not an alternating-series guess.
        if self.abs_upper() > SCALE:
            raise ValueError('Trig certificate only supports [-1,1]')
        p = I.rational(0)
        for n in range(terms):
            exponent = 2*n + (0 if cosine else 1)
            p += I.rational((-1)**n, factorial(exponent)) * self**exponent
        power = 2*terms - (1 if cosine else 0)
        m = self.abs_upper()
        radius = ceildiv(m**power, SCALE**(power-1)*factorial(power))
        return I(p.lo-radius,p.hi+radius)

    def sin(self):
        return self.trig(False)

    def cos(self):
        return self.trig(True)

    def decimal_bounds(self, digits=24):
        # Both strings are outward roundings, not nearest formatting.
        scale = 10**digits
        low = self.lo*scale//SCALE
        high = ceildiv(self.hi*scale,SCALE)
        def fmt(v):
            sign = '-' if v < 0 else ''
            v = abs(v)
            return f'{sign}{v//scale}.{v%scale:0{digits}d}'
        return [fmt(low),fmt(high)]

    def record(self):
        return {'lo_dyadic':str(self.lo),'hi_dyadic':str(self.hi),
                'denominator_power_of_two':BITS,'outward_decimal':self.decimal_bounds()}

    def diagnostic_midpoint(self):
        return (self.lo+self.hi)/(2*SCALE)


def spherical1(z):
    s,c = z.sin(),z.cos()
    j = s/z**2-c/z
    jp = s/z+2*c/z**2-2*s/z**3
    jpp = -2*jp/z-(1-2/z**2)*j
    dj = jp+j/z
    djp = -jp/z-j+j/z**2
    djpp = -jpp/z-jp+2*jp/z**2-2*j/z**3
    y = -c/z**2-s/z
    yp = -c/z+2*s/z**2+2*c/z**3
    dy = yp+y/z
    dyp = -yp/z-y+y/z**2
    return j,jp,jpp,dj,djp,djpp,y,yp,dy,dyp


def modal(eps, x):
    """Exact-real enclosures of q, h and their eps/x derivatives, lossless."""
    m = eps.sqrt(); u=m*x
    j,jp,jpp,dj,djp,djpp,_,_,_,_ = spherical1(u)
    jx,jpx,_,dx,dpx,_,yx,ypx,dyx,dpyx = spherical1(x)
    n = m*j*dx-jx*dj
    d = m*j*dyx-yx*dj
    mp = 1/(2*m); up = x*mp
    np = (mp*j+m*jp*up)*dx-jx*djp*up
    dp = (mp*j+m*jp*up)*dyx-yx*djp*up
    q=n/d
    qp=(np*d-n*dp)/d**2
    nx=m*m*jp*dx+m*j*dpx-jpx*dj-jx*djp*m
    ddx=m*m*jp*dyx+m*j*dpyx-ypx*dj-yx*djp*m
    qx=(nx*d-n*ddx)/d**2
    den=(1+q**2).sqrt()
    h=q/den
    hp=qp/den**3
    hx=qx/den**3
    return {'q':q,'qprime':qp,'denominator':d,'h':h,'hprime':hp,'hx':hx}
