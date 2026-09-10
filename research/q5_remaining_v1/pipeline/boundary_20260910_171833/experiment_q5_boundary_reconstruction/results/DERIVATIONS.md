# DERIVATIONS — formulas used in increment 2 (A5 sections 5-8)

Notation: `eps` = dielectric contrast, `x = ka`, `m = sqrt(eps)`, `z = m x`, `j = j1`, `y = y1`, `D_j(z) = (j1(z) + z j1'(z))/z`, `D_y` similarly.

## 0. Primitives (exact at mp.dps=50; the same expressions inside mpmath.iv at dps=35)
```
j1(z)  = sin z / z^2 - cos z / z          j1'(z) = sin z / z + 2 cos z / z^2 - 2 sin z / z^3
y1(z)  = -cos z / z^2 - sin z / z         y1'(z) = -cos z / z + 2 sin z / z^2 + 2 cos z / z^3
D_j(z) = sin z / z + cos z / z^2 - sin z / z^3
D_y(z) = -cos z / z + sin z / z^2 + cos z / z^3
D_j'(z)  = -j1'(z)/z - j1(z) + j1(z)/z^2        (also = -D_y(z) + 2 D_j(z)/z)
D_j''(z) = -j1''(z)/z - 2 j1'(z)/z + 2 j1'(z)/z^2 - 2 j1(z)/z^3
j1''(z)  = -(2/z) j1'(z) - (1 - 2/z^2) j1(z)     (spherical Bessel ODE, n = 1)
```
## 1. Reactance and first derivative (A5 section 8)
```
N(eps,x)   = m j1(z) D_j(x) - j1(x) D_j(z)
Den(eps,x) = m j1(z) D_y(x) - y1(x) D_j(z)          q(eps,x) = N / Den
```
Chain rule (`m = sqrt(eps)`, `m' = 1/(2m)`, `z' = dz/deps = x/(2m)`, `z'' = -x/(4 m^3)`):
```
A  = m j1(z)      A' = j1(z)/(2m) + (x/2) j1'(z)
N'   = A' D_j(x) - j1(x) D_j'(z) z'
Den' = A' D_y(x) - y1(x) D_j'(z) z'
q'   = (N' Den - N Den') / Den^2                                          [V1, used by mie_reactance]
```
Equivalent arrangements of the same scalar function (differing only as interval enclosures):
```
V2: replace D_j'(z) by its directly differentiated closed form (above).
V3: q' = q * (N'/N - Den'/Den).
V4: W = j1(x) D_y(x) - y1(x) D_j(x) (constant in eps);
    q' = W * ( D_j(z) A' - (x/2) j1(z) D_j'(z) ) / Den^2     (no N/Den interval products).
V5: V1 at iv.dps = 25 and 50.   V7: substitute eps = mid + delta as a small interval, then V1.
V6: q'(box midpoint, dps=50) - c * (halfwidth) * |q''|_upper(box),  c = 0.5 and 1.0.
```
## 2. Second derivative q'' (used only to localise the infimum)
```
A'' = -j1(z)/(4 m^3) + x m'^2 j1'(z) + (x^2/2) m' j1''(z)      with m' = 1/(2m)
d2/de2 [ D_j(z) ] = D_j''(z) z'^2 + D_j'(z) z''
N''   = A'' D_j(x) - j1(x) d2/de2[D_j(z)]        Den'' = A'' D_y(x) - y1(x) d2/de2[D_j(z)]
F  = N' Den - N Den' ;  q'' = [ (N'' Den - N Den'') Den - 2 F Den' ] / Den^3
```
`q'' < 0` on both intervals => `q'` strictly decreasing => `inf q'` at the right endpoint.
## 3. Amplitude and its derivative
```
h(eps,x) = |t| = q / sqrt(1 + q^2)           h'(eps,x) = q' / (1 + q^2)^{3/2}
```
`h` is strictly increasing where `q' > 0`, so `h^{-1}` exists on `I_i` and is Lipschitz with constant `<= 1/m_i`, `m_i = inf_{I_i} h'`.
## 4. Section 5: two-world modal difference
```
t = i q / (1 - i q)   (energy conservation Re t = -|t|^2)      q = -i t / (1 + t)
t(eps2) = t(c q1),  g2 = g1 / c
delta_mu = g2 t(c q1) - g1 t(q1) = -( g1 (c-1) q1^2 ) / ( (1 - i c q1)(1 - i q1) )
|delta_mu| / |g1 t(q1)| = |c-1| |q1| / sqrt(1 + c^2 q1^2)  <=  |c-1| |q1|
```
(Linear terms cancel in the common denominator; the survivor is `i q1 * i(c-1) q1 = -(c-1) q1^2`.)
## 5. Section 6: binary risk
```
D^2 = sum_i |delta_mu_i|^2 / sigma_i^2          (Sigma = diag(sigma_i^2), proper complex)
whitened projection: mean gap D, variance 1/2  =>  p_* = Phi(-D/sqrt2) = 0.5 * erfc(D/2)
Equal-prior optimal test: decide world 1 iff 2 Re sum_i conj(delta_i)(Y_i - mu_0,i)/sigma_i^2 > D^2.
```
## 6. Section 7: uniform material error bound
```
Z = a + xi,  a = |g| >= a_min > 0,  xi ~ N(0, sigma_a^2);   Y_i = g t_i + eta_i, eta_i proper complex, std sigma_i
good event E: |eta_i| <= b_i and |xi| <= b_a,  b_i = sigma_i sqrt(log 1000),  b_a = sigma_a Phi^{-1}(0.9995)
||Y_i| - a h_i| <= b_i                                          (reverse triangle inequality)
=> | |Y_i|/Z - h_i(eps_i) | <= (b_i + h_i b_a)/|Z| <= (b_i + H_i b_a)/(a_min - b_a)
=> |eps_hat_i - eps_i| <= bound_i = (b_i + H_i b_a) / ( (a_min - b_a) m_i )
estimator: eps_hat_i = h_i^{-1}( clip_{h_i(I_i)}( |Y_i| / Z ) ); clipping is harmless since h_i(eps_i) in h_i(I_i).
Union bound: P(E^c) <= 2 exp(-b_i^2/sigma_i^2) + 2(1 - Phi(b_a/sigma_a)) = 2e-3 + 1e-3 = 3e-3  => coverage >= 0.997.
Constants: a_min = 0.75, sigma_a = 0.001, b_a = 3.2905267314919e-3,
m_i = (4.3261945994e-4, 1.3086616724e-4) supplied and (4.2929719577e-4, 1.2943650641e-4) mine,
H_i = supplied interval q upper bounds (2.7198118011157357e-3, 1.3096911596704335e-3).
```
