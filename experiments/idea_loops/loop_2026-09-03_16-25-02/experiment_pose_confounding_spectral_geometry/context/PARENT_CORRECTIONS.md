# Parent Codex correction gate

Date: 2026-09-03

The first Family-1 implementation task invented a numerical self-cell rule that
was not present in the source material and omitted the lower endpoint of the
Hankel antiderivative. For

\[
g(r,r')=\frac{i}{4}H_0^{(1)}(k_b\lVert r-r'\rVert),\qquad
a=h/\sqrt{\pi},
\]

the equal-area disk integral is

\[
\begin{aligned}
I_{\rm self}
&=\frac{i}{4}2\pi\int_0^a rH_0^{(1)}(k_br)\,dr\\
&=\frac{i\pi}{2}
  \left[\frac{r}{k_b}H_1^{(1)}(k_br)\right]_{0}^{a}\\
&=\frac{i\pi a}{2k_b}H_1^{(1)}(k_ba)-\frac{1}{k_b^2},
\end{aligned}
\]

because
$\lim_{r\downarrow0}rH_1^{(1)}(k_br)/k_b=-2i/(\pi k_b^2)$.
Consequently,

\[
[G_D]_{nn}=k_b^2I_{\rm self}
=\frac{i\pi k_ba}{2}H_1^{(1)}(k_ba)-1.
\]

The omitted $-1/k_b^2$ term is not cosmetic: without it the alleged cell
integral tends to $1/k_b^2$ rather than zero as $h\to0$.

The original implementation task also stated the source derivative with the
wrong sign. The admissible identity is

\[
\nabla_s g(z,s)=+\frac{ik_b}{4}H_1^{(1)}(k_bR)\frac{z-s}{R}
=-\nabla_zg(z,s).
\]

The worker detected and repaired that sign through the pose finite-difference
test. Both corrections must be present before any Family-1 result is admitted.
The corrected run must add a direct numerical quadrature check for
$I_{\rm self}$ and a grid-refinement diagnostic; centered derivative agreement
alone validates only internal code consistency.
