# Finite amplitude-window supplement and independent formula check

This independently checks the finite-secant formula seen in the existing branch report at commit `d5e4c2eb66a9a618fecb0928c4b6eaed70b4b130`. Its recovery arrays are not merged with this execution: its printed medians differ. This supplement is not a claim that the formula is absent from the literature. It applies to a **separate bounded amplitude experiment**, not a necessary limit on the full coherent experiment.

## Exact criterion, including its necessity assumptions

Let $W_i=a h_i(u_i)+e_i$, $Z=a+e_a$, with independent Cartesian allowed error intervals $|e_i|\le B_i$, $|e_a|\le B_a$, gain $a\in[a_-,a_+]$, and material box $u_i\in[L_i,U_i]$. Each $h_i$ is positive, strictly increasing and concave on its interval, as certified for class M. Suppose $2\delta<U_i-L_i$ and

$$h_i(U_i)/h_i(L_i)\ge a_+/a_-$$

for both channels. This last condition ensures that the other material can compensate any allowed gain ratio; omitting it would invalidate the claimed necessity for the joint two-channel experiment. Put $d_a=\min(2B_a,a_+-a_-)$. A uniform estimator with both coordinate errors at most $\delta$ exists exactly when

$$a_-h_i(U_i)-(a_-+d_a)h_i(U_i-2\delta)\ \ge\ 2B_i,\qquad i=1,2.$$

**Proof.** Two gain references can agree exactly when $|a-a'|\le d_a$. For material values $u>v$ with gap $s$, the smallest positive amplitude separation after allowable gain compensation is the positive part of $a_-h_i(u)-(a_-+d_a)h_i(v)$. For fixed $s$, concavity makes this expression decrease as $v$ increases, so its infimum is at $u=U_i,v=U_i-s$. It strictly increases with $s$. Therefore the boundary value at $s=2\delta$ bounds every strictly task-bad pair, and the displayed non-strict inequality prevents error-interval overlap for every gap greater than $2\delta$. The feasible-set coordinate midpoint proves sufficiency.

If the displayed inequality fails, continuity gives a gap strictly greater than $2\delta$ whose smallest achievable absolute amplitude separation is below $2B_i$. Use gain $a_-$ for the larger material and choose the other gain in $[a_-,a_-+d_a]$: the upper endpoint applies when the compensated difference stays positive; otherwise an intermediate gain makes the difference zero. This intermediate choice is essential when the endpoint compensation overshoots. The ratio assumption lets the other material be chosen so its two noiseless amplitudes coincide. The reference intervals also overlap. These two admissible worlds share an observation but have disjoint $\delta$-success intervals in coordinate $i$. This proves necessity **for this amplitude experiment**. Additional phase data may distinguish them, so this argument is not a full-complex impossibility proof.

## Covered radial-window example: genuine splitting, no internal optimum

For the electric dipole's ideal tangential modal projection at $z=kR$, the radial factor is

$$c(z)^2=|D_{h_1}(z)|^2=z^{-2}-z^{-4}+z^{-6},$$
$$[c(z)^2]'=-2[(z^2-1)^2+2]/z^7<0.$$

Fix raw projected absolute error $b=10^{-6}$, amplitude-reference error $B_a=0.003291$, $a_-=3/4$, $a_+=5/4$ and $\delta=0.1$. Normalization produces $B_i(z)=b/c(z)$; it does NOT hold normalized noise fixed. The exact criterion above, certified endpoint Mie amplitudes and rational bisection yield:

| Task | Admissible common standoff, ideal amplitude experiment |
|---|---|
| First material | $0.2<kR\le z_1$, $z_1\in[26.73931028381939,26.73931028382011]$ |
| Second material | $0.2<kR\le z_2$, $z_2\in[6.56834737102779,6.56834737102851]$ |

The lower standoff excludes the first sphere's interior. The two certified task windows differ although both materials are isotropic. The joint window is their intersection. No multiple disconnected radial peaks or finite interior maximum occurs in this example. It is an ideal complete-modal-projection example, **not** a certification of the 18-direction receiver over all radii. The latter is separately established only at $kR=2$ in `MANUSCRIPT.md`.

## The inherited small bias budget is not fundamental

The allowed total normalized amplitude errors from the exact finite criterion are at least $2.6758002774787\times10^{-5}$ and $6.6439925866700\times10^{-6}$. Using upper-rounded inherited complex sigmas $2.666408132\times10^{-6}$ and $1.804798621\times10^{-6}$ and radius $\sigma\sqrt{6.908}$ gives each complex tail probability below $0.001$. The real reference with sigma $0.001$ has tail probability below $0.001$ at $0.003291$; `finite_windows.py` proves this with an interval integral of the Gaussian and rational pi brackets.

The second channel therefore permits model bias at least $1.9004269431\times10^{-6}$, over three times the inherited approximately $5.78\times10^{-7}$ sufficient margin. This comparison has **no reference drift allowance** beyond the stated total reference error; drift consumes that margin. The three-event union guarantee is at least $0.997$. This is a less conservative finite certificate, not evidence that a laboratory already meets it, and not a new acquisition algorithm.
