# Finite-error calibration windows for coherent sphere material recovery

**Replacement technical sections, not a submission-ready paper.** Base: immutable A5 and the parent-corrected Q5 V1 reports. No SOM-specific contribution or acquisition-policy advantage is claimed. The two-cycle selector remains withdrawn. Theorems below have different domains; the separately interrogated lossless-sphere certificate must not be transferred to simultaneous lossy spheres.

## 1. Experiments, constitutive domains, and error conventions

Use exp(-i omega t), outgoing spherical Hankel functions of the first kind, nonmagnetic isotropic media, and relative permittivity. In class M, the two spheres are interrogated **separately**. Their exact size parameters are x1=1/5 and x2=3/20, their real permittivities belong to [3/2,4] and [2,5], and a common complex electronics factor satisfies 3/4 <= |g| <= 5/4. Sphere geometry and illumination are known. The retained complex electric-dipole outputs are

$$Y_i=g t_i(\epsilon_i)+e_i+\eta_i,\qquad Z=|g|+e_a+\xi.$$

Here eta_i is proper complex normal with E|eta_i|^2=sigma_i^2; xi is real normal with variance sigma_a^2. These are distinct conventions. The modulus of a noisy complex reference is NOT a real Gaussian observation. A separate IQ-reference implementation may instead bound ||g+xi|-|g|| <= |xi|, using the complex-noise event directly; it uses a different raw measurement budget.

Class C has two simultaneous spheres, centers (-.06,0,0) and (.055,.02,0) m, radii .035 and .025 m, k=18 m^-1, and permittivities epsilon1+i*.03 and epsilon2+i*.05 in the same real intervals. The common receiver-frame x displacement lies in [-2,2] mm. Twelve three-component receivers occupy a radius-.6 m Fibonacci sphere. Four plane-wave direction/polarization pairs are (z,x),(z,y),(x,y),(x,z). The gain annulus is unchanged. The noisy complex reference has sigma_r=.01. A single scenario has one fixed sigma; comparisons and competing worlds must use that same covariance. Real whitening is sqrt(2) times realification divided by complex sigma.

No uncertainty in shape, anisotropy, or SLAM trajectory is included. Structural discrepancy is bounded only when an independent valid bound is supplied, not by fitting or by comparing solvers.

## 2. Attribution and finite calibration information windows

### Proposition 1: exact finite-error criterion in a declared intervention class

Let Theta be compact, T:Theta -> R^q continuous, and let Y_d(theta) be the compact set of observations allowed at design d, including specified model and noise errors. Assume a closed observation graph and measurable extrema. Estimators may output arbitrary vectors in a rectangular task space; requiring the estimate itself to be a jointly feasible nonlinear parameter would be a different problem. Define

$$\Omega_{j,d}=\sup\{|T_j(\theta)-T_j(\theta')|:\mathcal Y_d(\theta)\cap\mathcal Y_d(\theta')\ne\varnothing\}.$$

A uniform coordinatewise deterministic guarantee |T_j-hat T_j| <= delta_j exists if and only if Omega_j,d <= 2 delta_j for every coordinate.

**Proof.** At a shared observation, one output must lie within delta_j of both task values, which proves necessity. Conversely, form the feasible parameter fiber at y, and return the midpoint between the infimum and supremum of each T_j over that fiber. Every true task value is at most half that coordinate range away. Compactness gives extrema; the stated graph assumptions give a measurable estimator. Empty fibers are rejected. This is a set-membership result, not a new Fisher-information theorem.

For full Euclidean error balls with common whitened radius r, two observation sets intersect exactly when ||mu(theta)-mu(theta')|| <= 2r. A covered separation greater than 2r for every task-bad pair is sufficient. An infimum over the strict set |T_j-T_j'|>2delta_j need not attain its threshold; do not silently replace the exact intersection condition by an inequivalent strict-infimum assertion.

Define the task-specific window

$$\mathcal W_j(\delta_j)=\{d:\Omega_{j,d}\le 2\delta_j\},\quad \mathcal W_{\rm joint}=\bigcap_j\mathcal W_j.$$

A design includes standoff, aperture, polarization, frequency, exposure/repeats, receiver noise, reference drift, and forward-model fidelity. **Window splitting** means these task/mode-specific admissible design sets differ. It does not assert disconnected radial intervals, two peaks, or anisotropy. Such shapes require a separate demonstration. A4's fixed-relative-SNR kR=1 radial geometry optimum is not a material optimum.

### Proposition 2: when a reference can repair a finite ambiguity

Restrict the permitted intervention to a specified independent channel R_d(theta)+B_ref, with a full error ball B_ref of radius r_ref and unchanged main-data experiment. It repairs the target guarantee exactly when, for every task-bad pair with overlapping main-data sets,

$$\|R_d(\theta)-R_d(\theta')\|>2r_{\rm ref}.$$

**Proof.** The augmented observation sets are Cartesian products. Their intersection is nonempty precisely when both component intersections are nonempty. Apply Proposition 1. With state-dependent or correlated error sets, use the actual joint set, not this product simplification.

This necessity-and-sufficiency statement belongs to the declared error-set and intervention class. It does not establish that one amplitude reference is universally minimal among all electromagnetic acquisitions, additional polarizations, frequencies, or priors. Minimizing the number/cost of references additionally requires covering all task-bad pairs by admissible interventions. That optimization has not been certified for class C.

Correctly stacked designs share the SAME material, geometry, and electronics parameters. Their conflict relation is the intersection of the per-design relations. Giving every block an unrelated freely fitted gain changes the statistical experiment and can destroy the benefit of stacking. Adding observations cannot enlarge the feasible fiber when the original data, errors, and nuisance model are retained; reallocating a fixed power/time budget is not covered by that monotonicity.

### Local diagnostic, not a substitute for Propositions 1--2

For real-whitened data write J=[A B C], with actual Maxwell material, geometry, and electronics derivatives. For a scalar common gain, C consists of realifications of F and iF. For receiver translation B must contain g times the spatial derivative of the Maxwell observation field, not an arbitrary nuisance column. The material diagnostic is Q_[B,C] A; geometry uses Q_[A,C] B. For one material component the other material columns are nuisance too.

These are quotient directions, not three unique orthogonal physical compartments: parameter effects can overlap. A discrepancy tangent to their span can bias parameters while leaving little residual. A reference adds electronics information under a correct model, but changing the projection can increase or decrease a misspecified estimator's bias. A generic local expression is

$$\delta m_{\rm bias}=(A^TQ_NA)^{-1}A^TQ_N b,$$

on its identifiable support. It neither certifies a basin nor a global material error. The finite windows, not this matrix alone, determine the acceptance claim.

## 3. Closing the restricted-modal interval proof

Let a be the outgoing electric Mie coefficient, t=-a, and q=i a/(1-a). For lossless spheres in class M, q is positive real and |t|=h=q/sqrt(1+q^2). Introduce entire functions, with s=epsilon*x^2 and u=x^2,

$$J(s)=\sum_{n\ge0}\frac{(-1)^n(2n+2)s^n}{(2n+3)!},\quad D(s)=\sum_{n\ge0}\frac{(-1)^n(2n+2)^2s^n}{(2n+3)!},$$
$$C(u)=\sum_{n\ge0}\frac{(-1)^n u^n}{(2n)!},\quad S(u)=\sum_{n\ge0}\frac{(-1)^n u^n}{(2n+1)!},$$
$$Y=-C-uS,\quad Z=(1-u)C+uS,$$
$$n=\epsilon J(s)D(u)-J(u)D(s),\quad d=\epsilon J(s)Z(u)-Y(u)D(s),\quad q=x^3 n/d.$$

Here j1(sqrt(s))=sqrt(s)J(s), D_j(sqrt(s))=D(s), y1(x)=Y/x^2, and D_y(x)=Z/x^3. This formulation avoids cancellation and uncertified square-root inputs. The original Mie real denominator equals d/x^2.

### Theorem 3: covered monotonicity and justified endpoint bounds

On all of class M, q>0, q'>0, q''<0. Consequently h'>0 and h''<0. Safe rounded constants are

| Sphere | m <= min h' | H >= max h |
|---|---:|---:|
| 1 | 0.00045514789 | 0.00268777635 |
| 2 | 0.00013998065 | 0.00129313278 |

**Computer-assisted proof.** `interval_certificate.py` tiles the two intervals with 2500 and 3000 closed rational boxes of exact width 1/1000. Adjacent endpoints are checked equal and the last endpoint equals the domain endpoint. Every arithmetic operation rounds outward on the 10^-36 integer lattice. The sizes are rational 1/5 and 3/20 throughout, never binary approximations. The first 12 series terms and their first two derivatives are enclosed. For n>=12 and derivative order r<=2, the successive absolute term ratios, before multiplication by s<=1, are

$$\frac{n+1}{(n+1-r)(2n+2)(2n+1)}\ (C),\quad
\frac{n+1}{(n+1-r)(2n+3)(2n+2)}\ (S),$$
$$\frac1{2(n+1-r)(2n+5)}\ (J),\quad
\frac{n+2}{2(n+1)(n+1-r)(2n+5)}\ (D).$$

Each is below 1/2, so twice the first omitted absolute term is a valid symmetric remainder. Rational differentiation gives q',q''. The complete box ledger establishes positive denominator, positive q,q', and negative upper q''. In particular the largest q'' upper bounds are -0.00014921660... and -0.00003966935.... Since

$$h''=q''(1+q^2)^{-3/2}-3q(q')^2(1+q^2)^{-5/2}<0,$$

the minimum h' occurs at the upper permittivity endpoint. The exported lower bound q'/(1+q^2)^2 is smaller than h' and requires no uncertified root. The upper h bound follows from h<=q. Both inherited conservative slope/cap pairs are thereby independently certified. This is not 5500 independent theorems or a proof-assistant verification.

The corrected Bessel identities, independently regression-tested, are

$$D_j'= -j'/z-j+j/z^2,\qquad D_j''=-j''/z-j'+2j'/z^2-2j/z^3.$$

Any midpoint-Lipschitz alternative must use the complete halfwidth. No midpoint bound is needed in this certificate. The fixed-size theorem does not silently include machining, frequency, or position tolerance; those enter a separate validated error budget.

## 4. A finite Maxwell modal readout and a nonempty conditional budget

The next result is a constructive sufficient readout, not a claim that spherical harmonic measurement or Mie inversion is new.

### Theorem 4: eighteen-direction readout with an all-order leakage bound

Illuminate a single known sphere by a unit plane wave along x, polarized along z. Measure the two tangential complex field components on its centered sphere kR=2. Use the three Gauss--Legendre z nodes and six equally spaced azimuths (any common azimuth offset). Let P_lm=grad_S Y_lm/sqrt(l(l+1)), and use M=z_l C_lm, N=-(sqrt(l(l+1))z_l Y_lm*rhat/(kr)+D_z P_lm). The incident N_10 coefficient is c=-sqrt(6pi). Define

$$\widehat t=-\frac{\sum_j w_j P_{10}(n_j)^*\cdot E_s(Rn_j)}{c D_{h_1}(2)}.$$

For 0<x<=1/5 and 3/2<=epsilon<=5, the deterministic quadrature leakage obeys

$$|\widehat t-t|\le \frac{20}{9}\frac{100}{99}\,T_5(x),$$
$$T_l(x)=\frac{2x^{2l+1}}{2^l(2l-1)!!}\left(8+10l+\frac{40}{2l-1}\right).$$

The respective bounds at x=.2 and .15 are below 1.899e-10 and 8.019e-12, before multiplication by g.

**Proof of low-order exactness.** Tangential projection makes P_10 dot P_lm a spherical polynomial of degree at most l+1, and P_10 dot C_lm one of degree at most l. The 3-by-6 quadrature integrates spherical polynomials through degree five exactly. All modes l<=4 are therefore orthogonal to the target except N_10. Finite implementation errors in irrational nodes and receiver positions are NOT covered by ideal quadrature exactness; they belong to the readout-error budget.

**Proof of the all-order tail.** For l>=1, x<=.2, and sqrt(epsilon)x<=.5, the alternating spherical-Bessel series imply

$$0.975\frac{z^l}{(2l+1)!!}\le j_l(z)\le\frac{z^l}{(2l+1)!!},$$
$$0.95\frac{(l+1)z^{l-1}}{(2l+1)!!}\le D_{j_l}(z)\le\frac{(l+1)z^{l-1}}{(2l+1)!!}.$$

The finite outgoing-Hankel polynomial and its reversed-coefficient bound give -y_l(x)>=.7(2l-1)!!/x^(l+1), |h_l(x)|<=(5/4)(2l-1)!!/x^(l+1). The recurrence D_y=y_(l-1)-l*y_l/x then gives D_y>=.65*l*(2l-1)!!/x^(l+2)>0. The electric and magnetic real-denominator ratios consequently satisfy

$$|a_l|\le\frac{10x^{2l+1}}{(2l-1)!!(2l+1)!!},\qquad
|b_l|\le\frac{4x^{2l+1}}{(2l-1)!!(2l+1)!!}.$$

For clarity, a Hankel reversed coefficient r positions from its leading term has modulus at most x^r/r!. Thus -y_l >= [cos(x)-(exp(x)-1)] times the leading magnitude; cos(x)>=.98 and exp(x)<=1.25 justify .7. At the readout radius,

$$|h_l(2)|\le\frac{8(2l-1)!!}{2^{l+1}},\quad
|D_{h_l}(2)|\le\frac{8(2l-1)!!}{2^{l+2}}\left(l+\frac4{2l-1}\right),$$

using exp(2)<8 and the same recurrence. The incident coefficients have squared norm 2pi(2l+1) in each electric/magnetic family. The addition theorem gives sum_m||P_lm||^2=sum_m||C_lm||^2=(2l+1)/(4pi). Cauchy--Schwarz, with 1/sqrt(2)<=1, bounds the tangential l-th field by T_l. For l>=2, T_(l+1)/T_l <= x^2(l+1)/(2l(2l+1)) < 1/100; sum the geometric majorant. Finally |D_h1(2)|^2=13/64 and |D_h1(2)|>=9/20. Positive weights and 4pi*sup||P10||/|c|=1 prove the asserted alias bound. The proof includes every l>=5; a computed L-to-L+1 difference is not used.

### Instrument sensitivity and noise

If each measured tangential vector has a relative complex operator error with norm <=eta, its contribution to the retained modal error is <=g_max*eta*C_i, where the conservative explicit bound is

$$C_i=3H_i+4B_i+\frac{20}{9}\frac{100}{99}T_2(x_i).$$

B_i bounds the magnetic dipole; `readout_certificate.py` encloses its rational Mie ratio over the same complete material boxes. It gives B1<=2.146777e-5 and B2<=6.788701e-6. The resulting C1<=.013097449 and C2<=.005080797 include all higher modes, not just the target amplitude. This relative-operator assumption is an instrument acceptance specification; it cannot absorb an arbitrary geometry tangent without a separate bound.

For independent, equal-variance proper complex noise on the two orthonormal tangential components, the squared norm of the projection row is exactly

$$\|w_{\rm read}\|^2=\frac{224}{1053}.$$

Indeed the GL3 weights are (5/9,8/9,5/9), the polar squared sines are (2/5,1,2/5), and the six azimuth weights give sum_j w_j^2 sin^2(theta_j)=56pi^2/81. Combining this with |c|^2=6pi and |D_h1(2)|^2=13/64 yields the result. Per-component complex standard deviation <=2*sigma_i is therefore sufficient for modal deviation <=sigma_i.

### Theorem 5: uniform finite material recovery with the declared readout budget

If |e_i|<=b_model,i, |e_a|<=b_ref, and on a noise event |eta_i|<=b_noise,i, |xi|<=b_refnoise, put B_a=b_ref+b_refnoise<a_min. Clip |Y_i|/Z to the known h_i interval and invert h_i. Then

$$|\widehat\epsilon_i-\epsilon_i|\le
\frac{b_{\rm model,i}+b_{\rm noise,i}+H_iB_a}{(a_{\min}-B_a)m_i}.$$

**Proof.** Reverse triangle inequality bounds ||Y_i|-a*h_i| by the complex output error. The ratio denominator is at least a_min-B_a. Subtract h_i, bound the reference term by H_i B_a, and use clipping nonexpansiveness followed by the global mean-value lower slope m_i. This is finite and uniform on the full material intervals; local full rank is not its justification.

A target delta is certified whenever

$$b_{\rm model,i}+(H_i+\delta m_i)b_{\rm ref}
\le \delta(a_{\min}-b_{\rm refnoise})m_i-b_{\rm noise,i}-H_i b_{\rm refnoise}.$$

This inequality is sufficient for this estimator/budget. Its failure is NOT necessary failure of material identification or of another estimator.

The frozen readout realization uses 16 independent repeats, eta=.0005, absolute reference bias/drift <=.00025, and other independently bounded projected errors <=2e-6 per sphere. The latter must include finite-node positioning, incident-field and background subtraction, and size/frequency uncertainty not already accounted for. It is **not** an experimentally established bound. Put

$$b_{\rm model,i}=1.25\,[b_{\rm alias,i}+\eta C_i]+2\times10^{-6}.$$

The single-repeat modal sigmas remain 2.6664081e-6 and 1.8047986e-6; the real-reference sigma is .001. Average 16 times, and choose b_noise,i=3*sigma_i/4, b_refnoise=3.5*.001/4. Common bias is not divided by four. The union failure probability is bounded by

$$2e^{-9}+\frac{\sqrt{2/\pi}}{3.5}e^{-49/8}
<\frac{191}{252000}<0.001.$$

The last rational bound follows from exp(9)>8000, exp(49/8)>450, and sqrt(2/pi)<.8. It does not require independence between the three error events; independent repeats and specified per-component covariance are needed for variance reduction.

| Quantity | Sphere 1 | Sphere 2 |
|---|---:|---:|
| Model/instrument bias bound | 1.018615e-5 | 5.175508e-6 |
| Averaged noise-event radius | 1.999807e-6 | 1.353599e-6 |
| Total reference error B_a | .001125 | .001125 |
| Uniform material-error upper bound | .0446230 | .0761618 |
| Remaining absolute output margin at delta=.1 | 1.887519e-5 | 2.498920e-6 |

Thus the **mathematical readout/specification class is nonempty**. Its raw budget is 18 directions times two tangential complex components times 16 repeats times two separately read spheres: **1152 complex field samples plus 16 real reference samples**. It does not establish a laboratory's ability to meet the .05% relative-channel and 2e-6 other-error specifications. Correlated drift, probe loading, mutual coupling, imperfect plane waves, background subtraction, and channel covariance require independent acceptance measurements. A practical IQ reference uses 16 complex, not 16 real, reference samples and must use its own noise event. This explicitly paid repetition budget does not retroactively change the inherited approximately 5.78e-7 single-read certificate.

## 5. A Maxwell-specific route to a valid model-error radius

A solver discrepancy is not an error bound. The following residual estimate supplies a possible certifier and explains exactly what remains unimplemented for class C.

Let D be the union of the two known balls, chi=epsilon_r-1, P=chi E, and

$$\mathcal A_\chi P=(\chi^{-1}-\mathcal T)P=E^{\rm inc},\quad
\mathcal T P=k^2\int_D G_k(\cdot,r)P(r)\,dr.$$

Use the standard outgoing Maxwell volume operator, including its real principal-value/depolarization term. Its radiation identity gives Im<P,T P>>=0 for supported currents. This sign assumes the stated exp(-i omega t) convention.

### Proposition 6: passive residual enclosure

For positive loss and alpha=ess inf_D Im(chi)/|chi|^2>0, and any square-integrable approximate polarization P_h for which the continuous residual is defined,

$$r_h=E^{\rm inc}-\mathcal A_\chi P_h,\qquad
\|P-P_h\|_{L^2(D)}\le\alpha^{-1}\|r_h\|_{L^2(D)}.$$

**Proof.** For u=P-P_h, Im<u,A_chi u> <=-alpha||u||^2 by passivity. The equation A_chi u=r_h and Cauchy--Schwarz yield alpha||u||^2<=||u||||r_h||. No small-Born or weak-coupling assumption is used. This standard coercivity/energy argument is not claimed as an original optical theorem; compare Miller et al. [L3].

For receiver operator S,

$$\|gS(P-P_h)\|\le g_{\max}\|S\|\alpha^{-1}\|r_h\|.$$

For N three-component point receivers at least d_min from D,

$$\|S\|^2\le\frac{N|D|}{16\pi^2d_{\min}^2}
\left(2k^4+\frac{2k^2}{d_{\min}^2}+\frac6{d_{\min}^4}\right).$$

This follows by integrating the squared Frobenius norm of k^2 G. Its transverse and longitudinal eigenvalues give exactly the expression in parentheses. For several illuminations, stack residual L2 norms in quadrature; do not multiply the per-illumination bound by the number of illuminations twice.

In the declared class C, alpha >= .05/(16+.0025)=.00312451179..., d_min>=.503 m and |D|=4pi(.035^3+.025^3)/3. A conservative output amplification is g_max||S||/alpha<1583. These are analytic upper/lower constructions, while the printed decimals in `finite_pair.json` are ordinary evaluations, not outward-certified constants. A safe rational policy can use alpha>=.0031245, ||S||<=3.956 and multiplier<=1583.

To extend a validated polarization at a box center to a material box B, use the residual identity

$$r_\chi=r_{\chi_0}+(\chi_0^{-1}-\chi^{-1})P_h.$$

Then the supremum norm of the second term over B is bounded directly, before applying Proposition 6. Receiver-shift uncertainty acts on S, not on illumination in the declared class; its additional radius follows from a validated supremum of the spatial Green derivative. These bounds can feed a finite box-cover algorithm.

**Unresolved implementation:** there is no certified continuous residual L2 enclosure for the present multipole/DDA solutions, no certified spatial-derivative enclosure over every parameter box, and no complete class-C pair cover. In particular neither the DDA grid differences nor multipole-order differences may be inserted for ||r_h||. The current simulator checks are valuable diagnostics but do not execute Proposition 6 as a certificate.

## 6. Independent nonlinear checks and finite competitors

The V2 stream was registered before execution, using seeds 2026091107 and 2026091108. It contains 12 new scenes, not reused V1 final evidence. Generation uses a 312-cell vector DDA at spacing .011 m, fill quadrature 4. Inversion uses a separate standard vector-spherical-wave cluster solver, L=3, GL16-by-32 translation quadrature. Three fixed starts and 80 objective evaluations per start were used. The common gain is exactly profiled over its annulus; no log determinant was added.

| Frozen comparison | Material successes, both errors <=.1 |
|---|---:|
| No extra reference | 2/12 |
| One noisy complex reference, joint GLS | 7/12 |
| One fixed extra EM complex scalar | 2/12 |
| One random extra EM complex scalar | 2/12 |

The extra scalar count is matched, not the hardware cost. Fixed/random arms share the same candidate noise pool. These counts do not establish an adaptive-method advantage or a population success probability.

Posthoc, explicitly labelled **oracle diagnostics**, not additional fair baselines, give: matched-model/no-reference 11/12; matched-model/reference 12/12; independent data with true geometry and noisy reference 7/12; true complex gain 7/12; true gain plus true geometry 7/12; noise-free independent data with both oracles 8/12. Exact gain-amplitude-only anchoring also gives 7/12. Therefore the gain-amplitude scale is important, but calibration error alone cannot explain the remaining failures. These local-optimizer diagnostics do not exclude optimization error as an additional cause; they do not prove that every failure is solely due to structural discrepancy. The unchanged geometry-oracle count does not prove geometry never matters. Between .605 and .941 of the model-discrepancy norm projects into the local joint parameter tangent in these scenes; this is a bias diagnostic, not a validated discrepancy bound.

Across the first four scenes, L3-to-L4 relative differences are about 2.25e-6--2.44e-6, L4-to-L5 about 4.78e-8--9.08e-8, and quadrature differences about 2.7e-15--4.1e-15. DDA-to-L5 differences are .64%--1.52% for the original grid and .63%--1.67% for the finer grid. The finer-grid error does not monotonically decrease. Boundary continuity and modal quadrature are separately regression-tested. None of these comparisons is a continuum error bound.

An exploratory finite competitor fixes world A at (epsilon1,epsilon2,shift_mm)=(3.8,4.8,.4), g_A=exp(.2i), and world B's second material at 4.59. With a shared complex sigma=6.566412640846138e-5, optimizing the remaining admissible parameters gives the following numerical candidates:

| Candidate | World-B epsilon1 | World-B shift (mm) | Real-whitened distance, L7 | Gaussian binary risk evaluated at numerical means |
|---|---:|---:|---:|---:|
| No reference | 3.67332848 | .38975161 | 3.39773051 | .04467229 |
| sigma_r=.01 reference | 3.68182273 | .38016996 | 5.20127714 | .00465252 |

Their .1 material-success boxes are disjoint because epsilon2 differs by .21. The L5/L7 changes are small, but the optimization values remain **upper-distance numerical candidates**, not covered lower bounds and not certified Maxwell counterexamples.

For a future valid combined real-whitened mean-error radius beta, the corresponding binary risk is at least Phi(-(d_num+beta)/2). Thus the referenced pair would disprove a uniform .997 success guarantee if a continuous-model certificate established beta < -2 Phi^-1(.003)-5.20127714 (about .29). This is a concrete certification target, not an already proved impossibility. It does not show failure of lower-confidence recovery, more accurate references, or other acquisitions.

## 7. What changes relative to Born under the SAME prior

For known support and geometry, the first Born field has the form F_B=sum_j B_j(chi_j). With the fixed losses ell=(.03,.05), restrict chi_j=ell_j(q+i). Then F_B=(q+i)sum_j ell_j B_j, and an unknown common gain produces an exact ambiguity between q=30 and q'=40, g=1 and g'=(30+i)/(40+i). Both gains lie in the annulus and the real permittivities are (1.9,2.5) and (2.2,3.0). The already corrected q'=45 example must not be used without rescaling.

The known loss is itself constitutive prior information. Away from that proportional subfamily, fixed nonzero losses can provide a scale/phase anchor even in Born. It is false that all Born problems with unknown gain are nonidentifiable. Conversely, estimating losses as free unknowns changes the ambiguity class and cannot be compared as the same-prior experiment.

Distinguish five mechanisms. (1) **Internal response:** exact single-sphere Mie coefficients depend nonlinearly on material, even without any other sphere. (2) **Relative electric/magnetic response:** a_l/b_l can provide gain-invariant structure only when the channels are jointly calibrated and strong enough. (3) **Higher spatial modes:** l>=2 can add angular structure but may lie below the noise/model floor. (4) **Object diversity:** two separately read known spheres can have different nonlinear material response. (5) **Inter-object rescattering:** the cluster equation c=T(b+Uc) modifies illumination through other objects. Its effect is not synonymous with any of the preceding mechanisms, and is not guaranteed to improve finite recovery.


A mechanism ablation must hold losses, electronics, sensor geometry, physical noise scale, and reconstruction tolerance fixed. This round implements the full cluster model and the separate class-M readout; it has NOT established a controlled causal recovery advantage for each of the five mechanisms.

## 8. Antenna and receiver design implications

Use a task window, not maximal raw amplitude, to select standoff and polarization. A design favorable for spatial derivatives may remain poor at separating material scale from gain. A4's radial window, class-M modal-readout window, and class-C material window are different questions. At each candidate design, first state aperture/power/time/noise conventions, then remove only the permitted nuisance, then test finite task-bad pairs with valid error sets. Report a failed certificate as unresolved unless an explicit admissible competitor is certified.

The eighteen-direction construction is a concrete antenna-readout specification: two calibrated tangential polarizations, a centered scan kR=2, shared coherent gain, and controlled inter-channel errors. kR=2 is a constructive point, **not a proved optimal standoff**. A physical electric-field probe requires a calibrated probe-response operator. Its loading, polarization leakage and electronics drift must fit the stated error set. An electronics reference before the antenna does not by itself calibrate phase-center motion, illumination mismatch, or antenna--object coupling.

A known invertible digital gain A does not improve discrimination: transforming both means and covariance preserves their Mahalanobis distance. For an analog amplitude gain G before a later noise stage,

$$D^2(G)=\Delta s^*\big(\Sigma_{\rm pre}+G^{-2}\Sigma_{\rm post}\big)^{-1}\Delta s.$$

It increases with G when the stage noises are positive semidefinite and independent in the stated model, and saturates at the pre-amplifier-noise limit when that covariance is positive definite. It cannot remove an exact material/electronics ambiguity. This is ordinary cascade-noise analysis, not a newly identified theory called "Amplify" or "Gain-F"; those names remain unspecified in the request. A noise figure, if used, must be defined at a reference source temperature and bandwidth, not confused with the unknown complex calibration gain.

## 9. Implementable algorithms and scientific gates

### Annulus-constrained noisy-reference GLS (implemented)

```
v = concatenate(F(theta)/sigma, 1/sigma_r)
b = concatenate(y/sigma, z/sigma_r)
g0 = inner_product(v,b)/inner_product(v,v)
g  = clip(abs(g0), .75, 1.25) * exp(i*arg(g0))
r  = sqrt(2) * concatenate(real(b-v*g), imag(b-v*g))
optimize theta from each frozen start; keep the smallest objective
report BOTH material errors, geometry error, gain error, and all failed starts
```

At g0=0 every phase is optimal; choosing phase zero is valid. A marginalization logdet changes the objective. Annulus profiling is not always equivalent to unrestricted GLS formulas.

### Certified class-M recovery (proof above, quadrature and bounds tested)

```
check domain, probe calibration, covariance, bias and drift acceptance records
acquire both tangential components at 18 directions, for 16 repeats
apply the known linear modal row; average COMPLEX samples coherently
average real amplitude references (or use a separately specified IQ error event)
compute B_model, B_noise, B_reference and the uniform task bounds
if any bound exceeds .1: return UNRESOLVED, not a material guarantee
clip abs(Y_i)/Z into [h_i(eps_min), h_i(eps_max)]
invert monotone h_i by bracketing; retain phase data for diagnostics
```

### Class-C covering certifier (specification; Maxwell enclosure oracle missing)

```
input compact parameter domain, target tolerances, permitted error sets
request certified continuous forward enclosures over a box partition
if a model/roundoff enclosure is absent: return UNRESOLVED
for each box pair that can contain task-separated parameters:
    bound all possible observation distances from BELOW
    if the lower bound exceeds the error-set overlap radius: discard pair
    else subdivide a box; never use optimization distances as lower bounds
    if resources exhausted: return UNRESOLVED with surviving boxes
if no task-bad pair survives: issue covered finite separation certificate
```

A verified competitor can instead upper-bound a distance using validated point enclosures; it need not solve the pair-distance optimization globally. No such class-C certificate is claimed by the current code.

## 10. Submission boundary

Established: restricted-modal continuum enclosure; a finite Maxwell readout with all-order quadrature-tail control; conditional uniform material recovery at an explicit, nonzero instrument budget; exact finite-window logic; a passive-Maxwell route to a legitimate forward-error bound; and independently generated nonlinear diagnostics.

Not established: the instrument specifications in a real measurement system; a continuous-model enclosure for the simultaneous lossy pair; covered class-C finite separation or its certified finite impossibility/remedy; a new acquisition algorithm's advantage; or novelty sufficient for a TAP method paper. These missing links are substantive, not formatting tasks. The closest-paper comparison in `LITERATURE_AUDIT.md` prevents presenting classical Mie inversion, spherical quadrature, Schur profiling, or passivity as new contributions.


## Supplement: sharper amplitude criterion and a covered splitting example

`WINDOW_PROOF.md` supplies the exact finite-secant criterion in the explicitly restricted amplitude experiment, its necessity construction, a rigorous radial-window split, and a less conservative inherited single-read bias allowance. `finite_windows.py` verifies the endpoint Mie values, reference noise event and radial root brackets with outward integer intervals. This supplement does not extend the 18-direction implementation to every radius or turn amplitude-only indistinguishability into full-complex impossibility. The existing branch report containing the same formula was read and independently checked, not claimed as new literature priority.


### Recorded execution provenance

The verified raw file `results/independent_v2.json` has SHA256 `f30ffa7e185704d165e2012f425d5103c409fe53d87d7968dd8337ae41bc8f1a`. Its material-error medians are (0.2356033,0.4859658), (0.0554202,0.0693530), (0.2308569,0.4769725), and (0.2358079,0.4857988), respectively. Runtime was 11.7522 seconds and recorded process peak RSS was 295184 KiB in the stated environment; these are not cross-machine budgets. A fresh replay produced identical raw observations and fitted parameters. An initial orchestration call timed out after completing main recovery; the remaining diagnostic steps were resumed without changing data or settings. The original working branch acquired another report with different printed medians; that report is neither overwritten nor merged into this raw stream. Submission therefore uses `research/q5-readout-verification-v2`, branching from the immutable modal registration commit.
