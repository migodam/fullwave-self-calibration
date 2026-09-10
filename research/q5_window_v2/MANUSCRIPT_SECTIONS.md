# Q5 V2: replacement mathematical sections, not a submission-ready manuscript

Baseline: 418ef53a7d7d2f79d6228d023baf8b101e3b3e2e. The A5 theory and Q5 V1 parent corrections remain immutable. These sections replace overextended window, modal-readout, model-error, and numerical-validation claims. Their English is intended for eventual manuscript use; the accompanying PR and research report are in Chinese.

## 1. Scope, conventions, and claim hierarchy

We use time dependence $e^{-i\omega t}$, outgoing $h_l^{(1)}$, positive imaginary permittivity for passive materials, and relative permeability one. A proper complex noise variable $n\sim\mathcal{CN}(0,\sigma^2)$ has $\mathbb E|n|^2=\sigma^2$ and real/imaginary variances $\sigma^2/2$. Real whitening is $\sqrt2[\Re n;\Im n]/\sigma$. The latter factor must not be silently changed between testing bounds and least squares.

Three noninterchangeable model classes are used.

**M (resolved modes).** Two separately acquired lossless homogeneous-sphere electric dipole coefficients, with exact size parameters $x_1=1/5$, $x_2=3/20$, materials $\epsilon_1\in[3/2,4]$, $\epsilon_2\in[2,5]$, and a shared complex gain $3/4\le|g|\le5/4$. Shape, geometry, constitutive loss, and channel normalization are known unless explicitly allocated as bounded readout error. An optional E/M variant measures both electric and magnetic dipole coefficients of each sphere. Separately acquiring these modes is a substantive intervention, not an automatic property of class C.

**D (projected dipole tensor).** A single known sphere, $ka=1/10$, and an isolated electric dipole channel. The radial observation direction is known; range $R$ and a common scalar amplitude/phase are unknown. Three independent electric-dipole excitations and vector receiving components permit the transverse/radial tensor readout used below. Other multipoles are excluded by the measurement definition or bounded leakage. The exterior outgoing dipole field is exact in this channel. This is not an assertion that a plane-wave-driven sphere has no other modes.

**C (simultaneous lossy spheres).** Sphere centers $(-.06,0,0)$ and $(.055,.02,0)$ m, radii $.035,.025$ m, $k=18$ m$^{-1}$, materials $[1.5,4]+.03i$ and $[2,5]+.05i$, receiver-frame translation $s\in[-.002,.002]$ m along $x$, and $3/4\le|g|\le5/4$. Twelve three-component receivers lie on a nominal radius $.6$ m. Four illuminations are $(d,p)=(e_z,e_x),(e_z,e_y),(e_x,e_y),(e_x,e_z)$. Both spheres are present, including inter-object rescattering. This is the main unresolved recovery class. No result is extended to unknown shape, anisotropy, or SLAM.

A conditional theorem, an exact rational certificate, a floating-point diagnostic, and an independent numerical recovery are different evidence levels. None entails the others without an explicit bridge.

## 2. Finite calibration-aware information windows

For an acquisition design $a$, let $\theta=(\epsilon,s,g)$ denote a complete admissible world, $\mu_a(\theta)$ the observation mean, and $\mathcal E_a(\theta)$ an explicitly justified error set. Define the observation fiber

$$\mathcal F_a(y)=\{\theta:y\in\mu_a(\theta)+\mathcal E_a(\theta)\}.$$

For a vector-valued task $\tau(\theta)$ and coordinatewise tolerances $\delta_j>0$, define its robust window by the existence of an estimator satisfying

$$|\widehat\tau_j(y)-\tau_j(\theta)|\le\delta_j\quad\text{for every admissible }(\theta,y).$$

### Proposition 1: exact fiber criterion

When estimates may lie in the coordinate box containing the task prior, such an estimator exists if and only if every nonempty fiber obeys

$$\sup_{\theta\in\mathcal F_a(y)}\tau_j(\theta)-\inf_{\theta\in\mathcal F_a(y)}\tau_j(\theta)\le2\delta_j\quad\text{for all }j.$$

**Proof.** Necessity follows because two task values compatible with the same datum must both be within $\delta_j$ of the same estimate. For sufficiency, take the midpoint of each coordinate's supremum and infimum. This is an existence statement, not a claim that a local optimizer returns that midpoint. If an estimate must also belong to a nonrectangular joint feasible set, this coordinatewise construction needs additional analysis.

For common full Euclidean error balls of radius $B$, the exact equivalent pointwise condition is

$$\|\mu_a(\theta)-\mu_a(\theta')\|>2B$$

for every pair for which at least one task-coordinate separation exceeds $2\delta_j$. Full balls intersect exactly when their centers are at distance at most $2B$. A strictly positive **uniform** margin above $2B$ is sufficient. Because the bad-pair set is defined by a strict parameter inequality, replacing the pointwise condition by a strict inequality for its infimum is not in general necessary.

If a model-error ball only encloses the actual discrepancy set, the criterion is sufficient, not necessary for the physical experiment. Gaussian noise is unbounded: robust conclusions apply on a stated probability event, not for all Gaussian realizations.

The generic fiber/separation principle is not an originality claim. Gilev et al. [R1], in particular, already define a parameter-tolerance-dependent minimum forward distance in inverse Mie scattering and distinguish discrete search from continuous guarantees.

### 2.1 What splitting means here

Define $\mathcal W_{\rm mat}$, $\mathcal W_{\rm geom}$ and $\mathcal W_{\rm gain}$ using different tasks on the **same** complete-world model and error budget. Their differences are task-wise information-window splitting. A joint coordinatewise guarantee requires all relevant conditions. This does not mean a single window must have two disconnected distance intervals.

The earlier A4 discussion also used splitting for tensor/polarization branches in anisotropic media. That hypothesis is not proved here and is outside the current scope. The finite example in Section 3 establishes only isotropic, task-wise splitting. No renaming is allowed to turn it into an anisotropic result.

For a finite catalogue of interventions, each acquisition has a confusable full-world-pair set. With independent product error sets and a common complete world across all acquisitions, uniform distinguishability is equivalent to removing every bad pair from the intersection of these sets. Thus a minimum-cost catalogue cover is a conditional necessity-and-sufficiency formulation **within that fixed catalogue**. Profiling a different gain or material independently for each acquisition destroys this equivalence when the physical experiment shares those parameters. This is a design formulation, not a new successful selector; the V1 selector remains withdrawn.

## 3. A covered finite range window and exact material ambiguity

After orthonormal compression of the transverse and radial components of the outgoing electric-dipole tensor, its projective representative is

$$v(z)=\begin{bmatrix}\sqrt2(z^2+iz-1)\\2(1-iz)\end{bmatrix},\qquad z=kR>ka.$$

The suppressed scalar includes the material coefficient, electronics, and the common $e^{iz}/R^3$ propagation factor. It is not suppressed from the absolute noise budget.

### Proposition 2: finite projective chord

For two ranges $z,w>0$,

$$\sin^2\theta(z,w)=\frac{2(z-w)^2[(z+w)^2+z^2w^2]}{(z^4+z^2+3)(w^4+w^2+3)}.$$

**Proof.** The squared norm is $2(z^4+z^2+3)$. Direct expansion gives $\det[v(z),v(w)]=2\sqrt2(z-w)(z+w-izw)$. In complex dimension two, squared sine of the projective angle is squared determinant magnitude divided by the two squared norms.

For $z,w\in[a,b]$ with $|z-w|\ge d$,

$$\sin^2\theta\ge\frac{2d^2(4a^2+a^4)}{(b^4+b^2+3)^2}.\tag{1}$$

This is a covering bound on all parameter pairs, not a grid of projective angles. For an actual signal $cv(z)$, minimizing over every complex competing scalar gives distance $\|cv(z)\|\sin\theta$. Relaxing an annular gain constraint this way can only lower that distance, so it is a valid lower bound for the constrained experiment.

**Certified finite example.** Let the nominal design $z_0$ vary continuously over $[.6,1.5]$, let the true range satisfy $|z-z_0|\le.05$, require error at most $.02$, and bound the total deterministic observation error by $.005 S_{\min}$, where $S_{\min}$ is a valid minimum signal norm over the entire chosen class. Ninety consecutive rational design intervals of width $.01$ and (1), each covering both complete range intervals, give

$$\sin^2\theta\ge0.0003168815744976>4(.005)^2.$$

Proposition 1 therefore supplies a finite range estimator on every design in this interval. An attempted single coarse enclosure of the whole interval fails; the coarse failure is retained. Refinement is a proof improvement, not evidence that the original physical class was impossible. The acquisition has a fixed-relative-error budget relative to $S_{\min}$; it is **not** a fixed-transmit-power optimum.

The local limit is

$$\lim_{w\to z}\frac{\sin^2\theta(z,w)}{(z-w)^2}=\frac{2z^2(z^2+4)}{(z^4+z^2+3)^2}.$$

Its derivative has the sign of $-(z^2-1)(z^2+1)(z^2+6)$, hence it peaks at $z=1$. This recovers the inherited local range feature but does not establish a universal optimum for material recovery, logarithmic range, antenna cost, or class C.

### Proposition 3: task-wise splitting despite exact internal scattering

In class D, take materials $\epsilon=2$ and $\epsilon'=2.3$, the same range, and electric-mode coefficients $t_E(\epsilon)$ computed by exact sphere theory. Choose $g=1$ and $g'=t_E(2)/t_E(2.3)$. Then the entire projected electric-dipole tensor observations are equal. Exact rational enclosures give

$$0.6832304596833457\le |g'|^2\le0.6832304596833458,$$

which lies strictly inside the allowed squared-gain annulus. Material separation $.3>2(.1)$ defeats any uniform $.1$ material guarantee. Range nevertheless has the nonempty finite window above under its stated noise budget. Thus geometry can be recoverable while material is not, even though the material coefficient includes exact internal scattering.

This is an exact class-D counterexample, not a Maxwell counterexample for all measured modes of two lossy interacting spheres.

### 3.1 How range calibration assists material readout

The dipole tensor norm, apart from its scalar material factor, is proportional to

$$Q(z)=\frac{\sqrt{2(z^4+z^2+3)}}{z^3},\qquad \left|\frac{d\log Q}{dz}\right|=\frac{z^4+2z^2+9}{z(z^4+z^2+3)}.$$

At $z=1$ the latter equals $2.4$. A range tolerance of $.02$ therefore need not be sufficiently tight for a sub-percent amplitude budget: even its first-order effect is about $4.8\%$. For a rigorous finite budget, use the supremum of the displayed derivative on the entire range interval and exponentiate its product with the range error. The benefit of a geometry window is to bound this nuisance effect, not to manufacture material information that the channel lacks.

## 4. Exact interval completion for the restricted modal theorem

Write $j=j_1$, $D_j=j'+j/z$. The corrected identities are

$$D_j'=-j'/z-j+j/z^2,\qquad D_j''=-j''/z-j'+2j'/z^2-2j/z^3.$$

We avoid square-root enclosure problems by using $s=\epsilon x^2$ and

$$J(s)=\sum_{n\ge0}\frac{(-1)^n(2n+2)s^n}{(2n+3)!},\qquad D(s)=\sum_{n\ge0}\frac{(-1)^n(2n+2)^2s^n}{(2n+3)!}.$$

Then $j_1(z)=zJ(z^2)$ and $D_j(z)=D(z^2)$. After cancelling the common square-root factor in the electric Mie coefficient, set

$$u=\epsilon xJ(\epsilon x^2),\quad v=D(\epsilon x^2),\quad N=uD_j(x)-j_1(x)v,\quad D_0=uD_y(x)-y_1(x)v.$$

The convention is

$$a=\frac{N}{N+iD_0},\qquad t_E=-a,\qquad q=\frac{N}{D_0}=\frac{ia}{1-a},\qquad h=|t_E|=\frac{q}{\sqrt{1+q^2}}.$$

All interval endpoints, coefficients and arithmetic are rational. Sixteen polynomial terms are followed by a tail bound of twice the first omitted absolute term. For differentiated series of order $d\le2$, $s\le1$, and $n\ge16$, the successive absolute-term ratio is

$$\frac{(2n+4)^r}{(2n+2)^r}\frac{n+1}{n+1-d}\frac{s}{(2n+5)(2n+4)}<\frac12,\quad r=1,2.$$

This proves the whole infinite tail bound. External sine and cosine use twenty Taylor terms with their decreasing-term remainders. Consecutive closed material boxes of width $1/1000$ exactly cover the domains, including shared endpoints. Every box proves $q,q',D_0>0$. We use the conservative bound

$$h'=\frac{q'}{(1+q^2)^{3/2}}\ge\frac{q'_{\rm lo}}{(1+q_{\rm hi}^2)^2}.$$

Integer division exports lower bounds downward and upper bounds upward; there is no float conversion in the proof kernel, no asserted endpoint minimum, and no midpoint-Lipschitz halfwidth issue.

| Mode | Exact domain | Boxes | $h'$ lower | $h$ upper |
|---|---|---:|---:|---:|
| $x=1/5$ | $[3/2,4]$ | 2500 | .0004549905509324 | .0026882355258206 |
| $x=3/20$ | $[2,5]$ | 3000 | .0001399328571476 | .0012933206367973 |

These constants improve the enclosure, not the underlying physical experiment. The inherited Gaussian two-world example is not rebranded as a new theorem.

## 5. Finite transverse material information from E/M relative response

Let $a=q_E/(q_E+i)$ and $b=q_M/(q_M+i)$ be the outgoing electric and magnetic coefficients; both measured $T$ coefficients have the same minus sign. Their ratio $r=t_M/t_E$ cancels a common complex gain. For lossless material,

$$\Re r=\frac{q_M}{q_E}\frac{1+q_Eq_M}{1+q_M^2}.$$

The same exact interval engine, with explicit outward compression to rational denominator $10^{30}$, covers every material box and proves

| Size | $d\Re r/d\epsilon$ lower | $|r|$ upper |
|---|---:|---:|
| $1/5$ | .0012718145885011 | .0079899593807010 |
| $3/20$ | .0007380750428085 | .0052511493011775 |

Let the two channels be $Y_E=gt_E+e_E$ and $Y_M=gt_M+e_M$, with $|e_E|\le B_E$, $|e_M|\le B_M$, $|t_E|\ge h_{E,\min}$ and $B_E<g_{\min}h_{E,\min}$. Then

$$\left|\frac{Y_M}{Y_E}-r\right|\le\frac{B_M+R_{\max}B_E}{g_{\min}h_{E,\min}-B_E}.$$

Clip the real part of this ratio to the range of the monotone real-ratio function, then invert. The mean value theorem supplies the **finite**, not merely differential, bound

$$|\widehat\epsilon-\epsilon|\le\frac{B_M+R_{\max}B_E}{(g_{\min}h_{E,\min}-B_E)m_r}.\tag{2}$$

Thus a material change creates a measured projective change that survives complex-gain profiling in this restricted known-geometry readout. Differential E/M channel gain, geometry error, leakage, and omitted modes must enter $B_E,B_M$; they cannot be treated as a common scalar.

For the homogeneous-sphere Born model at the same known geometry and constitutive prior, every linear channel is proportional to the same contrast $\chi$. The E/M channel ratio is independent of $\chi$ wherever defined. Consequently the derivative used by (2) vanishes in that Born model. This comparison does not change the loss prior or electronics to make Born artificially weaker.

The existence of an informative magnetic mode does not imply it is cheap to read. Under equal absolute proper-complex noise per E and M channel, zero model error, and the conservative rational noise-event bounds below, (2) is satisfied by 9528 and 23482 independent repeats respectively. These are sufficient counts for this estimator and budget, not necessary sample lower bounds or optimal hardware costs.

## 6. Reference readout, leakage, and physically interpretable budgets

For class M, the inherited amplitude-normalization estimator obeys

$$|\widehat\epsilon_i-\epsilon_i|\le\frac{B_i+H_iB_a}{(a_{\min}-B_a)m_i}.\tag{3}$$

This follows from reverse triangle inequality, a positive reference denominator, clipping, and inverse Lipschitz continuity. It is not a new GLS identity.

A physically meaningful electronic reference can be $Z=g+\eta+d_g$, with known coupling/normalization. Then $||Z|-|g||\le|\eta|+|d_g|$. This replaces an ideal real amplitude oracle by a measured complex channel. It assumes the reference and target channels share the gain being calibrated; a different cable/path response requires a separate bounded calibration factor.

With sixteen independent repeats of each of the two modal channels and the reference, take upper bounds on single-read complex standard deviations $\sigma_1\le2.666409\cdot10^{-6}$, $\sigma_2\le1.804800\cdot10^{-6}$, $\sigma_Z\le.001$. The reference noise radius is $.001(2.629)/4$. A positive rational Taylor partial sum proves $\exp(2.629^2)>1000$, so each averaged complex noise lies inside its radius with probability at least $.999$. A union bound gives the three-channel event probability at least $.997$ without requiring inter-channel independence; within-channel averaging does require the stated noise variance reduction.

An explicit sufficient specification is:

| Contribution | Required bound |
|---|---:|
| Fixed shared-reference drift $|d_g|$ | .0005 |
| Relative modal channel gain error $|\kappa_i|$ | .0005 |
| Cross-mode leakage coefficient $|\lambda_{ij}|$ | .0003 |
| Size/constitutive allocation, each normalized amplitude | $10^{-6}$ |
| Unresolved-mode allocation, each | $10^{-7}$ |
| Background-subtraction bias, each | $10^{-6}$ |

The deterministic mode budget is bounded by

$$b_i\le g_{\max}(.0005H_i+.0003H_j)+10^{-6}+10^{-7}+10^{-6}.$$

The rationally exported bounds are $b_1\le4.2651424425\cdot10^{-6}$ and $b_2\le3.9164137202\cdot10^{-6}$. Total reference error is at most $.00115725$. Equation (3) then yields material errors at most **.0267923733035967** and **.0629779399380477** on the probability event.

A separate rational size certificate covers $|x-x_0|/x_0\le1/20000$ and the full material domains. It proves $|\partial h/\partial x|\le.0407206984353109$ and $.0260702567554705$. Gain-weighted size biases are at most $5.090087305\cdot10^{-7}$ and $2.444086571\cdot10^{-7}$, inside the stated $10^{-6}$ allocations when the remaining constitutive assumptions are exact. The size tolerance applies to the product $ka$: radius and frequency uncertainties must share it. At a nominal 10 mm radius and exact frequency, it corresponds to 0.5 micrometre radius error; this is a required specification, not demonstrated manufacturing accuracy.

Under the **inherited real-amplitude reference noise event**, without repeats or extra reference drift, the new rigorous modal constants permit a conservative second-channel bias of at least $1.4483933265\cdot10^{-6}$, rather than the earlier sufficient $5.78\cdot10^{-7}$. The older certificate is restrictive but not refuted; the older constants and files are preserved. Changing to the complex reference is a separately declared statistical experiment.

**Attainability boundary.** The inequalities define a nonempty numerical specification, with 48 complex readouts for the three-channel sixteen-repeat experiment. They do not establish that an actual antenna/probe/reference apparatus meets it. Full-angular modal projection and incident/outgoing separation can define the intended channels mathematically. Finite aperture, probe perturbation, channel normalization, reference path drift, finite modal truncation and realistic loss must be independently bounded. Repetition reduces neither common drift nor common model bias. Hardware attainability and minimum total physical cost remain unresolved.

## 7. Uniform class-C forward stability: a legitimate route to model bounds

Let $p=\chi E\in L^2(D;\mathbb C^3)$, with $D$ the union of the two balls, and use the distributional Maxwell VIE

$$A_\chi p=(\chi^{-1}-G_k)p=E^{\rm inc},\qquad G_k=(k^2I+\nabla\nabla)\frac{e^{ik|\cdot|}}{4\pi|\cdot|}*.$$

The distributional local term is included, not omitted as in an unqualified point-dipole sum. Standard VIE formulations and numerical-range analysis are discussed in [R4,R5]. We claim a domain-specific explicit bound, not novelty for VIE coercivity or the residual principle.

### Proposition 4: uniform coercivity on the declared two-sphere class

For every allowed lossy material pair,

$$\Re\langle p,A_\chi p\rangle\ge\frac{43}{500}\|p\|_{L^2(D)}^2,\qquad \|A_\chi^{-1}\|\le\frac{500}{43}.\tag{4}$$

**Proof.** Extend $p$ by zero. The static part $G_0=\nabla\nabla(4\pi r)^{-1}*$ has Fourier multiplier $-\xi\xi^T/|\xi|^2$, so it is self-adjoint nonpositive. The Hermitian real part of the dynamic difference is an integrable symmetric real kernel. Its transverse and radial eigenvalues, multiplied by $4\pi r^3$, are respectively

$$f_t(t)=(t^2-1)\cos t-t\sin t+1,\quad f_r(t)=2(\cos t+t\sin t-1),\quad t=kr.$$

Since $(\cos t+t\sin t-1)'=t\cos t$, $|f_r(t)|\le t^2$. Also $f_t'(t)=t\cos t-t^2\sin t$, giving $|f_t(t)|\le t^2/2+t^4/4\le t^2$ when $t^2\le2$. Globally, writing $f_t=t^2\cos t-(\cos t+t\sin t-1)$ gives $|f_t|\le3t^2/2$.

Within either ball, $2ka_i<\sqrt2$. Newton's potential formula and the Schur test bound each real dynamic self block by $k^2a_i^2/2$. For distinct balls of center distance $d$, the same test gives

$$\|\Re G_{k,12}-G_{0,12}\|\le\frac{k^2(a_1a_2)^{3/2}}{2\sqrt{(d-a_1)(d-a_2)}}.$$

The factor $3/2$ in the cross-kernel bound and the ball volume $4\pi a^3/3$ are both included. The lower center-distance bound $d>.1167$ is verified by squaring exact rationals. The cross norm is below $.0485$.

The real parts of inverse contrasts are minimized at the upper endpoints of the real contrast intervals because those intervals exceed their fixed losses. After subtracting the self-block bounds, the two diagonal lower bounds are $.1348500033330000$ and $.1487109436025620$. The two-by-two quadratic-form lower bound, or Gershgorin's inequality, is greater than $\min(d_1,d_2)-.0485>.086=43/500$. All last inequalities are checked as rational inequalities. The bounded sesquilinear form is coercive, so Lax-Milgram gives existence, uniqueness and (4).

### Corollary: residual-to-observation transfer, not a completed certificate

For any actual continuous trial polarization $\widetilde p$ and its actual continuum residual $r=E^{\rm inc}-A_\chi\widetilde p$,

$$\|p-\widetilde p\|_{L^2}\le\frac{500}{43}\|r\|_{L^2}.$$

All receivers stay at least $.503$ m from either ball. For $N=12$ three-component receivers and total volume $V$, the observation operator satisfies

$$\|S\|\le\frac{\sqrt{3NV}}{4\pi}\left(\frac{k^2}{d_{\min}}+\frac{2k}{d_{\min}^2}+\frac2{d_{\min}^3}\right)<6.$$

The factor $\sqrt3$ safely bounds the dyadic Frobenius norm by its spectral norm. Squaring the expression and using $\pi>3.14159$ makes the bound $6$ an exact rational check. Hence, including $|g|\le1.25$,

$$\|gS(p-\widetilde p)\|<88\|r\|_{L^2}.\tag{5}$$

For four illuminations use the product/Frobenius norm; the same block-diagonal operator constant applies, not an extra unexplained factor four. Whitening multiplies (5) by the appropriate noise metric norm.

**What remains missing.** The present DDA grid difference is not $\|r\|_{L^2}$. The multipole order difference is not $\|r\|_{L^2}$ either. We have not certified a continuum residual or a finite material separation over the class-C world pairs. Equation (5) supplies the previously missing stable transfer constant, but does not close those two obligations.

A concrete next residual construction is available: for an approximate multisphere multipole field, use exact single-sphere interior/exterior solutions. Within ball $i$, the residual equals incident field plus the other sphere's outgoing field minus the regular incident expansion used by the approximate solution. A rigorous translation/quadrature remainder and a high-mode tail bound would then bound an actual residual without estimating a singular volume integral by two solver differences. This construction is a proof target, not an implemented certificate in this version.

## 8. Classification, antenna decisions, and the limits of amplification

For a numerical diagnostic, the real whitened local sensitivities are material columns $A$, actual receiver-shift derivative $B$, and two electronic columns $C=[\mathcal R(f),\mathcal R(if)]$. Only $B=\partial_s\mathcal R(gf)$ from the physical field model represents geometry. A generic extra vector does not. Projection against $[B,C]$ diagnoses local material visibility; it does not supply Propositions 1 or 4 or label an arbitrary observed discrepancy causally.

Distinct physical mechanisms must remain separate:

1. Internal scattering changes each sphere's material response. If the only observed effect is one scalar coefficient, a free common gain can still absorb it (Proposition 3).
2. Electric/magnetic relative response can change the projective observation and cancel a common gain (Section 5), but the magnetic channel can be too weak at the stated absolute noise.
3. Higher multipoles can add different material responses; their measurable magnitude, differential calibration, and truncation bounds must be budgeted.
4. Object diversity can defeat a common-gain ambiguity, but known constitutive loss already contributes such information and is not an electronic reference.
5. Inter-object rescattering couples the two objects. It is present in class C, but no universal benefit follows from that fact alone.

For the fixed-loss Born model, $\mu=gB(u+i\ell)$ with $\ell=(.03,.05)$. The inherited legal ambiguity uses $u=30\ell$, $u'=40\ell$, gains $1$ and $(30+i)/(40+i)$, hence materials $(1.9,2.5)$ and $(2.2,3.0)$. Their Born means are identical for the same geometry and $B$. Away from the parallel class, known loss can already constrain scale; declaring all Born material information absent would be false.

Antenna design should therefore choose geometry-sensitive vector channels, material-sensitive relative modes, and a physically shared reference under the same noise and bias budget. A useful design criterion is a **proved** bad-world separation margin minus twice the error radius, divided by a declared acquisition cost. This is a research objective, not a validated new policy. A fixed or random acquisition with the same complex scalar count remains a mandatory baseline; hardware cost parity requires more than count parity.

For invertible linear processing $L$, transforming both signal and covariance preserves the Mahalanobis distance: $(L\Delta)^*(L\Sigma L^*)^{-1}L\Delta=\Delta^*\Sigma^{-1}\Delta$. An amplifier cannot remove an exact gain/material ambiguity. In the explicitly assumed scalar chain $Y=A(S+N_{\rm pre})+N_{\rm post}$, output SNR is $|A|^2P_S/(|A|^2P_{\rm pre}+P_{\rm post})$: gain helps against downstream noise, saturates against pre-amplifier noise, and amplifies common bias. No undefined `Gain-F` formalism is invoked. Any proposed gain/noise-factor figure of merit must first specify this chain and its constraints.

## 9. Independent numerical results and honest acceptance

The registered V2 stream is separate from V1, with seed 2026091107, fixed three starts, and a DDA generator distinct from the numerical vector-spherical-wave inverse implementation. The main inverse uses $l_{\max}=3$. Single-sphere boundary conditions, incident expansion, quadrature refinement and multipole refinement were checked before recovery. These checks are diagnostics, not a proof of class-C continuum accuracy.

| Method | Material successes / 12 | Median maximum material error | Median sensor relative residual |
|---|---:|---:|---:|
| No reference | 2 | .533669 | .011996 |
| Noisy-reference GLS | 7 | .062018 | .012742 |
| Fixed extra EM scalar | 2 | .524025 | .011998 |
| Random extra EM scalar | 2 | .533840 | .011998 |

Success requires both material errors at most .1. Reference GLS reduces gain/material tradeoff but may increase sensor-only residual because its objective also fits an independent reference. It neither directly fixes receiver geometry nor makes a structurally incorrect forward model correct. Fixed losses supply prior material information independently of the reference. For scene 8, reference GLS has material errors approximately (.2264,.2941), despite a sensor relative residual of only .01547. Such a result fails the material task.

The finer DDA grid is not monotonically closer to the high-order sphere solver in the first four diagnostics. This is evidence against using either difference as an error upper bound. Raw observations, all starts, optimizer failures, hashes, the interrupted checkpoint and its deterministic resume are retained in the evidence archive. This is a 12-scene development comparison, not a final test, not a population success rate and not proof of a new algorithmic advantage.

## 10. Submission gate

The restricted rational proofs, task-wise finite example, explicit class-C coercivity and independent numerical diagnostic are usable research materials. A strong TAP claim about simultaneous lossy-sphere material recovery still needs: a certified continuum residual/forward enclosure; covered finite material separation or certified realistic competing worlds under that uncertainty; verified physical readout error budgets; and a locked independent final recovery comparison. The nearest-prior originality audit is not closed for all directly relevant papers. No title or abstract should claim these obligations have already been met.

## Appendix A. Sharp finite amplitude boundary and a second splitting example

A stronger finite claim was present in a pre-existing draft report on the research branch. It was not accepted from that report alone. This appendix independently supplies its proof and a new rational concavity/secant audit. Numerical bounds here use an outward noise radius and therefore need not equal the draft's last decimals.

Consider the explicitly restricted **amplitude statistical experiment**

$$W_i=a h_i(\epsilon_i)+e_i,\quad |e_i|\le B_i,\qquad Z=a+e_a,\quad |e_a|\le B_a,$$

with rectangular material prior $[L_i,U_i]$, $h_i>0$ strictly increasing and concave, $a\in[a_-,a_+]$, product interval error sets, and $0<2\delta<U_i-L_i$. Set $d_a=\min(2B_a,a_+-a_-)$. Assume, for every other channel needed in a counterexample,

$$a_-h_j(U_j)\ge(a_-+d_a)h_j(L_j).\tag{A1}$$

This is a quantitative cross-channel amplitude-overlap condition, not an assumed generic nuisance direction.

### Proposition 5: finite-secant necessary and sufficient condition

A uniform coordinatewise material estimator with error at most $\delta$ exists for this experiment if and only if, for every channel $i$,

$$G_i:=a_-h_i(U_i)-(a_-+d_a)h_i(U_i-2\delta)\ge2B_i.\tag{A2}$$

**Sufficiency.** Reference-compatible worlds satisfy $|a-a'|\le d_a$. For a bad material coordinate $\epsilon_i>\epsilon'_i+2\delta$, the smallest possible signed amplitude separation is bounded below by

$$a_-\{h_i(\epsilon_i)-h_i(\epsilon'_i)\}-d_a h_i(\epsilon'_i).$$

Concavity makes the fixed-length increment smallest at the upper material endpoint, and increasing $h_i$ makes the subtracted term most adverse there. Strict material separation makes this quantity strictly larger than $G_i$. Thus (A2) prevents amplitude intervals from intersecting. Apply Proposition 1. The equality case in (A2) is sufficient because the bad-material inequality is strict.

**Necessity.** If $G_i<2B_i$, take $\epsilon_i=U_i$, $\epsilon'_i=U_i-2\delta-\eta$ for sufficiently small positive $\eta$, $a=a_-$, and choose $a'\in[a_-,a_-+d_a]$ so the absolute amplitude separation is at most $2B_i$. For $B_i>0$ it can be made strictly smaller; for $B_i=0$ and $G_i<0$, choose the exact zero crossing of the compensation. Continuity supplies these choices. In another channel choose $\epsilon'_j=L_j$ and $h_j(\epsilon_j)=(a'/a_-)h_j(L_j)$, possible by (A1). Those channels match exactly. Reference intervals overlap because $|a-a'|\le2B_a$. Their observation midpoints are compatible with both worlds, whose material boxes are disjoint. No uniform estimator exists.

The need for an amplitude-only experiment is essential. Taking magnitudes of coherent data yields a sufficient procedure, but collisions after that processing do not prove impossibility for the original complex experiment. An electronic-reference-only intervention likewise cannot resolve any bad pair with identical gain; a perfect full-complex-gain reference is sufficient exactly when all remaining fixed-gain bad pairs are distinguishable under the base error sets. Necessity and sufficiency must always name the intervention and error class.

### A.1 Concavity and the actual checked endpoint quantity

The new `certify_secant.py` covers all 5500 material boxes, uses $u''=x^3(2J'+sJ'')$, $v''=x^4D''$, and bounds

$$q''=\frac{N''D_0-ND_0''}{D_0^2}-\frac{2D_0'(N'D_0-ND_0')}{D_0^3}.$$

It proves $q''\le-.0001492172302260$ and $q''\le-.0000396694431295$ for the two modes. Combined with $q,q'>0$ this gives

$$h''=\frac{q''}{(1+q^2)^{3/2}}-\frac{3q(q')^2}{(1+q^2)^{5/2}}<0.$$

Endpoint evaluations now follow from a proved concavity statement, not a grid observation. Rational square-root bounds use integer square roots and explicitly verified outward endpoints. With $B_a=.003291$ and $\delta=.1$, the checked finite separations satisfy

$$G_1\in[.0000535160055495,.0000535160055496],$$
$$G_2\in[.0000132879851733,.0000132879851734].$$

Condition (A1) is checked at both other-channel endpoint ranges. Using the inherited modal noise, upward standard-deviation bounds and $\sqrt{\log1000}<2.629$, (A2) leaves additional second-channel amplitude bias at least $1.8991770120\cdot10^{-6}$. This is a different, stronger finite-set estimator bound than the inverse-slope certificate (3). Neither number certifies an apparatus.

For $2B_a<a_+-a_-$, the exact allowable reference precision can equivalently be written

$$B_a\le\min_i\frac{a_-[h_i(U_i)-h_i(U_i-2\delta)]-2B_i}{2h_i(U_i-2\delta)},$$

provided (A1) and the nonnegative right-hand side hold. A negative right side means that even perfect amplitude calibration cannot meet this amplitude-experiment specification at that noise, not that all coherent acquisitions are impossible.

### A.2 Fixed-raw-error material windows

For an ideal full-angular tangential electric-dipole projection, the outgoing radial factor has squared magnitude

$$c(z)^2=z^{-2}-z^{-4}+z^{-6},\quad [c(z)^2]'=-2[(z^2-1)^2+2]/z^7<0.$$

Let geometry be known, $z=kR>1/5$, both raw amplitude-channel errors be bounded by $B_0=10^{-6}$, and reference error remain $.003291$. Dividing by the known radial factor gives $B_i(z)=B_0/c(z)$. Equation (A2) defines the exact amplitude-task window by $c(z)G_i\ge2B_0$. Rational bisection propagates both endpoints of the $G_i$ enclosures and locates the two respective cutoffs near

$$z_1^*=26.7393102838,\qquad z_2^*=6.5683473710.$$

Every range below the certified lower endpoint is covered; above the upper endpoint the amplitude-experiment necessity construction applies. The tiny transition interval is retained in `secant_certificate_final.json`. The two material tasks therefore have different distance windows, and their simultaneous task uses the intersection. There is no interior optimum in this fixed-raw-error example. This is distinct both from the class-D range window at fixed relative error and from a class-C unknown-geometry material guarantee.
