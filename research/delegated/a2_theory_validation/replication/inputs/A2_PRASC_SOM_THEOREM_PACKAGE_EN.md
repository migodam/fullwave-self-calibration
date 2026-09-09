# Rank, Bias, and Acquisition Certificates for Phase-Preserving Full-Wave Self-Calibration

**Research replacement manuscript and theorem supplement. 5 September 2026.**

**Status:** finite-dimensional theoretical results and an implementable, conditional control algorithm. No new experiments were run. No superiority over a properly tuned direct joint full-wave solver is established. The inherited numerical results are source-reported negative controls, not validation of this algorithm.

## Abstract

Unknown array geometry can corrupt coherent inverse-scattering measurements through a low-dimensional, physically structured phase error. Retaining phase is potentially preferable to discarding it, but current-space extensions can also absorb the signatures needed for calibration. We develop a finite-dimensional accounting framework that separates physical state elimination from free-current nuisance enlargement. In particular, the deterministic cutoff of the original subspace optimization method (SOM) is distinguished from the dimension of its freely optimized current space. After whitening, realification, and joint gauge reduction, map-information retention and pose visibility are expressed as the two normalized Schur complements of a common tangent experiment. We prove information contraction under matched phaseless observation even after nuisance elimination, and characterize the exact loss caused by nested free-current spaces. We then give calibration-neutral admission conditions for additional current directions, an exact shared-map acquisition innovation formula, and a bias-aware local risk bound. Together, these results provide a checkable information budget for rank enlargement and acquisition changes. An a posteriori state-residual certificate connects reduced current computations to the physical forward model without treating a state penalty as an independent measurement. A safeguarded, phase-preserving rank-adaptive algorithm is proposed, with conditional descent guarantees and explicit failure outputs. The results do not establish a universal nonlinear calibration basin or an advantage exclusive to SOM. They support a limits-and-design contribution whose method-level value must be tested against matched direct inversion.

## 1. Introduction

The response of an electromagnetic imaging system depends jointly on material and acquisition geometry. In coherent measurements, a displacement changes the propagation phase, receiver sampling, and, when the transmitter moves, the induced current. These effects are not interchangeable. A receiver-side sampling perturbation may admit a data-equivalent current representation without being a physical change of the induced current. Such a representation is automatic if an unrestricted current-to-data matrix is surjective and therefore cannot establish self-calibration.

The central question is not whether pose and material can be included in one optimizer. It is whether a current-space approximation can be controlled so that its computational or approximation benefit does not consume the information needed to estimate geometry. SOM provides a natural family of current coordinates, but its ranks must be interpreted correctly. Increasing the number of data-determined singular components in classical SOM decreases the dimension of the complementary ambiguous current space. Increasing the dimension of an independently optimized correction space does the opposite. A monotonicity result for the latter cannot be attributed to the former without a separate derivation.

We retain complex coherent measurements and compare them with phaseless observations only through a matched statistical experiment. We separately analyze the physical reduced-state model, a free-current diagnostic envelope, and approximate state solves. This distinction makes it possible to state both useful certificates and impossibility results. In particular, a zero visibility score for an enlarged diagnostic envelope need not imply nonidentifiability of the physical model, while a positive score is a conservative first-order certificate when the physical nuisance space is contained in that envelope.

The contribution is not the use of projections, Schur complements, canonical correlations, or continuation. It is a typed, bias-aware control problem with explicit admission and acquisition conditions. Its scientific status remains conditional: mathematical safety for a declared tangent experiment is not the same as nonlinear recovery from an arbitrary initialization.

## 2. Model, metrics, and gauge

Let the shared material parameter be $\chi\in\mathbb R^q$, the shared trajectory parameter be $x\in\mathbb R^p$, and the frequency-dependent contrast be $\chi_\ell=f_\ell(\chi)\in\mathbb C^n$. The map $f_\ell$ is a declared material/dispersion model. Each acquisition has its own current $j_\ell\in\mathbb C^n$:

$$
y_\ell=S_\ell(x)j_\ell+\varepsilon_\ell,\qquad
M_\ell j_\ell=b_\ell,\qquad
M_\ell=I-D_{\chi_\ell}D_\ell,\quad b_\ell=D_{\chi_\ell}e_\ell(x).
$$

The current and internal-field spaces are physically different spaces, even when both discretizations have dimension $n$. On a world-fixed grid with a fixed homogeneous background, $D_xD_\ell=0$.

Writing $T_\ell=Df_\ell(\chi)$ and $G_{\chi,\ell}=D_{E_{\rm tot,\ell}}T_\ell$, the independent-current tangent is

$$
\delta y_\ell=S_\ell\delta j_\ell+H_{S,\ell}h,
\qquad
M_\ell\delta j_\ell-G_{\chi,\ell}\delta\chi-H_{D,\ell}h=0.
$$

Here $H_{S,\ell}h=(D_xS_\ell[h])j_\ell$ and $H_{D,\ell}h=D_{\chi_\ell}D_xe_\ell[h]$. If $M_\ell$ is invertible, the exact physical reduced derivatives are

$$
A_\ell=S_\ell M_\ell^{-1}G_{\chi,\ell},\qquad
B_\ell=H_{S,\ell}+S_\ell M_\ell^{-1}H_{D,\ell}.
$$

There is no additional free-current block in this exact reduced model. A model with $K_rc+A\delta\chi+Bh$ must declare $c$ as an additional correction/nuisance or explicitly label the expression as an envelope.

For proper complex noise with $\mathbb E[\varepsilon\varepsilon^*]=\Sigma$, let $W^*W=\Sigma^{-1}$. The real embeddings of real and complex parameter Jacobians are respectively $\mathcal E(WA)$ and $\mathcal R(WK)$. These embeddings have noise covariance $I/2$; multiply them by $\sqrt2$ for unit real covariance. Otherwise their Grams are half the Fisher information. Improper or cross-acquisition correlated noise requires full real-covariance whitening, not independent block whitening.

Let $\mathcal G\subset\mathbb R^q\oplus\mathbb R^p$ be the joint gauge tangent. For a positive definite parameter metric $M_\theta$, choose a complement $Z$ with $Z^TM_\theta Z=I$ and first form $[A\ B]Z$. An absolute pose is not a target on the quotient unless it is anchored. More generally, let $T=[A\ B]Z$, let $L$ be a well-defined target on quotient coordinates, let $R$ be a right inverse of $L$, and let $N$ span $\ker L$. The target and nuisance tangents are $TR$ and $TN$. Replacing $R$ by $R+NC$ does not change $P_{\operatorname{Ran}(TN)^\perp}TR$. Only after this operation or explicit gauge fixing may separate map and pose blocks be used.

For a current metric $M_c\succ0$ and data metric $N_y\succ0$,

$$
S^{\dagger}_{M_c,N_y}
=M_c^{-1/2}(N_y^{1/2}SM_c^{-1/2})^\dagger N_y^{1/2},
\qquad
P_U^{M_c}=U(U^*M_cU)^{-1}U^*M_c.
$$

The data projector is $SS^{\dagger}_{M_c,N_y}$ and is orthogonal in $N_y$. For a retained basis $U$, set $G=U^*M_cU$ and $T=N_y^{1/2}SUG^{-1/2}$. The minimum-current-norm least-squares coefficient for data $v$ is $G^{-1/2}T^\dagger N_y^{1/2}v$. This follows by the change of variable $a=G^{1/2}c$ and ordinary Euclidean least squares. The resulting receiver lift must be applied to $H_S$, not to the total $B$ while retaining $H_D$ again.

## 3. Information contraction after nuisance elimination

### Theorem 1: matched observation and efficient information

Assume a dominated statistical model differentiable in quadratic mean, square-integrable mean-zero scores, and a parameter-independent channel from $Y$ to $Z$. The target and nuisance models, including gauge restrictions, are the same in both experiments. Then

$$
s_Z=\mathbb E[s_Y\mid Z],\qquad J_Z\preceq J_Y,
\qquad J^{\rm eff}_{Z}\preceq J^{\rm eff}_{Y}.
$$

The last inequality holds for either map or pose as target, including singular nuisance information. In a target direction $v$, equality in the efficient inequality holds if and only if the coherent efficient score in that direction is measurable with respect to $Z$.

**Proof.** For every bounded test function $g(Z)$, differentiating its expectation through the channel gives $\mathbb E[g s_Z]=\mathbb E[g s_Y]$. The defining property of conditional expectation gives the score identity. Conditional variance gives $J_Y-J_Z=\mathbb E[\operatorname{Cov}(s_Y\mid Z)]\succeq0$.

Let $T$ denote conditional expectation onto functions of $Z$. Let $\mathcal H_Y$ be the closed nuisance-score space and $\mathcal H_Z=\overline{T\mathcal H_Y}$. For a target score $s_v$, let $e_Y=(I-P_{\mathcal H_Y})s_v$. Since the removed nuisance component maps into $\mathcal H_Z$, the transformed efficient score is $e_Z=(I-P_{\mathcal H_Z})Te_Y$. Orthogonality gives

$$
\|e_Y\|^2-\|e_Z\|^2
=\|e_Y-Te_Y\|^2+\|P_{\mathcal H_Z}Te_Y\|^2\ge0.
$$

This proves the matrix inequality by all directional quadratic forms. Equality implies $e_Y=Te_Y$. Conversely, if this measurability holds, then for $n\in\mathcal H_Y$, $\langle e_Y,Tn\rangle=\langle e_Y,n\rangle=0$, so the second term is also zero. For finite nuisance parameters, minimizing $\|s_v-b^Ts_\eta\|^2$ gives the generalized Schur complement. A null nuisance score has zero covariance with every target score; therefore the range condition needed for the pseudoinverse Schur formula is automatic. No inverse Fisher matrix is required. QED.

**Counterexamples and boundaries.** If $Y=a e^{ix}+\varepsilon$ with proper circular Gaussian noise and known $a>0$, $|Y|$ is independent of $x$, whereas coherent pose information is $2a^2/\sigma^2$. If an unknown nuisance phase $\eta$ changes the mean to $ae^{i(x+\eta)}$, the efficient information about $x$ is zero in both experiments. Positive equality also occurs for a zero-mean circular Gaussian scale family: its score is a function of $|Y|^2$. An independently engineered intensity sensor with a different noise variance is not ordered by this theorem.

The raw experiment matters. Phaseless SOM may use total-field intensity. If the incident field depends on unknown pose, $Z=|Y^{\rm sca}+e_{\rm rx}(x)|^2$ is not a parameter-independent channel from scattered-field data. Use raw coherent total-field data as the parent experiment, or explicitly include the reference measurement and its noise. A known reference can retain sensitivity to the phase of the scattered field: for real nonzero a and c, |c+a exp(ix)|^2 has derivative -2ac sin(x). Thus a phase-only scattered-field perturbation is not necessarily a phase-only perturbation of the total field.

## 4. Phase precision versus a local branch

For the outgoing convention, $G_2(r)=iH_0^{(1)}(kr)/4$ and $G_3(r)=e^{ikr}/(4\pi r)$. Large-argument Hankel expansions and their differentiated expansions give

$$
\partial_r\log G_2=ik-\frac1{2r}+O((kr^2)^{-1}),\qquad
\partial_r\log G_3=ik-\frac1r.
$$

Thus the relative phase derivative is $k$ to leading order. The far-field Maxwell dyadic has the leading transverse factor $I-\widehat r\widehat r^T$ multiplying $G_3$; the same conclusion requires nonvanishing polarization coupling. The normalized pose Jacobian is order $k$ only when path derivatives and relative amplitudes stay bounded and cancellation is excluded. At fixed coherent SNR the raw Fisher scale is order $k^2$. Efficient information needs, additionally, a uniformly positive nuisance-separation margin. At fixed absolute noise, spreading and frequency-dependent source/material amplitudes can change these powers.

### Proposition 2: a conditional phase-basin statement

For the resolved single-path model

$$
f(h)=\sum_iw_i[1-\cos(k\ell_i^Th)],\qquad w_i>0,
$$

on a convex region with $|k\ell_i^Th|\le\gamma<\pi/2$ for every $i$,

$$
\nabla^2 f(h)\succeq k^2\cos\gamma\,L^TWL.
$$

If $L$ has full column rank, $f$ is strongly convex on this region.

**Proof.** Twice differentiating each cosine gives $k^2w_i\cos(k\ell_i^Th)\ell_i\ell_i^T$. Summing and bounding every cosine below by $\cos\gamma$ gives the claim. QED.

An individual path repeats after displacement $2\pi/(k|\ell_i^Tv|)$ along direction $v$. Multiple paths may remove exact common aliases but can also create destructive interference and extra minima. The proposition is not a full-wave global basin theorem.

If a genuine error confidence set obeys $h^T\Sigma_x^{-1}h\le c^2$ and path Taylor remainders are bounded by $C_{\rm path}\|h\|^2/2$, a sufficient phase gate is

$$
k_{\rm next}\left[c\max_i\sqrt{\ell_i^T\Sigma_x\ell_i}
+\frac{C_{\rm path}c^2\lambda_{\max}(\Sigma_x)}2\right]
+\delta_{\phi,\rm model}\le\gamma<\pi/2.
$$

This follows from Cauchy-Schwarz in the ellipsoid and the Taylor bound. The square root of $\lambda_{\max}(L\Sigma_xL^T)$ is a conservative replacement for the maximum row standard deviation. For a genuinely centered Gaussian error, $c^2$ may be the appropriate chi-square quantile. An inverse Gauss-Newton matrix is not generally a covariance upper bound; using it here without a coverage argument produces a surrogate, not a certificate. A wrong phase branch can have zero residual and arbitrarily small local covariance.

## 5. Dual tangent spectra

All matrices in this section are real, unit-noise whitened, and gauge-fixed. Set

$$
C=\operatorname{Ran}K_r,\quad \Pi=P_{C^\perp},\quad A_c=\Pi A,\quad B_c=\Pi B,
$$
$$
K_0=A_c^TA_c,\quad K_e=A_c^T(I-P_{\operatorname{Ran}B_c})A_c,
\quad J_x=B_c^T(I-P_{\operatorname{Ran}A_c})B_c.
$$

### Theorem 3: singular-support Schur duality

Let $a=\operatorname{rank}A_c$ and $b=\operatorname{rank}B_c$, and use only the supports of $K_0$ and $B_c^TB_c$. Let $Q_A,Q_B$ be orthonormal range bases and let $c_i$, $i\le\min(a,b)$, be the singular values of $Q_A^TQ_B$. The normalized map and pose information matrices have spectra

$$
\{1-c_i^2\}_{i=1}^{\min(a,b)}\cup\{1\}^{a-\min(a,b)},
\qquad
\{1-c_i^2\}_{i=1}^{\min(a,b)}\cup\{1\}^{b-\min(a,b)}.
$$

The common zero multiplicity is $\dim(\operatorname{Ran}A_c\cap\operatorname{Ran}B_c)$. The unit multiplicities are $a-\operatorname{rank}(Q_A^TQ_B)$ and $b-\operatorname{rank}(Q_A^TQ_B)$. The raw spectra of $J_x$ and $K_e$ need not coincide.

**Proof.** On their supports, write $A_c=Q_AR_A$ and $B_c=Q_BR_B$ with invertible reduced factors. Congruence by these factors reduces the two Schur forms to $I_a-CC^T$ and $I_b-C^TC$, where $C=Q_A^TQ_B$. Their nonzero correlation eigenvalues are the squared singular values of $C$, with the stated dimension padding. Equality $\|Q_Av\|=\|P_{\operatorname{Ran}B_c}Q_Av\|$ characterizes an intersection direction; zero projection characterizes an orthogonal direction. These observations give the zero and unit multiplicities. Congruence, unlike orthogonal similarity, does not preserve raw eigenvalues, proving the last caution. QED.

The generalized map eigenvalues are therefore squared sines of principal angles, but generalized eigenvalues on $\ker K_0$ are undefined rather than assigned zero or one.

A pose prior $\Lambda_x\succeq0$ yields

$$
K_{e,\Lambda}=K_0-A_c^TB_c(B_c^TB_c+\Lambda_x)^\dagger B_c^TA_c.
$$

Indeed, its quadratic form is the minimum of $\|A_cu-B_ch\|^2+\|\Lambda_x^{1/2}h\|^2$. This proves $K_e\preceq K_{e,\Lambda}\preceq K_0$ and monotonicity in the prior. Exact angle geometry is recovered in the augmented space using $[A_c;0]$ and $[B_c;\Lambda_x^{1/2}]$, not by pretending that a nonprojector in the original data space is a projector. With both map and pose priors, use columns $[A_c;\Lambda_\chi^{1/2};0]$ and $[B_c;0;\Lambda_x^{1/2}]$ and normalize by the prior-inclusive diagonal blocks. Deterministic regularization curvature is not independent data Fisher information.

### Proposition 4: visibility, dimensions, and the nonlinear boundary

Let $N=\operatorname{Ran}[K_r,A]$, $d=\dim N$, and $b=\operatorname{rank}B$. Then

$$
B_v=P_{N^\perp}B,\qquad
\dim\ker B_v=p-b+\dim(\operatorname{Ran}B\cap N),
$$
$$
\max(p-b,p+d-m)\le\dim\ker B_v\le p-b+\min(b,d).
$$

These bounds are sharp over matrices with the prescribed dimensions. Pose is identifiable in the linearized experiment exactly when $\operatorname{rank}B_v=p$. On the observable support its least-squares error is bounded by the data/model-error norm divided by the smallest positive singular value.

**Proof.** Restrict projection to $\operatorname{Ran}B$. Its kernel is the intersection with $N$, so rank-nullity gives the identity. Intersections of subspaces of dimensions $b,d$ in $\mathbb R^m$ range from $\max(0,b+d-m)$ to $\min(b,d)$; coordinate subspaces and their rotations attain both endpoints. For the error bound, apply the norm of $B_v^\dagger$ after nuisance projection. QED.

This is a first-order statement. The map $x\mapsto x^3$ is locally injective at zero despite zero derivative. Conversely, $F(u,h)=h-u^2$ has $A=0,B=1$ at zero but admits the zero-data curve $h=u^2$. Constant-rank regularity of the nuisance chart, or a full-rank derivative on a reduced identifiable parameter chart, is needed for a smooth nonlinear local inverse. No global basin follows from a positive singular value.

## 6. Rank information, neutral admission, and bias

### Theorem 5: exact nested-nuisance loss

Fix the linearization, physical target, noise metric, gauge, and frequency. Suppose $C_r\subset C_{r+1}$ and set $N_r=C_r+\operatorname{Ran}A$, $E_r=N_{r+1}\cap N_r^\perp$. Then

$$
J_x(r)-J_x(r+1)=B^TP_{E_r}B\succeq0.
$$

Equality holds exactly when $P_{E_r}B=0$. Direction $h$ loses positive information exactly when $P_{E_r}Bh\ne0$. The decrease is positive definite exactly when $P_{E_r}B$ is injective. Moreover,

$$
\operatorname{rank}B_v(r)-\operatorname{rank}B_v(r+1)
=\dim(\operatorname{Ran}B_v(r)\cap E_r).
$$

**Proof.** The orthogonal decomposition $N_{r+1}=N_r\oplus E_r$ gives $P_{N_r^\perp}-P_{N_{r+1}^\perp}=P_{E_r}$. Congruence gives the loss. Every directional loss is $\|P_{E_r}Bh\|^2$, proving equality and strictness. Projection from $\operatorname{Ran}B_v(r)$ onto $E_r^\perp$ has kernel its intersection with $E_r$, proving the rank-drop identity. QED.

The absolute map information $A^TP_{(C_r+\operatorname{Ran}B)^\perp}A$ also decreases by the same argument. Additional free current improves expressiveness, not information in the same fixed experiment. The relative map spectrum need not be monotone: take $A=(1,1,0)^T$, $B=(1,0,1)^T$, and successively $C_0=0$, $C_1=\operatorname{span}e_1$, $C_2=\operatorname{span}\{e_1,e_2+e_3\}$. The sole supported map retention is respectively $3/4,1,0$.

### Theorem 6: calibration-neutral admission

Let $N$ be the existing nuisance space, $V=P_{N^\perp}B$, and $Q=P_{N^\perp}K_+$ the residualized candidate-current matrix. Admitting a coefficient subspace $Z$ preserves all existing pose information if and only if

$$
V^TQZ=0.
$$

If $Q$ has $d$ independent columns, the largest real coefficient space with this property has dimension $d-\operatorname{rank}(V^TQ)\ge d-p$.

**Proof.** The new nuisance is $N\oplus\operatorname{Ran}(QZ)$. Theorem 5 says equality is equivalent to orthogonality between $\operatorname{Ran}(QZ)$ and $V$. This is precisely the displayed equation. The maximal coefficient space is the kernel of $V^TQ$ and rank-nullity gives its dimension. Injectivity of $Q$ makes coefficient dimension equal added data-space dimension. QED.

For truly complex free coefficients, arbitrary real kernels are not automatically admissible complex subspaces. Let $J_c$ be multiplication by $i$ in real coordinates and set $F=V^TQ$. The largest complex-linear admissible coefficient space is

$$
\ker F\cap\ker(FJ_c).
$$

It is $J_c$-invariant because $J_c^2=-I$, and every $J_c$-invariant subspace contained in $\ker F$ lies in this intersection. With $d$ candidate complex coefficients, its complex dimension is at least $d-p$, provided the same injectivity condition holds. It can still exclude physically necessary current corrections; therefore neutral admission must be paired with a bias certificate.

For an explicitly declared positive semidefinite current-approximation utility $T$, the best $r$-dimensional subspace within the safe kernel maximizes $\operatorname{tr}(Z^TTZ)$. It consists of the leading eigenvectors of the compressed operator $P_{\rm safe}TP_{\rm safe}$. To prove this, expand an orthonormal $Z$ in the compressed eigenbasis: the trace is $\sum_i\lambda_iw_i$ with $0\le w_i\le1$ and $\sum_iw_i=r$, maximized by placing unit weights on the largest eigenvalues. This utility is an approximation/prior utility, not an information claim.

### Theorem 7: bias-aware calibration risk

In the fixed linear model

$$
y=Na+Bh+D_rz+\varepsilon,\qquad \|z\|\le1,\quad \varepsilon\sim N(0,I),
$$

assume $B_v=P_{N^\perp}B$ has full column rank and let $\widehat h=B_v^\dagger y$. For any pose metric $M_x\succ0$,

$$
\sup_{\|z\|\le1}\mathbb E\|\widehat h-h\|_{M_x}^2
=\operatorname{tr}(M_xJ_x^{-1})
+\|M_x^{1/2}B_v^\dagger D_r\|^2.
$$

**Proof.** Nuisance projection annihilates $Na$, and $B_v^\dagger B=I$. The error is the sum of deterministic bias $B_v^\dagger D_rz$ and mean-zero noise $B_v^\dagger\varepsilon$. The noise covariance is $J_x^{-1}$ and the cross expectation vanishes. Maximizing the squared bias over the unit ball is the squared operator norm displayed. QED.

For a deterministic error ball of radius $\delta+\beta_r$, one instead has

$$
\|\widehat h-h\|_{M_x}\le
\frac{\delta+\beta_r}{\sigma_{\min}(B_vM_x^{-1/2})}.
$$

This deterministic statement can be used simultaneously over data-adaptive rank choices if the noise/model-error bounds hold uniformly. Fixed-rank Gaussian covariance formulas cannot simply be conditioned on a rank selected from the same noise realization.

### Proposition 8: state feasibility, soft filters, and moving projectors

For exact state-constrained response sets $\mathcal F_r(h)$, nonempty for every target direction, set inclusion $\mathcal F_r(h)\subseteq\mathcal F_{r+1}(h)$ implies nonincreasing profiled information by minimization over a larger set. Bare nesting of current bases is insufficient for approximate or projected state equations. For example, let $M=I$, $S=I$, $H_D=e_2$, $U_1=e_1$, and $U_2=I$. Galerkin responses are $0$ and $e_2h$, so their pose Grams increase from zero to one. Exact state feasibility fails at rank one for $h\ne0$.

A nonzero diagonal rescaling of free nuisance coordinates leaves their range unchanged. Thus a soft filter used only as a coordinate rescaling does not change pose information. A genuine ridge nuisance model gives

$$
\min_c\{\|z-Kc\|^2+\lambda\|c\|^2\}
=z^T[I-K(K^TK+\lambda I)^{-1}K^T]z.
$$

The information increases with $\lambda$ because the penalty imposes more nuisance control. It is not the Gram obtained by squaring the residual-filter matrix. In a singular direction, current bias is $\lambda/(\sigma^2+\lambda)$ times the true coefficient and noise amplification is $\sigma/(\sigma^2+\lambda)\le1/(2\sqrt\lambda)$.

For a constant-rank matrix $T$ and $P=TT^\dagger$,

$$
\dot P=(I-P)\dot T T^\dagger+(T^\dagger)^T\dot T^T(I-P).
$$

Differentiate $P^2=P$ and $P^T=P$ to obtain zero diagonal derivative blocks and transposed off-diagonal blocks. Differentiating $PT=T$ determines the off-diagonal block $(I-P)\dot TT^\dagger$, proving the formula. Consequently $\|\dot P\|\le2\|\dot T\|/\sigma_{\min}^+(T)$.

For a leading spectral projector of $H=S^*S$, the off-diagonal coefficients are $v_j^*\dot H v_i/(\lambda_i-\lambda_j)$, $i\le r<j$. A positive separating gap is required; multiplicities within the retained cluster are harmless. At $H(t)=\operatorname{diag}(1+t,1-t)$ the top-one projector jumps at zero without a rank change. A scalar hard threshold also jumps when its retained rank changes. These are different events.

The smooth filter $F_\lambda(H)=H(H+\lambda I)^{-1}$ has derivative $\lambda(H+\lambda I)^{-1}\dot H(H+\lambda I)^{-1}$, hence Lipschitz derivative bound $\|\dot H\|/\lambda$. Smoothing stabilizes coordinates but does not create information.

## 7. Shared-map acquisitions and an information budget

After removing each acquisition's own current nuisance, write its real white tangent as $(a_\ell,b_\ell)$. Stack them into $(A,B)$, sharing one map increment and one trajectory parameter. The profiled pose information is

$$
J=B^T(I-P_{\operatorname{Ran}A})B.
$$

### Theorem 9: shared-map compensation and innovation

A pose direction is hidden in the stack if and only if there is one map increment $u$ with $a_\ell u=b_\ell h$ for every acquisition. If every $a_\ell$ is injective and every individual acquisition hides every pose direction, let $T_\ell=a_\ell^\dagger b_\ell$. Stacked pose visibility is complete if and only if $[T_2-T_1;\ldots;T_L-T_1]$ is injective. In that special regime $(L-1)q\ge p$ is necessary and is sufficient for generic unconstrained compensation matrices, not universally for a physical array.

If $G=A^TA\succ0$, define $H=G^{-1}A^TB$ and append an acquisition $(a,b)$. Then

$$
J_{\rm new}=J+(b-aH)^T(I+aG^{-1}a^T)^{-1}(b-aH).
$$

Equality holds if and only if $b-aH=0$. The new hidden space is $\ker J\cap\ker(b-aH)$.

**Proof.** Zero projected data is equivalent to $Bh\in\operatorname{Ran}A$, giving the common compensation criterion. Injectivity of each $a_\ell$ makes its compensator unique, namely $T_\ell h$, which proves the difference-matrix criterion and dimension count.

For the innovation, complete the old square:

$$
\|Bh-Au\|^2=h^TJh+(u-Hh)^TG(u-Hh).
$$

Set $z=u-Hh$ and $V=b-aH$. The added term is $\|Vh-az\|^2$. Minimizing the sum of this term and $z^TGz$ gives $h^TV^T[I-a(G+a^Ta)^{-1}a^T]Vh$. Multiplying out verifies that the bracket equals $(I+aG^{-1}a^T)^{-1}$. It is positive definite, so equality and the kernel identity follow. QED.

A legitimate fixed map prior can be included as extra rows before this derivation. In the singular case, the variational/common-compensator criterion remains exact, but the inverse formula must not be applied without its hypothesis.

Individual profiling permits a different map per frame. Therefore $J_{\rm stack}\succeq\sum_\ell J_\ell$. Equality as matrices holds exactly when, for every $h$, the per-frame least-squares minimizing map sets have a common point. This follows directly by comparing the common minimum with the sum of independent minima; equality of nonnegative excess terms forces every frame to attain its own minimum.

The smallest example is $a_1=a_2=1$, $b_1=1$, $b_2=-1$: each frame has zero pose information but the stack has information two. Identical compensation matrices add no new rank. If each new frame also introduces unrestricted current nuisance spanning all its data rows, its cleaned acquisition is zero and it adds no information.

### Corollary 10: rank-acquisition information budget

At one fixed linearization, first enlarge the old-data nuisance space, incurring loss $L_r=B^TP_EB$. Then append a new acquisition, with innovation $I_a$ computed by Theorem 9 for the enlarged old model. The final pose information relative to the original model is exactly

$$
J_{\rm final}-J_{\rm original}=I_a-L_r.
$$

Thus $I_a\succeq L_r$ is a sufficient and necessary Loewner condition for that combined change not to degrade pose information. The proof is subtraction of the two exact identities. This is a local information-budget certificate, not a recovery theorem.

Greedy pose-Schur log-determinant design is not generally submodular: the two scalar acquisitions above have zero separate gain and positive joint gain, contradicting diminishing returns even after adding a positive pose-information offset inside the logarithm. A positive but sufficiently small map prior preserves the counterexample by continuity. Full joint-information log-determinant under fixed-dimensional independent observations has different guarantees; those cannot be transferred to the profiled pose criterion.

## 8. A posteriori numerical-physics certificates

### Theorem 11: certifying approximate states and derivatives

For $Mj=b$ with invertible $M$, let $\widetilde j$ be a computed current and $z_0=M\widetilde j-b$. If $\|M^{-1}\|\le\kappa$ is certified, then

$$
\|\widetilde j-j\|\le\kappa\|z_0\|,\qquad
\|WS(\widetilde j-j)\|\le\|WS\|\kappa\|z_0\|.
$$

For a real parameter direction $v$, let $t_v=D_vj$ and an approximate derivative $\widetilde t_v$ have residual

$$
z_v=M\widetilde t_v-b_v+M_v\widetilde j.
$$

Then

$$
\|\widetilde t_v-t_v\|\le\kappa(\|z_v\|+\|M_v\|\kappa\|z_0\|),
$$

and the data derivative error is at most $\|WS_v\|\kappa\|z_0\|+\|WS\|\|\widetilde t_v-t_v\|$ for fixed whitening.

**Proof.** Subtracting the state equations gives $\widetilde j-j=M^{-1}z_0$. Differentiating $Mj=b$ gives $Mt_v=b_v-M_vj$. Subtract this from the approximate derivative residual to obtain $\widetilde t_v-t_v=M^{-1}[z_v-M_v(\widetilde j-j)]$. The norm inequalities follow. QED.

One implementable sufficient certificate is $\|I-CM\|\le\eta<1$, for a declared preconditioner $C$. The Neumann series gives $M^{-1}=(CM)^{-1}C$ and $\kappa\le\|C\|/(1-\eta)$. A heuristic inverse-norm estimate is not a verified bound. Failure to certify does not prove singularity.

For approximate nuisance and pose matrices with certified errors $\epsilon_N,\epsilon_B$, suppose their exact and approximate nuisance ranks agree, the approximate smallest nonzero singular value is $s>\epsilon_N$, and $\widetilde P_N$ is the approximate projector. A conservative projector error is $2\epsilon_N/(s-\epsilon_N)$. To see this, write $P-\widetilde P=(I-\widetilde P)P-\widetilde P(I-P)$ and bound each term by the matrix perturbation times a pseudoinverse; the smallest exact nonzero singular value is at least $s-\epsilon_N$. Hence

$$
\|B_v-\widetilde B_v\|
\le\epsilon_B+\frac{2\epsilon_N}{s-\epsilon_N}\|\widetilde B\|.
$$

Subtract this bound, with parameter-metric scaling, from the approximate smallest pose singular value to obtain a lower certificate. Approximate small singular values cannot certify exact zeros without additional structural rank information. All certificates are local to the declared parameter region. State residuals at a wrong parameter do not prove the parameter is correct.

## 9. Safeguarded PRASC control algorithm

### Objective and state approximation

Use the physical, coherent objective

$$
\Phi(\chi,x)=\frac12\sum_\ell\|\sqrt2\mathcal E(W_\ell[F_\ell(\chi,x)-y_\ell])\|^2+R(\chi,x),
$$

with a fixed declared prior and a gauge-fixed parameter chart. For numerical state approximation, one may use $\widetilde j_{\ell,r}=j_{\det,L}(x;y)+U_{\ell,r}(x)c_{\ell,r}$, where $c_{\ell,r}$ minimizes the full state residual at fixed $(\chi,x)$. This is an approximate state solve, not an additional independent observation and not a free-current physical Fisher model. At a complete basis the state solution coincides with the physical state. A direct solver is the mandatory fallback when truncation cannot meet the numerical-physics certificate.

On a fixed-rank smooth chart, derivatives include both $D_xj_{\det,L}$ and $(D_xU_r)c_r+U_rD_xc_r$. Internal rotations of an SVD cluster must be handled by projectors or aligned bases. A spectral gap is required between retained and discarded clusters, not among all retained singular values.

### Pseudocode

```
Declare dispersion, raw coherent likelihood, gauges, metrics, priors,
accuracy budgets, admissible acquisitions, rank families, and stopping rules.
Initialize a finite set of coarse anchored candidates at the lowest informative band.
Do not interpret a zero-scatterer Jacobian as evidence that pose is observable.

For each fixed active frequency set:
    Rebuild S and e at the current geometry; cache pose-independent D.
    Compute a physical state or a certified approximate state and derivatives.
    Compute physical map/pose information separately from free-current envelopes.
    For each candidate rank or safe block admission:
        evaluate state/derivative error certificates;
        evaluate absolute pose/map strength and relative retention;
        evaluate numerical subspace uncertainty and truncation/model-bias bounds;
        reject unsupported directions; a saturated free-current envelope fails
        that envelope certificate, not the physical reduced-state model.
    Choose the least costly admissible approximation, or use the declared
    approximation utility on the bias-information Pareto set.
    If no rank is admissible:
        evaluate acquisitions with a shared-map Schur complement;
        use the innovation/rank-loss budget, including pair look-ahead if needed;
        otherwise report a local certificate failure or use direct inversion.
    Compute a joint damped step, including all moving-basis derivatives.
    Accept only against the physical objective or certified objective intervals.
    Reject steps that invalidate the declared chart, phase region, or error bounds.
    Handle rank, retention-threshold, and gap-closure events separately.
    Advance frequency only with a stated covered-error gate or a labeled surrogate.
Stop with supported local estimates and explicit uncalibrated/ambiguous modes,
or return a reason-specific failure. Never equate optimizer status with calibration.
```

The numerical state dimension r_num is not the rank r_free of an additional freely varying current nuisance. In the exact physical objective r_free=0, even when a large numerical basis is used. The nested-nuisance theorem controls r_free or a separately declared conservative envelope; it does not imply that solving the same physical state more accurately reduces its Fisher information. A complete numerical-state fallback is therefore tested against physical information, not an artificial saturated free-current envelope.

A useful gate includes absolute information: $\sigma_{\min}(B_vM_x^{-1/2})$ must exceed both its numerical/model uncertainty and the level required by the target error budget. Relative pose retention is a generalized eigenvalue against $B^TB$ on its support; a missing pose direction fails an all-pose gate. Map targets are declared physically in advance, or selected on separate data with explicit multiplicity/selection control. An empty or newly discarded map target set is not a passed gate. Relative ratios with an arbitrary epsilon denominator are not substitutes for these requirements.

### Theorem 12: safeguarded descent with inexact states

On a fixed active-data objective, assume $\Phi$ is bounded below, its gradient is Lipschitz on the relevant level set, and a positive definite step matrix satisfies $h_-I\preceq H_k\preceq h_+I$. Let $\widetilde g_k$ satisfy

$$
\|g_k-\widetilde g_k\|\le\kappa_g\|\widetilde g_k\|,\qquad
\kappa_g<h_-/h_+.
$$

Set $d_k=-H_k^{-1}\widetilde g_k$. Armijo backtracking against the physical objective gives descent and $\|g_k\|\to0$, provided all iterates remain in that level set. Every accumulation point is stationary. This does not imply convergence to the true scene.

**Proof.** The directional derivative obeys

$$
g_k^Td_k\le-\left(h_+^{-1}-\kappa_g h_-^{-1}\right)\|\widetilde g_k\|^2=-c\|\widetilde g_k\|^2,
$$

where $c>0$, and $\|d_k\|\le h_-^{-1}\|\widetilde g_k\|$. The descent lemma guarantees Armijo decrease for every sufficiently small step, with a uniform positive threshold depending on the Lipschitz constant, $c$, and $h_-$. Backtracking therefore has a uniform positive accepted step lower bound. Summing the resulting decreases and using the lower objective bound proves $\sum_k\|\widetilde g_k\|^2<\infty$. The relative error condition then gives $g_k\to0$. Continuity proves stationarity of accumulation points. QED.

An approximate objective requires certified error intervals small enough to preserve the same decrease test. With a fixed nonvanishing state-error floor, only an approximate-stationarity statement is justified. If a fixed-stratum least-squares residual has a zero at the target, full-column-rank Jacobian with minimum singular value at least $s$, and Jacobian Lipschitz constant $L$, an exact Gauss-Newton step satisfies $\|e_{k+1}\|\le L\|e_k\|^2/(2s)$ while remaining in that chart. This follows by integrating the Jacobian along the error segment and applying its pseudoinverse. The chart radius is additionally limited by the state resolvent and retained spectral gap.

Hysteresis alone does not imply finitely many events: a low-rank gate may demand expansion while the expanded-rank gate demands immediate contraction at the identical parameter. A finite-event implementation uses a finite monotone promotion ladder and forbids immediate reversal, or charges every reversible event a fixed decrease $\delta>0$ in one common bounded-below merit function. The latter allows at most $(\Phi_0-\inf\Phi)/\delta$ such events. Changing the frequency objective breaks a common-merit argument and must be treated stagewise.

### Complexity and failure outputs

A reduced state normal solve has dimension $r$, but forming $MU_r$, updating geometry-dependent bases, computing derivative certificates, and evaluating a full-physics acceptance objective can dominate its cost. A dense implementation may cost $O(nr^2+r^3)$ after operator products; direct factorization costs $O(n^3)$ but can be reused for multiple illuminations and derivatives. Iterative full-wave solves can reverse the comparison. Report operator applications, factorizations, reuse, SVD work, forward/adjoint solves, memory, and wall time. There is no unconditional cost theorem from $r<n$ alone.

Separate failure labels should include unresolved gauge, no informative starting band, unsupported pose directions, empty bias-information feasible set, uncertain numerical rank, unavailable resolvent certificate, phase-branch ambiguity, model-mismatch residual, and exhausted optimization budget. A zero envelope visibility score is an envelope robustness failure, not automatically a physical nonidentifiability certificate.

## 10. Fair-comparison theorem and verdict

Exact state elimination with an invertible $M$ and a complete invertible current coordinate map is a bijection of feasible models. Consequently the physical predictions, minimizers, and target Fisher information are identical to direct joint inversion. Linearized constrained Gauss-Newton also agrees after exact elimination; damping/trust-region metrics must be pulled back consistently. Different coordinate-Euclidean damping only changes preconditioning.

A minimal conditional benefit and failure pair is

$$
y_1=h+c+\varepsilon_1,\qquad y_2=\epsilon c+\varepsilon_2,
\qquad \varepsilon_i\sim N(0,\sigma^2),\quad\epsilon>0.
$$

Full least squares has $\widehat h=y_1-y_2/\epsilon$ and risk $\sigma^2(1+\epsilon^{-2})$. Freezing $c=0$ gives $\widehat h=y_1$ and risk $\sigma^2+c_*^2$, strictly better if $c_*^2<\sigma^2/\epsilon^2$. A declared bound on $c_*$ can make this a nonoracle conditional policy. Conversely, with negligible noise and nonzero $c_*$, full inversion is exact and truncation has pose bias $c_*$. Its leftover residual is only $\epsilon c_*$ and can be arbitrarily small. Adding an independent map observation gives unit map retention without repairing this bias.

The conditional benefit is regularization/finite-iteration control, not information creation. A direct solver permitted the same prior and subspace policy can reproduce it. The proposed algorithm therefore has no established universal or exclusive SOM advantage. What survives is a checkable local rank/acquisition safety framework with explicit bias and numerical-physics conditions. It is not yet a demonstrated TGRS imaging method or TAP calibration method.

## 11. Originality firewall and verified SOM connection

Classical SOM determines stable current coefficients from the data and optimizes the complementary ambiguous current. In Chen's 2018 author monograph, the coefficients are given in Eq. (6.57), pp. 150-152; the independent-current objective is Eq. (6.64), p. 153; and the complementary-projector implementation is Eq. (6.69), p. 161. These imply that the deterministic cutoff $L$ is not the free nuisance rank.

The same monograph defines $V_D^+$ using the right singular vectors of $G_D$, pp. 161-163. Eq. (6.70) presents an intersection model. Eq. (6.71), p. 163, explicitly distinguishes its practical sequential projection $P_{S^-}V_D^+$ from the exact intersection. The original article's low-dimensional-to-high-dimensional continuation is prior art; the monograph also explicitly describes it on p. 163 and the NFFT continuation on p. 168.

A generic exact three-way intersection decomposition must not be inferred. In $\mathbb R^2$, take $S^+=\operatorname{span}e_1$, $S^-=\operatorname{span}e_2$, and $D^\pm=\operatorname{span}(e_1\pm e_2)$. Both intersections of $S^-$ with $D^\pm$ vanish, so their sum with $S^+$ is not the whole space. Meanwhile $P_{S^-}D^+=\operatorname{span}e_2$ is nonzero. Exact intersection formulas need additional commuting/reducing-space assumptions. This is a mathematical qualification of the interpretation, not a claim that the published numerical TSOM implementation is invalid.

The project state-defect operator $D_U=H_D+MUC_U$ maps real pose increments into state-current residuals. The domain propagation operator $G_D$ maps complex current into internal field. They are not the same operator or the same SVD problem. The P-side belongs to a pushed-forward joint data tangent, not a third peer current-space fold.

Phaseless SOM retains the deterministic/ambiguous partition. Chen's monograph, Eqs. (8.28)-(8.30), pp. 221-222, shows that it uses total-field intensity and modifies recovery of the major coefficients. It is not an unconstrained phaseless current inversion.

### Selected primary references and access boundary

- X. Chen, *Computational Methods for Electromagnetic Inverse Scattering*, 2018. The user-supplied author monograph was read at equation level in the cited sections. https://www.wiley-vch.de/en/areas-interest/engineering/electrical-electronics-engineering-10ee/electromagnetic-theory-10ee4/computational-methods-for-electromagnetic-inverse-scattering-978-1-119-31198-0
- X. Chen, "Subspace-Based Optimization Method for Solving Inverse-Scattering Problems," TGRS 48, 42-49, 2010. https://doi.org/10.1109/TGRS.2009.2025122 . Publisher metadata/abstract screened; equation attribution above is to the monograph, not an unaccessed original equation number.
- Y. Zhong and X. Chen, "Twofold Subspace-Based Optimization Method for Solving Inverse Scattering Problems," Inverse Problems 25, 085003, 2009. https://doi.org/10.1088/0266-5611/25/8/085003 . Original PDF was indexed but direct retrieval failed; original equation numbers were not certified. The author's monograph provides the verified implementation account.
- L. Pan, Y. Zhong, X. Chen, and S. P. Yeo, "Subspace-Based Optimization Method for Inverse Scattering Problems Utilizing Phaseless Data," TGRS 49, 981-987, 2011. https://doi.org/10.1109/TGRS.2010.2070512 . Publisher/institutional abstract screened; full derivation verified through the author monograph rather than the original article.
- L. Bellomo et al., "An Improved Antenna Calibration Methodology for Microwave Diffraction Tomography in Limited-Aspect Configurations," TAP, 2014. https://doi.org/10.1109/TAP.2014.2308534 . Author manuscript read; Section IV-C and Eq. (25), manuscript p. 6, explicitly treat phase-center correction. Calibration modes and incident/receiving fields already have spectral and physical constraints.
- G. Huang, R. Nammour, and W. Symes, "Full-Waveform Inversion via Source-Receiver Extension," Geophysics 82, R153-R171, 2017. https://doi.org/10.1190/geo2016-0301.1 . Primary abstract/indexed author text screened; full PDF retrieval failed. Detailed equation-level equivalence to this work remains unverified.
- L. Metivier and R. Brossier, "Receiver-Extension Strategy for Time-Domain Full-Waveform Inversion Using a Relocalization Approach." https://doi.org/10.1190/geo2020-0922.1 . Primary repository record screened; full-text retrieval blocked. It is adjacent extension/relocalization prior art, not evidence of the present guarantees.
- Y. Li, K. Lee, and Y. Bresler, "Identifiability in Bilinear Inverse Problems With Applications to Subspace or Sparsity-Constrained Blind Gain and Phase Calibration," TIT, 2017. https://doi.org/10.1109/TIT.2016.2637933 . Author preprint arXiv:1501.06120 read; Section 2.1 and Definition 2.7 formalize transformation-group identifiability.
- A. J. Weiss and B. Friedlander, "Array Shape Calibration Using Sources in Unknown Locations--A Maximum Likelihood Approach," IEEE TASSP 37, 1958-1966, **1989**. https://ieeexplore.ieee.org/document/45542 . Original pages 1958-1959 visually verified, including joint least squares Eqs. (6)-(8) and missing absolute orientation without reference geometry.
- M. Cetin, O. Onhon, and S. Samadi, "Handling Phase in Sparse Reconstruction for SAR: Imaging, Autofocusing, and Moving Targets," EUSAR 2012. https://sites.rochester.edu/sdis/wp-content/uploads/2023/03/cetin_EUSAR12.pdf . Section 3.1, Eqs. (4)-(5), author manuscript p. 3, verified joint image/phase-error objective. Such a joint objective is not new here.
- Z. Idriss and R. G. Raj, "Data-Driven Calibration Technique for Quantitative Radar Imaging," arXiv:2503.07316, 2025. https://arxiv.org/abs/2503.07316 . Full HTML read; Section III, Eqs. (7)-(8), already combine multifrequency SOM and unknown complex transmitter calibration factors. These are gain/phase factors, not the present pose-rank safety theorem.
- C. Ye, R. Zhang, W. Wu, and B. Shim, "Self-Calibration DOA Estimation for Movable Antenna Systems with Antenna Position Errors," arXiv:2605.23140, 2026. https://arxiv.org/abs/2605.23140 . Version 2 HTML read; Sections II-III use an anchored steering/noise-subspace model and alternating DOA/position-error estimation. This is array/MUSIC self-calibration, not current-space SOM.
- Q. Hu, B. Zhang, and H. Zhang, "Convergent and Efficient Iteratively Regularized Contrast Source Inversion-Type Methods for Inverse Medium Scattering Problems," arXiv:2512.10260v3, 29 July 2026. https://arxiv.org/abs/2512.10260 . Full HTML read; Theorem 9 and Eq. (52) establish convergence to an epsilon-stationary point for their regularized scheme. Generic SOM convergence is not unoccupied novelty territory.
- NIST DLMF, Eqs. 10.17.5 and 10.17.11. https://dlmf.nist.gov/10.17 . Primary mathematical reference for the Hankel function and derivative asymptotics.

The bibliography is a bounded screen, not an exhaustive novelty certificate. Some important neighboring original PDFs could not be read. No assertion of firstness is justified. No new experimental result appears in this document.

## 12. Limitations and conclusion

The proved results are finite-dimensional, local, and conditional on a declared statistical and physical model. A relative spectral gate does not control truncation bias. A local inverse-information matrix does not certify the phase branch. A soft coordinate filter does not reduce a free nuisance range. A state constraint does not constitute a second observation of the same data, and eliminating an exact state returns the original physical inverse problem.

The strongest retained contribution is an information-budget framework: current admissions have an explicit loss, acquisitions have an explicit shared-map innovation, and their combination can be accepted only after bias and numerical-physics checks. This is scientifically more precise than a new name for joint inversion. Whether it becomes a useful SOM method depends on obtaining informative, inexpensive certificates and demonstrating a predeclared accuracy, basin, or computational advantage against direct joint inversion with identical data, priors, initialization, and tuning resources. That advantage remains open.

## Supplement A. A physical-state counterexample to envelope nonidentifiability

Consider the scalar discrete state $j=a(1+d j)$ with real material $a$, fixed complex $d$, and $1-ad\ne0$. Let the coherent observation be $F(a,x)=e^{ikx}j(a)$, $k\ne0$. This is the exact rational state response $j(a)=a/(1-ad)$, not a truncation of its multiple-interaction series. Its real-parameter derivatives are

$$
A=\frac{e^{ikx}}{(1-ad)^2},\qquad B=ik e^{ikx}\frac{a}{1-ad}.
$$

The determinant of the two-column real embedding is

$$
\det[\mathcal E(A),\mathcal E(B)]
=\operatorname{Im}(\overline A B)
=\frac{k a(1-a\operatorname{Re}d)}{|1-ad|^4}.
$$

It is nonzero whenever $a\ne0$ and $1-a\operatorname{Re}d\ne0$. The physical map and pose are then locally identifiable on a phase branch. However, adding an unrestricted complex current correction gives the one-column complex block $K=e^{ikx}$. Its realification spans all of $\mathbb R^2$, so the envelope map and pose information are both zero. Therefore envelope nonidentifiability is not physical nonidentifiability. This is a typed discrete counterexample, not a claim that a one-voxel model represents all Maxwell scattering.

## Supplement B. Rank uncertainty and ordinary versus certified gaps

A small computed singular value is not an exact null certificate. A structural algebraic identity, such as zero nuisance codimension, can certify a null; floating point thresholding normally certifies only detectability relative to the specified uncertainty scale. Let $\widetilde H=H+E$ be Hermitian with $\|E\|\le\epsilon$ and let the retained gap of $H$ be $\delta=\lambda_r(H)-\lambda_{r+1}(H)>\epsilon$. For old discarded and new retained eigenvectors, $X=V_-^*\widetilde V_+$ satisfies a Sylvester equation whose spectral intervals are separated by at least $\delta-\epsilon$. Writing its solution as an exponential integral gives $\|X\|\le\epsilon/(\delta-\epsilon)$. Combining the two off-diagonal projector blocks gives the conservative bound $\|P-\widetilde P\|\le2\epsilon/(\delta-\epsilon)$. This justifies an absolute gap gate tied to operator uncertainty, rather than a scale-free visual knee alone.

## Supplement C. No oracle rank in a perturbation certificate

The same-rank hypothesis in the projector bound is essential and is not itself certified by discarding small computed singular values. A verifiable sufficient case is a complete nuisance generator with exactly d columns, a proved absence of additional nuisance columns, and a certified smallest column singular value exceeding its operator error. Both the computed and true generators then have full column rank d. Alternatively, exact dependencies must be established structurally and removed by a proved factorization. If neither condition is available, the projector perturbation theorem is a conditional statement, not an implementable certificate of the exact free-nuisance information.

The obstruction already occurs with $N(t)=t e_1$ and $B=e_1$: at $t=0$ the pose information is one; at every nonzero $t$ it is zero. Thus arbitrarily small matrix error can destroy positive free-nuisance information at a rank-deficient generator. Numerical resolvable rank is not algebraic rank for a model with unbounded nuisance coefficients. A declared positive nuisance prior, a bounded correction/trust-region model, or a structurally fixed-rank chart can remove this discontinuity, but each is an additional assumption. The algorithm must label a result without such support as a local diagnostic or an unavailable certificate, rather than insert an oracle rank.


## Supplement D. Three algorithmic ranks and an implementable state chart

Keep the original data-determined cutoff L_det, the numerical state approximation dimension r_num, and the number r_free of freely varying correction coordinates distinct. A projected TSOM candidate family can be formed from the columns of P_(S-) V_(D,+), followed by rank-revealing orthogonalization in the declared current metric. Its actual dimension need not equal either the domain cutoff or the formal intersection dimension. Nestedness is checked at one fixed operator and cutoff, not presumed across geometry updates.

On a fixed chart, write R=M U, d=b-M j_det, and choose a fixed positive definite numerical residual metric M_s. If R has full column rank, the least-squares state coefficient and its directional derivative are

$$
c=(R^*M_sR)^{-1}R^*M_s d,\qquad
Dc=(R^*M_sR)^{-1}\{D(R^*M_s d)-D(R^*M_sR)c\}.
$$

Together with Dj_det, DU, and DM, this supplies an implementable full derivative of the approximate state. If the metric varies, its derivative must also be included. These normal equations are a numerical state solve, not an additional stochastic observation. Because the fixed data-derived major current includes measurement noise, a purely complementary ambiguous basis need not converge to the exact physical state as its dimension grows. The algorithm must then admit declared residual-current corrections to the major subspace, as a numerical approximation device, or fall back to an exact physical state solve. It must not silently claim that the original complementary-only SOM chart can represent every noisy physical state.

A truly free model-error current may instead be declared through F_aug=F_phys+S U c, equivalently M j=b+M U c. This is a different statistical model with an explicit state defect, correction bound or prior. The same free correction and prior must be offered to the direct baseline. A numerical state chart, a physical current, and a model-error current are not interchangeable.
