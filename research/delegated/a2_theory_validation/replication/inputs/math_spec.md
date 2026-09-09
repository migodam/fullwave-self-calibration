# A2 PRASC/SOM Mathematical Specification (implementation extract)

Extracted 2026-09-05 from `inputs/A2_PRASC_SOM_THEOREM_PACKAGE_EN.md` (EN theorem supplement), with seed ranges from `inputs/A2_PRASC_SOM_VALIDATION_PROTOCOL.md`. Statements and displayed formulas are preserved verbatim or near-verbatim so implementation can be checked against the source. The EN supplement labels the second numbered statement "Proposition 2"; the requested section headings below follow the parent numbering (Theorem 2).

## Section 2 conventions (model, metrics, and gauge)

Shared material parameter $\chi\in\mathbb R^q$, shared trajectory parameter $x\in\mathbb R^p$, frequency-dependent contrast $\chi_\ell=f_\ell(\chi)\in\mathbb C^n$; each acquisition has its own current $j_\ell\in\mathbb C^n$:

$$
y_\ell=S_\ell(x)j_\ell+\varepsilon_\ell,\qquad
M_\ell j_\ell=b_\ell,\qquad
M_\ell=I-D_{\chi_\ell}D_\ell,\quad b_\ell=D_{\chi_\ell}e_\ell(x).
$$

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

## Theorem 1 (matched observation / efficient information)

Assume a dominated statistical model differentiable in quadratic mean, square-integrable mean-zero scores, and a parameter-independent channel from $Y$ to $Z$. The target and nuisance models, including gauge restrictions, are the same in both experiments. Then

$$
s_Z=\mathbb E[s_Y\mid Z],\qquad J_Z\preceq J_Y,
\qquad J^{\rm eff}_{Z}\preceq J^{\rm eff}_{Y}.
$$

The last inequality holds for either map or pose as target, including singular nuisance information. In a target direction $v$, equality in the efficient inequality holds if and only if the coherent efficient score in that direction is measurable with respect to $Z$.

Counterexamples and boundaries: if $Y=a e^{ix}+\varepsilon$ with proper circular Gaussian noise and known $a>0$, $|Y|$ is independent of $x$, whereas coherent pose information is $2a^2/\sigma^2$. If an unknown nuisance phase $\eta$ changes the mean to $ae^{i(x+\eta)}$, the efficient information about $x$ is zero in both experiments. Positive equality also occurs for a zero-mean circular Gaussian scale family: its score is a function of $|Y|^2$. An independently engineered intensity sensor with a different noise variance is not ordered by this theorem.

The raw experiment matters. If the incident field depends on unknown pose, $Z=|Y^{\rm sca}+e_{\rm rx}(x)|^2$ is not a parameter-independent channel from scattered-field data. Use raw coherent total-field data as the parent experiment, or explicitly include the reference measurement and its noise. A known reference can retain sensitivity to the phase of the scattered field: for real nonzero $a$ and $c$, $|c+a\exp(ix)|^2$ has derivative $-2ac\sin(x)$. Thus a phase-only scattered-field perturbation is not necessarily a phase-only perturbation of the total field.

## Theorem 2 (phase-basin strong convexity / phase gate)

For the resolved single-path model

$$
f(h)=\sum_iw_i[1-\cos(k\ell_i^Th)],\qquad w_i>0,
$$

on a convex region with $|k\ell_i^Th|\le\gamma<\pi/2$ for every $i$,

$$
\nabla^2 f(h)\succeq k^2\cos\gamma\,L^TWL.
$$

If $L$ has full column rank, $f$ is strongly convex on this region. An individual path repeats after displacement $2\pi/(k|\ell_i^Tv|)$ along direction $v$. Multiple paths may remove exact common aliases but can also create destructive interference and extra minima. This is not a full-wave global basin theorem.

If a genuine error confidence set obeys $h^T\Sigma_x^{-1}h\le c^2$ and path Taylor remainders are bounded by $C_{\rm path}\|h\|^2/2$, a sufficient phase gate is

$$
k_{\rm next}\left[c\max_i\sqrt{\ell_i^T\Sigma_x\ell_i}
+\frac{C_{\rm path}c^2\lambda_{\max}(\Sigma_x)}2\right]
+\delta_{\phi,\rm model}\le\gamma<\pi/2.
$$

This follows from Cauchy-Schwarz in the ellipsoid and the Taylor bound. The square root of $\lambda_{\max}(L\Sigma_xL^T)$ is a conservative replacement for the maximum row standard deviation. For a genuinely centered Gaussian error, $c^2$ may be the appropriate chi-square quantile. An inverse Gauss-Newton matrix is not generally a covariance upper bound; using it here without a coverage argument produces a surrogate, not a certificate. A wrong phase branch can have zero residual and arbitrarily small local covariance.

## Theorem 3 and pose-prior paragraph (dual tangent spectra)

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

The generalized map eigenvalues are therefore squared sines of principal angles, but generalized eigenvalues on $\ker K_0$ are undefined rather than assigned zero or one.

A pose prior $\Lambda_x\succeq0$ yields

$$
K_{e,\Lambda}=K_0-A_c^TB_c(B_c^TB_c+\Lambda_x)^\dagger B_c^TA_c.
$$

Indeed, its quadratic form is the minimum of $\|A_cu-B_ch\|^2+\|\Lambda_x^{1/2}h\|^2$. This proves $K_e\preceq K_{e,\Lambda}\preceq K_0$ and monotonicity in the prior. Exact angle geometry is recovered in the augmented space using $[A_c;0]$ and $[B_c;\Lambda_x^{1/2}]$, not by pretending that a nonprojector in the original data space is a projector. With both map and pose priors, use columns $[A_c;\Lambda_\chi^{1/2};0]$ and $[B_c;0;\Lambda_x^{1/2}]$ and normalize by the prior-inclusive diagonal blocks. Deterministic regularization curvature is not independent data Fisher information.

## Proposition 4 (visibility / dimensions)

Let $N=\operatorname{Ran}[K_r,A]$, $d=\dim N$, and $b=\operatorname{rank}B$. Then

$$
B_v=P_{N^\perp}B,\qquad
\dim\ker B_v=p-b+\dim(\operatorname{Ran}B\cap N),
$$
$$
\max(p-b,p+d-m)\le\dim\ker B_v\le p-b+\min(b,d).
$$

These bounds are sharp over matrices with the prescribed dimensions. Pose is identifiable in the linearized experiment exactly when $\operatorname{rank}B_v=p$. On the observable support its least-squares error is bounded by the data/model-error norm divided by the smallest positive singular value.

This is a first-order statement. The map $x\mapsto x^3$ is locally injective at zero despite zero derivative. Conversely, $F(u,h)=h-u^2$ has $A=0,B=1$ at zero but admits the zero-data curve $h=u^2$. Constant-rank regularity of the nuisance chart, or a full-rank derivative on a reduced identifiable parameter chart, is needed for a smooth nonlinear local inverse. No global basin follows from a positive singular value.

## Theorem 5 (nested-nuisance loss)

Fix the linearization, physical target, noise metric, gauge, and frequency. Suppose $C_r\subset C_{r+1}$ and set $N_r=C_r+\operatorname{Ran}A$, $E_r=N_{r+1}\cap N_r^\perp$. Then

$$
J_x(r)-J_x(r+1)=B^TP_{E_r}B\succeq0.
$$

Equality holds exactly when $P_{E_r}B=0$. Direction $h$ loses positive information exactly when $P_{E_r}Bh\ne0$. The decrease is positive definite exactly when $P_{E_r}B$ is injective. Moreover,

$$
\operatorname{rank}B_v(r)-\operatorname{rank}B_v(r+1)
=\dim(\operatorname{Ran}B_v(r)\cap E_r).
$$

The absolute map information $A^TP_{(C_r+\operatorname{Ran}B)^\perp}A$ also decreases by the same argument. Additional free current improves expressiveness, not information in the same fixed experiment. The relative map spectrum need not be monotone: take $A=(1,1,0)^T$, $B=(1,0,1)^T$, and successively $C_0=0$, $C_1=\operatorname{span}e_1$, $C_2=\operatorname{span}\{e_1,e_2+e_3\}$. The sole supported map retention is respectively $3/4,1,0$.

## Theorem 6 (neutral admission, including complex-safe kernel)

Let $N$ be the existing nuisance space, $V=P_{N^\perp}B$, and $Q=P_{N^\perp}K_+$ the residualized candidate-current matrix. Admitting a coefficient subspace $Z$ preserves all existing pose information if and only if

$$
V^TQZ=0.
$$

If $Q$ has $d$ independent columns, the largest real coefficient space with this property has dimension $d-\operatorname{rank}(V^TQ)\ge d-p$.

For truly complex free coefficients, arbitrary real kernels are not automatically admissible complex subspaces. Let $J_c$ be multiplication by $i$ in real coordinates and set $F=V^TQ$. The largest complex-linear admissible coefficient space is

$$
\ker F\cap\ker(FJ_c).
$$

It is $J_c$-invariant because $J_c^2=-I$, and every $J_c$-invariant subspace contained in $\ker F$ lies in this intersection. With $d$ candidate complex coefficients, its complex dimension is at least $d-p$, provided the same injectivity condition holds. It can still exclude physically necessary current corrections; therefore neutral admission must be paired with a bias certificate.

For an explicitly declared positive semidefinite current-approximation utility $T$, the best $r$-dimensional subspace within the safe kernel maximizes $\operatorname{tr}(Z^TTZ)$. It consists of the leading eigenvectors of the compressed operator $P_{\rm safe}TP_{\rm safe}$. This utility is an approximation/prior utility, not an information claim.

## Theorem 7 (bias-aware risk)

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

For a deterministic error ball of radius $\delta+\beta_r$, one instead has

$$
\|\widehat h-h\|_{M_x}\le
\frac{\delta+\beta_r}{\sigma_{\min}(B_vM_x^{-1/2})}.
$$

This deterministic statement can be used simultaneously over data-adaptive rank choices if the noise/model-error bounds hold uniformly. Fixed-rank Gaussian covariance formulas cannot simply be conditioned on a rank selected from the same noise realization.

## Proposition 8 key facts (soft filters, ridge, moving projector derivative)

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

## Theorem 9 (shared-map compensation and innovation)

After removing each acquisition's own current nuisance, write its real white tangent as $(a_\ell,b_\ell)$. Stack them into $(A,B)$, sharing one map increment and one trajectory parameter. The profiled pose information is

$$
J=B^T(I-P_{\operatorname{Ran}A})B.
$$

A pose direction is hidden in the stack if and only if there is one map increment $u$ with $a_\ell u=b_\ell h$ for every acquisition. If every $a_\ell$ is injective and every individual acquisition hides every pose direction, let $T_\ell=a_\ell^\dagger b_\ell$. Stacked pose visibility is complete if and only if $[T_2-T_1;\ldots;T_L-T_1]$ is injective. In that special regime $(L-1)q\ge p$ is necessary and is sufficient for generic unconstrained compensation matrices, not universally for a physical array.

If $G=A^TA\succ0$, define $H=G^{-1}A^TB$ and append an acquisition $(a,b)$. Then

$$
J_{\rm new}=J+(b-aH)^T(I+aG^{-1}a^T)^{-1}(b-aH).
$$

Equality holds if and only if $b-aH=0$. The new hidden space is $\ker J\cap\ker(b-aH)$.

A legitimate fixed map prior can be included as extra rows before this derivation. In the singular case, the variational/common-compensator criterion remains exact, but the inverse formula must not be applied without its hypothesis.

Individual profiling permits a different map per frame. Therefore $J_{\rm stack}\succeq\sum_\ell J_\ell$. Equality as matrices holds exactly when, for every $h$, the per-frame least-squares minimizing map sets have a common point. This follows directly by comparing the common minimum with the sum of independent minima; equality of nonnegative excess terms forces every frame to attain its own minimum.

The smallest example is $a_1=a_2=1$, $b_1=1$, $b_2=-1$: each frame has zero pose information but the stack has information two. Identical compensation matrices add no new rank. If each new frame also introduces unrestricted current nuisance spanning all its data rows, its cleaned acquisition is zero and it adds no information.

## Corollary 10 (rank-acquisition budget)

At one fixed linearization, first enlarge the old-data nuisance space, incurring loss $L_r=B^TP_EB$. Then append a new acquisition, with innovation $I_a$ computed by Theorem 9 for the enlarged old model. The final pose information relative to the original model is exactly

$$
J_{\rm final}-J_{\rm original}=I_a-L_r.
$$

Thus $I_a\succeq L_r$ is a sufficient and necessary Loewner condition for that combined change not to degrade pose information. The proof is subtraction of the two exact identities. This is a local information-budget certificate, not a recovery theorem.

Greedy pose-Schur log-determinant design is not generally submodular: the two scalar acquisitions above have zero separate gain and positive joint gain, contradicting diminishing returns even after adding a positive pose-information offset inside the logarithm. A positive but sufficiently small map prior preserves the counterexample by continuity. Full joint-information log-determinant under fixed-dimensional independent observations has different guarantees; those cannot be transferred to the profiled pose criterion.

## Theorem 11 (a posteriori certificates, projector error bound)

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

One implementable sufficient certificate is $\|I-CM\|\le\eta<1$, for a declared preconditioner $C$. The Neumann series gives $M^{-1}=(CM)^{-1}C$ and $\kappa\le\|C\|/(1-\eta)$. A heuristic inverse-norm estimate is not a verified bound. Failure to certify does not prove singularity.

For approximate nuisance and pose matrices with certified errors $\epsilon_N,\epsilon_B$, suppose their exact and approximate nuisance ranks agree, the approximate smallest nonzero singular value is $s>\epsilon_N$, and $\widetilde P_N$ is the approximate projector. A conservative projector error is $2\epsilon_N/(s-\epsilon_N)$. To see this, write $P-\widetilde P=(I-\widetilde P)P-\widetilde P(I-P)$ and bound each term by the matrix perturbation times a pseudoinverse; the smallest exact nonzero singular value is at least $s-\epsilon_N$. Hence

$$
\|B_v-\widetilde B_v\|
\le\epsilon_B+\frac{2\epsilon_N}{s-\epsilon_N}\|\widetilde B\|.
$$

Subtract this bound, with parameter-metric scaling, from the approximate smallest pose singular value to obtain a lower certificate. Approximate small singular values cannot certify exact zeros without additional structural rank information. All certificates are local to the declared parameter region. State residuals at a wrong parameter do not prove the parameter is correct.

## Theorem 12 (safeguarded descent: inequality part)

On a fixed active-data objective, assume $\Phi$ is bounded below, its gradient is Lipschitz on the relevant level set, and a positive definite step matrix satisfies $h_-I\preceq H_k\preceq h_+I$. Let $\widetilde g_k$ satisfy

$$
\|g_k-\widetilde g_k\|\le\kappa_g\|\widetilde g_k\|,\qquad
\kappa_g<h_-/h_+.
$$

Set $d_k=-H_k^{-1}\widetilde g_k$. The directional derivative obeys

$$
g_k^Td_k\le-\left(h_+^{-1}-\kappa_g h_-^{-1}\right)\|\widetilde g_k\|^2=-c\|\widetilde g_k\|^2,
$$

where $c>0$, and $\|d_k\|\le h_-^{-1}\|\widetilde g_k\|$.

## Section 10 minimal conditional benefit/failure pair

Exact state elimination with an invertible $M$ and a complete invertible current coordinate map is a bijection of feasible models. A minimal conditional benefit and failure pair is

$$
y_1=h+c+\varepsilon_1,\qquad y_2=\epsilon c+\varepsilon_2,
\qquad \varepsilon_i\sim N(0,\sigma^2),\quad\epsilon>0.
$$

Full least squares has $\widehat h=y_1-y_2/\epsilon$ and risk $\sigma^2(1+\epsilon^{-2})$. Freezing $c=0$ gives $\widehat h=y_1$ and risk $\sigma^2+c_*^2$, strictly better if $c_*^2<\sigma^2/\epsilon^2$. A declared bound on $c_*$ can make this a nonoracle conditional policy. Conversely, with negligible noise and nonzero $c_*$, full inversion is exact and truncation has pose bias $c_*$. Its leftover residual is only $\epsilon c_*$ and can be arbitrarily small. Adding an independent map observation gives unit map retention without repairing this bias.

## Seed ranges mentioned in the validation protocol

- E1: 101-110
- E2: 201-212
- E3: 301-310
- E4: 1001-1020 (final scene/noise), with 1-10 reserved for tuning/cost calibration
- E5: 401-412
