# Phase-Preserving, Rank-Adaptive Self-Calibration for Full-Wave Subspace Optimization under Unknown Array Geometry

> Superseded on 6 September 2026 by [the A2-based physical-state and experimental revision](a2_research/PAPER_DRAFT_A2.md). This historical draft is retained for traceability; its rank/basin claims must not be used without the A2 corrections.

**Working theory-and-method draft — 5 September 2026**  
**Project name:** TriSpace SOM-SLAM  
**Provisional method name:** phase-preserving rank-adaptive self-calibrating SOM (PRASC-SOM)  
**Scientific status:** several finite-dimensional propositions are proved below; the
central algorithm, nonlinear basin advantage, and high-fidelity performance remain
to be validated. This document is a corrected manuscript foundation, not a
submission-ready claim of experimental success.

## Abstract

Coherent subspace optimization methods (SOM) exploit the complex scattered field
and can recover information that is weakened or removed by magnitude-only
measurements. The same phase that makes coherent full-wave inversion informative
also makes it highly sensitive to antenna-position error, especially as frequency
increases. Existing phaseless SOM methods address missing or unreliable phase, but
discarding phase is not the only possible response when its corruption is generated
mainly by a low-dimensional array pose. We formulate a phase-preserving alternative:
estimate that pose jointly with the scattering object while controlling how much
retained-current freedom SOM is allowed to use.

At a fixed full-wave linearization, we separate three spectra with different
meanings. The singular spectrum of the sensing operator determines which current
modes are data-supported. A pose-side projected Jacobian,
\(B_{\mathrm{vis}}=P_{\mathcal N_r^\perp}B\), measures which pose directions survive
retained-current and map nuisance. A map-side generalized spectrum,
\(K_{\mathrm{eff}}v_i=\rho_iK_0v_i\), measures the fraction of known-pose map
information that survives unknown pose after retained-current nuisance is removed.
We prove that, for nested retained nuisance spaces at the same linearization, the
pose information matrix \(J_x(r)=B^TP_{\mathcal N_r^\perp}B\) is nonincreasing in
the Loewner order as retained rank increases. Thus higher rank can recover richer
current content while weakening—and, under nuisance saturation, eliminating—local
pose information. We also state a matched-observation information-contraction result for
phaseless data and derive the familiar high-frequency tradeoff in a local
large-\(kr\) path model: pose sensitivity grows as \(O(k)\), and local Fisher
information as \(O(k^2)\) at fixed coherent SNR, while phase periodicity shrinks
the unambiguous displacement scale as \(O(k^{-1})\).

These facts motivate PRASC-SOM: begin with low frequency and conservative current
rank to obtain a broad calibration basin; retain complex phase; update pose and map
on fixed-rank strata; increase frequency and rank only when current support, pose
visibility, and map-survival margins are simultaneously acceptable; and redesign
the acquisition when the stacked quotient pose Jacobian is deficient. This is not
claimed to be novel merely because it estimates pose and image jointly. The
candidate contribution is the SOM-specific use of retained-rank, dual
\(\rho/B_{\mathrm{vis}}\) spectra, and phase-error-aware frequency continuation as a
calibration control law. Prior deterministic experiments support the tangent
semantics and rank-hiding obstruction but do not yet validate the proposed method.

**Keywords:** electromagnetic inverse scattering; subspace optimization method;
self-calibration; coherent phase; unknown antenna geometry; rank adaptation;
frequency continuation; local identifiability.

## 1. Introduction

Full-wave inverse scattering depends on geometry as well as material. Antenna
coordinates determine propagation distance, incident phase, receiver sampling,
polarization response, and sometimes the current induced inside the object. A small
coordinate error can therefore produce a large coherent residual at medium and high
frequency. Standard SOM, which builds reduced current coordinates from the sensing
operator, inherits a second difficulty: when geometry changes, the operator and its
retained singular subspace change too.

One response is to avoid phase. Phaseless SOM is established prior art and still
uses spectral information to divide the contrast source into data-determined and
ambiguous parts; it is incorrect to say that phaseless SOM has no current-space
constraint. The relevant limitation is more precise. Under a matched generative
model, intensity or magnitude is a many-to-one observation of the coherent complex
field. Phase-only parameter signatures can disappear, and the transformed
observation cannot contain more Fisher information than its coherent parent. If
phase corruption is caused mainly by a low-dimensional pose, discarding the entire
phase channel may sacrifice exactly the information that can estimate that pose.

This paper therefore starts from a different principle:

> Do not remove phase merely because geometry corrupts it; estimate the
> low-dimensional geometry that corrupts phase, and control the SOM current space so
> that it cannot absorb all pose evidence.

The principle alone is not novel. Joint image-and-geometry inversion,
self-calibration, autofocus, blind calibration, and source/receiver extension all
precede this work. The scientific question is narrower and harder:

> Does SOM provide a specific spectral mechanism that makes coherent
> self-calibration more controllable than an undifferentiated joint full-wave solve?

Our proposed answer has three parts.

1. **Dual observability spectra.** The current spectrum, the pose-visible spectrum
   of \(B_{\mathrm{vis}}\), and the map-survival spectrum \(\rho\) answer three
   different questions and must be monitored together.
2. **A rank--calibration theorem.** At a fixed linearization, nested retained-current
   spaces produce nested nuisance spaces, so pose Fisher/Gram information is
   monotonically nonincreasing with SOM rank. Rank is therefore a scientific control
   variable, not only a reconstruction truncation parameter.
3. **Dual continuation.** Low frequency provides a wider phase branch, while low
   rank reduces current nuisance and pose hiding. After pose uncertainty falls,
   higher frequency and higher rank can be unlocked to recover finer scattering
   structure.

The resulting method is phase preserving: it uses the complex field throughout.
The phrase “feature fusion” means a physically shared, whitened, multifrequency and
multi-acquisition objective—not an informal concatenation of incompatible
features. Frequencies share material parameters through an explicit dispersion
model; frames share pose or trajectory variables; and current nuisance may remain
acquisition-specific.

### 1.1 Contributions and epistemic status

The manuscript makes or targets the following claims.

1. **[Proved, finite dimensional]** A retained-current representation splits a
   receiver-geometry tangent into a minimum-norm SOM lift and an orthogonal leakage.
   The unrestricted full-row-rank lift is vacuous as an identifiability claim.
2. **[Proved, fixed linearization and nested effective nuisance]** Retained-rank
   growth makes the pose information matrix nonincreasing in the Loewner order;
   equality and strict-loss conditions are explicit.
3. **[Proved under a matched observation channel]** A phaseless transformation
   cannot increase the joint Fisher information of the coherent observation;
   equality requires the coherent score to be recoverable from the phaseless data.
4. **[Proved, finite dimensional]** The generalized eigenvalues \(\rho_i\) are
   squared sines of principal angles between current-cleaned map and pose tangent
   ranges, on the observable map support and without a pose prior.
5. **[Method proposal, open performance claim]** A rank- and frequency-adaptive
   coherent SOM algorithm uses those spectra as acceptance gates for calibration,
   imaging refinement, and acquisition adaptation.
6. **[Preliminary numerical evidence only]** Existing experiments verify tangent
   identities, rank events, nuisance saturation, and some source-diversity effects.
   They falsify a proposed state witness and a naive visible-coordinate solver. They
   do not establish an advantage for PRASC-SOM.

## 2. Prior art and the originality firewall

### 2.1 SOM, Twofold SOM, and phaseless SOM

Chen's SOM uses the singular system of the data operator to split the contrast
source into a deterministic component and an ambiguous component. Zhong and Chen's
Twofold SOM introduces an additional domain-side reduction. Pan, Zhong, Chen, and
Yeo extend SOM to phaseless data; their method still uses spectrum analysis and a
deterministic/ambiguous contrast-source partition, with the deterministic recovery
modified because phase is unavailable. Consequently:

- SOM, truncated current coordinates, and deterministic/ambiguous current
  partitions are prior art;
- “phaseless SOM loses all current constraints” is false and is not used here;
- the present question is whether retaining complex phase and adapting the SOM rank
  enables a distinct self-calibration mechanism.

The exact Twofold domain fold is not reconstructed from memory in this draft. It
must be inserted only after equation-level verification of the primary source.

### 2.2 Self-calibration and geometry extension

Antenna phase-center calibration, microwave-system calibration, joint
transmitter/image reconstruction, blind gain/phase calibration, radar autofocus,
and source/receiver-extension full-waveform inversion all establish that geometry
or instrument parameters can be estimated or extended jointly with a scene.
Therefore none of the following is claimed as new:

- adding pose variables to a nonlinear inverse problem;
- alternating between pose and material updates;
- using coherent phase to localize an antenna;
- using low-to-high frequency continuation;
- eliminating nuisance variables by a Schur complement;
- using singular values or principal angles as generic conditioning diagnostics.

### 2.3 The narrow candidate contribution

| Object | Status | Role in this work |
|---|---|---|
| coherent/phaseless inverse scattering | prior art | defines the information choice |
| SOM and Twofold SOM | prior art | supplies structured current coordinates |
| joint pose--map inversion | prior art | mandatory baseline |
| receiver/source extension | prior art | closest geometry-nuisance genus |
| \(K_{\mathrm{eff}}\), Schur complement, principal angles | textbook machinery | map-side information accounting |
| retained-rank monotonic loss of pose information | derivable linear algebra, physically specialized here | explains why SOM rank controls calibration |
| dual \(\rho/B_{\mathrm{vis}}\) gating | supported recombination, novelty candidate | protects map and pose simultaneously |
| phase-preserving rank/frequency continuation | method hypothesis, novelty candidate | converts phase sensitivity into calibration |
| state witness \(T_U\) as a solver | falsified in tested model | retained only as a negative constraint |

The paper is scientifically valuable only if the final theory or experiments show
that the SOM-specific control law yields at least one of: a larger reliable basin,
lower computational cost at matched accuracy, better pose/image accuracy at matched
cost, or a verifiable acquisition/rank policy unavailable to a static direct joint
solve. If it does none of these, the method claim fails even if all projector
identities are correct.

## 3. Typed coherent full-wave model

### 3.1 Independent-current and reduced-map equations

For acquisition/frequency index \(\ell\), let

\[
y_\ell=S_\ell(x)j_\ell+\varepsilon_\ell,\qquad
j_\ell-D_{\chi_\ell}\{e_\ell(x)+D_\ell j_\ell\}=0.
\tag{1}
\]

Here \(y_\ell\in\mathbb C^{M_\ell}\), current
\(j_\ell\in\mathbb C^n\), real material parameters
\(\chi\in\mathbb R^q\), and pose/trajectory parameters
\(x\in\mathbb R^p\). A physical multifrequency model maps shared material
parameters \(\vartheta\) to \(\chi_\ell=\chi(\omega_\ell;\vartheta)\); it must not
silently keep a frequency-scaled scattering potential fixed.

For a world-fixed domain grid and fixed homogeneous background, receiver motion
changes \(S_\ell\), transmitter motion changes \(e_\ell\), and the
domain-to-domain operator \(D_\ell\) is pose independent. With

\[
M_\ell=I-D_{\chi_\ell}D_\ell,\qquad
E_{\mathrm{tot},\ell}=e_\ell+D_\ell j_\ell,
\]

the local independent-current equations are

\[
\delta y_\ell=S_\ell\delta j_\ell+H_{S,\ell}h,
\qquad
H_{S,\ell}h=(D_xS_\ell[h])j_\ell,
\tag{2}
\]

\[
M_\ell\delta j_\ell
-D_{E_{\mathrm{tot},\ell}}\delta\chi_\ell
-H_{D,\ell}h=0,
\qquad
H_{D,\ell}h=D_{\chi_\ell}D_xe_\ell[h].
\tag{3}
\]

The receiver term \(H_{S,\ell}h\) is a re-sampling signature and can be
represented by a pseudo-current. The transmitter term
\(M_\ell^{-1}H_{D,\ell}h\), when it exists, is a physical induced-current
variation. Lifting the total pose Jacobian through \(S_\ell^\dagger\) while also
retaining (3) double counts the transmitter effect.

For hiding diagnostics, it is useful to write the augmented data-space envelope

\[
\delta y_\ell
=K_{\ell,r}\delta c_\ell
+A_{\chi,\ell}\delta\chi
+B_\ell h,
\qquad
K_{\ell,r}=S_\ell U_{\ell,r}.
\tag{4}
\]

The columns of \(U_{\ell,r}\) are retained current coordinates. Equation (4) is
not automatically the exact Jacobian of one physical parameterization. In a pure
independent-current formulation, the direct data block has
\(A_{\chi,\ell}=0\), and material enters through the state constraint (3). In a
pure reduced-state formulation, \(A_{\chi,\ell}\) and \(B_\ell\) are total
full-wave derivatives and the additional free block \(K_{\ell,r}\delta c_\ell\)
is absent unless an explicit model-error/current-correction variable is declared.
A method that retains all three blocks must state that augmented variable and its
state coupling. Otherwise it can double count physical current response. The
nuisance space used below deliberately allows all three blocks as a conservative
data-side envelope; its monotonicity theorem applies to that envelope, or to a
state-feasible model only when the resulting effective nuisance spaces are nested.

### 3.2 Whitening, realification, and gauge quotient

Let \(W_\ell^*W_\ell=\Sigma_{\varepsilon,\ell}^{-1}\). All projectors and
singular values are computed after whitening. For a complex matrix multiplying a
complex coefficient, define

\[
\mathcal R(A)=
\begin{bmatrix}
\Re A&-\Im A\\
\Im A& \Re A
\end{bmatrix}.
\tag{5}
\]

For a complex Jacobian multiplying a real vector, define

\[
\mathcal E(A)=
\begin{bmatrix}
\Re A\\
\Im A
\end{bmatrix}.
\tag{6}
\]

Thus

\[
\bar K_{\ell,r}=\mathcal R(W_\ell K_{\ell,r}),\quad
\bar A_\ell=\mathcal E(W_\ell A_{\chi,\ell}),\quad
\bar B_\ell=\mathcal E(W_\ell B_\ell).
\tag{7}
\]

This prevents complex current freedom from being confused with real map or pose
freedom.

In general, a map--platform rigid gauge is a *joint* tangent
\(\mathcal G\subset\mathbb R^q\times\mathbb R^p\), not the direct sum of an
independent map gauge and pose gauge. We therefore assume below that an external
anchor or an explicit local gauge-fixing constraint has first made the map and pose
charts separable. Only under that convention do \(Z_\chi\) and \(Z_x\) denote
orthonormal bases of the remaining map and pose coordinates, with
\(\widehat A_\ell=\bar A_\ell Z_\chi\) and
\(\widehat B_\ell=\bar B_\ell Z_x\). Without such an anchor, every block formula
must instead be pulled back to a complement of the *joint* gauge before target and
nuisance blocks are defined. Additional unanchored measurements do not remove the
joint global gauge.

## 4. The three coupled spectra

### 4.1 Current support and SOM margin

Let \(S_\ell=P_\ell\Sigma_\ell V_\ell^*\) and let
\(U_{\ell,r}=V_{\ell,1:r}\) for a hard truncated SOM basis. The quantities

\[
\eta_S(\ell,r)=\frac{\sigma_{\ell,r}}{\sigma_{\ell,1}},
\qquad
g_S(\ell,r)=
\frac{\sigma_{\ell,r}-\sigma_{\ell,r+1}}{\sigma_{\ell,1}}
\tag{8}
\]

answer whether the last retained current mode is data-supported and whether the
hard rank is separated from its complement. A small \(\eta_S\) amplifies
pseudo-current coefficients; a small gap makes the retained projector highly
sensitive to geometry. A soft spectral filter replaces a discontinuous rank event
with bias; it does not create new information. Equation (8) is used below the
chosen ambient maximum rank; at the terminal rank we set the next singular value
to zero by convention.

### 4.2 Pose-side visibility

For one stacked linearization, suppress \(\ell\) temporarily and define the
retained-current range and total pose nuisance

\[
\mathcal C_r=\operatorname{Range}(\bar K_r),
\qquad
\mathcal N_r=\mathcal C_r+\operatorname{Range}(\widehat A).
\tag{9}
\]

The pose signature that current and map cannot reproduce is

\[
B_{\mathrm{vis}}(r)=P_{\mathcal N_r^\perp}\widehat B,
\qquad
J_x(r)=B_{\mathrm{vis}}(r)^TB_{\mathrm{vis}}(r).
\tag{10}
\]

It answers:

> Which quotient pose directions remain unique after the declared current and map
> variations are allowed to explain the data?

Local quotient pose identifiability holds exactly when

\[
\operatorname{rank}B_{\mathrm{vis}}(r)=p-\dim\mathcal G_x.
\tag{11}
\]

The smallest positive singular value controls local conditioning, whereas the rank
alone only detects exact first-order loss.

### 4.3 Map-side survival and the generalized spectrum \(\rho\)

First remove retained-current nuisance:

\[
\Pi_r=P_{\mathcal C_r^\perp},\qquad
\widetilde A_r=\Pi_r\widehat A,\qquad
\widetilde B_r=\Pi_r\widehat B.
\tag{12}
\]

The known-pose map information and unknown-pose effective map information are

\[
K_0(r)=\widetilde A_r^T\widetilde A_r,
\tag{13}
\]

\[
K_{\mathrm{eff}}(r)
=\widetilde A_r^T
\bigl(I-P_{\operatorname{Range}(\widetilde B_r)}\bigr)
\widetilde A_r.
\tag{14}
\]

With a Gaussian pose prior of precision \(\Lambda_x\succeq0\), replace (14) by the
Schur information

\[
K_{\mathrm{eff},\Lambda}(r)
=K_0(r)
-\widetilde A_r^T\widetilde B_r
\bigl(\widetilde B_r^T\widetilde B_r+\Lambda_x\bigr)^\dagger
\widetilde B_r^T\widetilde A_r.
\tag{15}
\]

On the observable support of \(K_0(r)\), define

\[
K_{\mathrm{eff}}(r)v_i=\rho_i(r)K_0(r)v_i,
\qquad 0\le \rho_i(r)\le1.
\tag{16}
\]

The interpretation is:

> \(\rho_i\) is the fraction of known-pose information in map mode \(v_i\) that
> remains after unknown pose is eliminated.

This construction is generic nuisance-elimination mathematics. It is not specific
to inverse scattering. The inverse-scattering content lies in the full-wave
operators \(A,B\), the SOM current range removed before them, and the resulting
frequency/rank behavior.

**Proposition 1 (principal-angle form).** Suppose \(K_0(r)\) is positive definite
on a selected map subspace and no pose prior is used. The generalized eigenvalues
in (16) are the squared sines of the principal angles between
\(\operatorname{Range}(\widetilde A_r)\) and
\(\operatorname{Range}(\widetilde B_r)\), with unit eigenvalues for map directions
orthogonal to the pose range.

**Proof.** Take a thin QR factorization
\(\widetilde A_r=Q_AR_A\), with \(R_A\) nonsingular on the selected support. Under
the congruence \(w=R_Av\), the pencil becomes

\[
\bigl[I-Q_A^TP_{\operatorname{Range}(\widetilde B_r)}Q_A\bigr]w
=\rho w.
\]

The nonzero singular values of \(Q_A^TQ_B\), for an orthonormal basis \(Q_B\) of
the pose range, are the cosines of the principal angles. Therefore the eigenvalues
are \(1-\cos^2\theta_i=\sin^2\theta_i\), supplemented by ones when appropriate.
\(\square\)

### 4.4 Duality without conflation

\(B_{\mathrm{vis}}\) and \(\rho\) are opposite block-Schur views of the same
current-cleaned joint tangent:

- \(B_{\mathrm{vis}}\) asks how much pose information survives current and map;
- \(\rho\) asks how much map information survives current and pose.

The first eliminates current and map to assess pose; the second eliminates current
and pose to assess map. They are not the same matrix and need not have the same
dimension. Under full-rank, no-prior normalization, both block views are governed
by the same canonical correlations between the current-cleaned map and pose ranges;
their raw spectra are still different. The exact singular/prior-weighted extension
remains an open theorem target. The complete rank decision therefore requires three
gates:

\[
\boxed{
\text{current support }\sigma(S),\quad
\text{pose survival }\sigma(B_{\mathrm{vis}}),\quad
\text{map survival }\rho.
}
\tag{17}
\]

Importantly, \(\rho_i(r)\) is a ratio whose numerator and denominator both change
with \(r\). It is not generally monotone in rank. Only the pose-side Loewner
monotonicity proved next is unconditional under nested nuisance spaces.

## 5. Theory supporting phase-preserving self-calibration

### 5.1 Equivalent-current lift: useful coordinate, insufficient evidence

For a retained current basis \(U_r\), define

\[
C_r=(S U_r)^\dagger H_S,\qquad
Q_r=U_rC_r,\qquad
R_r=(I-P_{\operatorname{Range}(S U_r)})H_S.
\tag{18}
\]

Then

\[
H_S=S U_rC_r+R_r,\qquad
(S U_r)^*R_r=0,\qquad
\operatorname{rank}C_r\le p.
\tag{19}
\]

For each pose direction, \(Q_rh\) is the minimum-norm retained pseudo-current in
the Euclidean form shown when \(U_r\) is orthonormal; a non-Euclidean current norm
requires the corresponding weighted pseudoinverse and projector. If \(S\) is full
row rank and unrestricted current is allowed, \(R=0\) for every \(H_S\). Exact
unrestricted equivalence is therefore a negative control: it demonstrates
expressive redundancy, not geometry identifiability, non-scattering information,
or improved \(K_{\mathrm{eff}}\).

### 5.2 Coherent-to-phaseless information contraction

Let the real parameter vector \(\theta\) contain map, pose, and any declared
nuisance variables. Let \(Y\) be a coherent complex observation whose dominated
model is differentiable in quadratic mean, has a parameter-independent local
support, and has square-integrable score. Let \(Z\) be generated from \(Y\) by a
\(\theta\)-independent Markov channel, including the deterministic transformation
\(Z=|Y|^2\). Assume the derivative can be passed through the channel integral.

**Proposition 2 (Fisher information data processing).** Under the conditions above,

\[
J_Z(\theta)\preceq J_Y(\theta).
\tag{20}
\]

Equality holds if and only if the coherent score is measurable with respect to
\(Z\), almost surely, in every parameter direction under comparison.

**Proof.** Let \(s_Y=\nabla_\theta\log p_\theta(Y)\). The score identity for a
parameter-independent channel gives
\(s_Z=\mathbb E_\theta[s_Y\mid Z]\). Since the score has zero mean,

\[
J_Z
=\mathbb E[s_Zs_Z^T]
=\operatorname{Var}(\mathbb E[s_Y\mid Z])
\preceq\operatorname{Var}(s_Y)
=J_Y
\]

by conditional-variance decomposition. Equality is equivalent to
\(s_Y=\mathbb E[s_Y\mid Z]\) almost surely. \(\square\)

For noiseless differentiable fields,

\[
D_\theta |y|^2[h]
=2\operatorname{Re}\{\operatorname{diag}(\bar y)D_\theta y[h]\}.
\tag{21}
\]

An elementwise phase perturbation \(D_\theta y[h]=i\alpha\odot y\), with real
\(\alpha\), is annihilated by (21). Thus the observation most sensitive to small
path-length changes in coherent phase can be invisible to first-order intensity.

Equation (20) is a matched-model statement. A direct intensity sensor with a
different noise law, dynamic range, or calibration error is not automatically
comparable to a coherent receiver followed by a deterministic transformation.
Likewise, a nuisance-eliminated pose Fisher inequality requires the same parameter
and nuisance experiment; a rigorous efficient-information version is a priority
for the final theorem package.

### 5.3 Frequency increases local precision and global ambiguity

Away from the source singularity, the outgoing two-dimensional Helmholtz Green
function has the large-\(kr\) form

\[
G_k(r)=a(k,r)e^{ikr}\{1+O((kr)^{-1})\},
\qquad |a(k,r)|=O((kr)^{-1/2}).
\tag{22}
\]

For a path-length perturbation \(\delta r\),

\[
\delta\arg G_k=k\,\delta r
+O\!\left(\frac{|\delta r|}{k r^2}\right),
\qquad
\frac{\delta |G_k|}{|G_k|}
=O\!\left(\frac{|\delta r|}{r}\right).
\tag{23}
\]

Hence, at fixed coherent SNR and within one local branch, the phase-dominated pose
Jacobian scales as \(O(k)\) and its Gram/Fisher information as \(O(k^2)\). But a
\(2\pi\) phase change occurs after a path displacement of order \(2\pi/k\).
Medium/high frequency is therefore locally excellent for pose precision and
globally dangerous for initialization. Multiple scattering, limited aperture, and
multiple paths alter constants and can create additional minima; (23) is a local
asymptotic mechanism, not a global convexity theorem.

### 5.4 Retained-rank monotonicity

Assume a fixed pose, map, frequency, whitening, and gauge quotient. Also assume
that the admissible data-side current nuisance is exactly the retained range below,
or that state elimination produces effective nuisance spaces with the same nesting.
Let the retained spaces be nested:

\[
\operatorname{Range}(U_r)\subseteq\operatorname{Range}(U_{r+1}).
\tag{24}
\]

Then \(\mathcal C_r\subseteq\mathcal C_{r+1}\) and
\(\mathcal N_r\subseteq\mathcal N_{r+1}\).

**Theorem 1 (monotone loss of pose information).**

\[
P_{\mathcal N_{r+1}^\perp}\preceq P_{\mathcal N_r^\perp},
\qquad
J_x(r+1)\preceq J_x(r).
\tag{25}
\]

Consequently,

\[
\operatorname{rank}B_{\mathrm{vis}}(r+1)
\le\operatorname{rank}B_{\mathrm{vis}}(r).
\tag{26}
\]

Let
\(\mathcal E_r=\mathcal N_{r+1}\cap\mathcal N_r^\perp\) be the newly added
nuisance component. Then

\[
J_x(r)-J_x(r+1)=\widehat B^TP_{\mathcal E_r}\widehat B.
\tag{27}
\]

Equality holds exactly when \(P_{\mathcal E_r}\widehat B=0\). A pose direction
\(h\) loses strictly positive information exactly when
\(P_{\mathcal E_r}\widehat Bh\ne0\).

**Proof.** For nested finite-dimensional subspaces,
\(P_{\mathcal N_r^\perp}-P_{\mathcal N_{r+1}^\perp}
=P_{\mathcal E_r}\succeq0\). Congruence by \(\widehat B\) gives (25) and (27).
For positive semidefinite matrices \(0\preceq J_x(r+1)\preceq J_x(r)\),
\(\ker J_x(r)\subseteq\ker J_x(r+1)\), which gives (26). The equality and
strictness statements follow from the squared norm
\(h^T[J_x(r)-J_x(r+1)]h=\|P_{\mathcal E_r}\widehat Bh\|^2\).
\(\square\)

The theorem is stronger than a dimension count: it orders the full directional
pose information. It is a conservative theorem for free retained-current nuisance.
If the state equation restricts newly added current coefficients so that effective
nuisance ranges are not nested, (25) cannot be invoked without a new derivation. It
is also not a cross-frequency or cross-iteration theorem because \(A,B,S\),
whitening, and the linearization can change. If a hard threshold changes the
retained dimension, that is a rank event. If \(\rho_i\) crosses a decision threshold
while rank stays fixed, that is an observability event, not a rank event.

**Corollary 1 (dimension obstruction).** In a real data space of dimension \(m\),

\[
\dim\ker B_{\mathrm{vis}}(r)
\ge
\max\{0,\;p_g+\dim\mathcal N_r-m\},
\tag{28}
\]

where \(p_g\) is the dimension of the pose chart remaining after joint gauge fixing
or anchoring. If \(\mathcal N_r\) fills the data space, all quotient pose
information vanishes.

### 5.5 Stacked acquisition with shared map

For \(L\) acquisitions, let retained-current coefficients be frame-specific while
map and pose are shared. After realification and whitening,

\[
\bar K_r=\operatorname{blkdiag}(\bar K_{1,r_1},\ldots,\bar K_{L,r_L}),
\quad
\bar A=
\begin{bmatrix}\widehat A_1\\ \vdots\\ \widehat A_L\end{bmatrix},
\quad
\bar B=
\begin{bmatrix}\widehat B_1\\ \vdots\\ \widehat B_L\end{bmatrix}.
\tag{29}
\]

Remove frame-specific current nuisance first:

\[
\Pi_C=P_{\operatorname{Range}(\bar K_r)^\perp},\qquad
\widetilde A=\Pi_C\bar A,\qquad
\widetilde B=\Pi_C\bar B.
\tag{30}
\]

Then the exact stacked quotient pose information is

\[
J_{x,\mathrm{stack}}
=\widetilde B^T
\bigl(I-P_{\operatorname{Range}(\widetilde A)}\bigr)
\widetilde B.
\tag{31}
\]

The stacked system is locally pose-identifiable modulo gauge if and only if
\(\operatorname{rank}J_{x,\mathrm{stack}}=p_g\). This formula preserves the shared
map coupling; summing independently map-projected per-frame matrices can be wrong
because each frame would be allowed a different map compensator.

Stacking can recover a direction hidden in every single frame if the compensating
map perturbations required by those frames are mutually inconsistent under the
shared-map constraint. Conversely, adding measurements need not help when the
current nuisance grows in the same block directions, transmitters are symmetry
equivalent, or the global gauge remains unanchored.

## 6. PRASC-SOM: the proposed self-calibration design

### 6.1 Design principle

PRASC-SOM is not a new algebraic parameterization of the same static objective.
Its candidate value is a policy that changes the admissible current rank,
frequency, and acquisition using observable spectral margins. The method follows:

\[
\boxed{
\text{calibrate with conservative wave/current complexity}
\;\longrightarrow\;
\text{unlock high-frequency/high-rank scattering information}.
}
\tag{32}
\]

At frequency \(k\) and candidate rank \(r\), let \(M_x\succ0\) be a declared pose
parameter metric and define normalized gates

\[
\eta_x(r)=
\frac{\lambda_{\min}\!\left(M_x^{-1/2}J_x(r)M_x^{-1/2}\right)}
     {\lambda_{\max}\!\left(
       M_x^{-1/2}\widehat B^T\widehat B M_x^{-1/2}\right)+\epsilon},
\tag{33}
\]

\[
\eta_\rho(r)=
\min_{v_i\in\mathcal T_k}\rho_i(r),
\tag{34}
\]

where \(\mathcal T_k\) is the map subspace scheduled for update at the current
stage. A hard-rank feasible set is

\[
\mathcal F_k=
\left\{
r:
\eta_S(r)\ge\gamma_S,\;
g_S(r)\ge\gamma_{\mathrm{gap}},\;
\eta_x(r)\ge\gamma_x,\;
\eta_\rho(r)\ge\gamma_\rho
\right\}.
\tag{35}
\]

When a soft filter is used, the gap gate is replaced by an effective-degree and
projector-Lipschitz bound. Because \(\rho(r)\) need not be monotone, the algorithm
searches the feasible set rather than assuming that all ranks below a threshold
are feasible. Hysteresis prevents rank chatter.

### 6.2 Phase-error-aware frequency gate

Let \(\Sigma_x\) be a local pose covariance or inverse regularized information
estimate. For each candidate frequency, let \(L_{\mathrm{path}}(k)\) map pose
error to the dominant path-length perturbations. Advance to \(k_{\mathrm{next}}\)
only when

\[
k_{\mathrm{next}}
\sqrt{\lambda_{\max}
\left(
L_{\mathrm{path}}\Sigma_x L_{\mathrm{path}}^T
\right)}
\le\gamma_\phi,
\tag{36}
\]

where \(\gamma_\phi<\pi\) is a conservative phase-uncertainty budget chosen for the
acquisition. Equation (36) is an operational gate, not a universal basin theorem.
It turns estimated pose uncertainty into a frequency-continuation decision.

### 6.3 Algorithm

**Algorithm 1: phase-preserving rank-adaptive self-calibrating SOM**

**Input:** coherent complex measurements; noise covariance; material-dispersion
model; low-frequency initialization; candidate rank/filter family; map and pose
gauge anchors; thresholds and hysteresis; admissible acquisition actions.

1. Initialize at the lowest informative frequency with a conservative current rank.
   Estimate coarse pose and low-order map/current modes using the complete
   full-wave residual.
2. At the current iterate, rebuild \(S,e,D\), the full-wave state, and analytic or
   verified automatic-differentiation Jacobians. Whiten and realify them.
3. For each candidate rank/filter, compute the current support and gap, nuisance
   range, \(B_{\mathrm{vis}}\), \(J_x\), \(K_0\), \(K_{\mathrm{eff}}\), and
   target \(\rho_i\).
4. Select a feasible rank that maximizes a declared imaging utility, such as
   recoverable map dimension or expected reduction, subject to all gates in (35).
   If no rank is feasible, reduce rank/frequency, strengthen a legitimate prior,
   or request a new acquisition; do not call a rank-deficient update calibrated.
5. On a fixed-rank stratum, compute a joint full-wave Gauss--Newton or
   Levenberg--Marquardt step. Eliminate current/map blocks by exact Schur
   algebra only when this preserves the declared objective and state model.
6. Accept the step using the complete coherent data residual, state residual where
   applicable, trust-region ratio, pose/map change, gauge consistency, and spectral
   stability. Optimizer status alone is not an acceptance certificate.
7. If rank changes, treat it as a rank event: rebuild the local model, transport
   projectors rather than arbitrary SVD vectors, and reset incompatible curvature
   history. If \(\rho_i\) crosses a threshold at fixed rank, treat it as an
   observability event and update the active map-mode set.
8. When (36) holds and the next frequency adds supported current modes, increase
   frequency. Increase rank only after pose and map gates remain satisfied.
9. If the stacked \(B_{\mathrm{vis}}\) is deficient or ill-conditioned, choose a
   transmitter, directional illumination, frame, aperture, frequency, trajectory,
   or anchor that improves the quotient criterion in (31). Recompute rather than
   assuming that more samples are independent.
10. Stop when full-physics residuals are noise-consistent, pose and map increments
    are stable, all claimed modes pass their information gates, and no unresolved
    gauge remains.

### 6.4 Acquisition objective

A practical local design score is a robust combination of

\[
\lambda_{\min}(J_{x,\mathrm{stack}}),\qquad
\log\det(J_{x,\mathrm{stack}}+\Lambda_x),\qquad
\min_{i\in\mathcal T}\rho_i,\qquad
\operatorname{tr}K_{\mathrm{eff}}.
\tag{37}
\]

The first two protect pose; the latter two protect map information. The design must
be evaluated after shared-nuisance elimination. Multiple copies of a symmetric
source can leave an orientation null, and no acquisition design can remove an
unanchored global rigid gauge.

### 6.5 What would make the method genuinely SOM-specific?

Direct joint full-wave inversion can use the same raw measurements and parameters.
PRASC-SOM is scientifically distinct only through testable consequences of its
current-subspace structure:

1. the nested-rank theorem predicts and controls when current freedom absorbs pose;
2. the dual spectrum identifies which map modes can be safely unlocked after pose
   calibration;
3. early low-rank steps reduce current nuisance and computational dimension;
4. rank and frequency events define a continuation path tied to measured
   observability rather than a fixed schedule.

These mechanisms are hypotheses of algorithmic benefit, not proof of benefit. At
matched forward physics, initialization, priors, and stopping rules, the method
fails its main claim if it does not improve basin size, cost, or accuracy relative
to a well-tuned direct joint inversion.

## 7. Preliminary falsification evidence inherited from the prior loop

The existing deterministic 2D scalar Helmholtz loop is retained as a theory audit.
It was not designed for the corrected PRASC-SOM method and must not be presented as
its validation.

1. **Unrestricted-lift vacuity.** A full-row-rank sensing matrix reproduces a
   receiver-pose perturbation to relative residual \(5.5\times10^{-16}\). A
   rank-four retained basis leaves representative leakage \(0.784\). This supports
   the coordinate decomposition, not pose recovery.
2. **Moving SOM coordinates.** Near a spectral degeneracy, a hard retained
   projector jumps by \(1.424\), while a soft filter changes by \(0.127\). This
   supports event-aware coordinates, not improved imaging.
3. **Typed pose derivative.** Receiver re-sampling and transmitter-induced
   physical-current derivatives agree with centered finite differences, and their
   sum matches the total derivative to numerical precision.
4. **Rank saturation.** After correcting a relative-threshold bug, the nuisance
   range fills the real data space at \(r\ge M-1\) in the tested model, so visible
   pose rank is zero. The restricted solver is either a reparameterization of the
   direct objective or freezes hidden pose coordinates and performs worse.
5. **State-witness failure.** The tested \(T_U\), hard state equality, and soft
   variable projection do not discriminate correct pose or recover hidden
   directions. State constraints are therefore secondary hypotheses, not the
   present method's information source.
6. **Acquisition diversity.** Separated or directional transmitters remove a
   single isotropic point-source orientation deficiency in the tested setting, but
   do not remove the global map--pose gauge.

These results explain why the corrected method relies primarily on coherent phase,
rank control, frequency continuation, and independent acquisition diversity rather
than on unrestricted equivalent currents or a scalar state penalty.

## 8. Five theory-validation experiments for the next phase

No result from this section is claimed yet. The experiments are deliberately
falsification-first.

### E1. Coherent phase information versus phaseless information

Compute matched coherent and intensity likelihoods from the same generative model.
Across frequency and pose error, compare full and efficient Fisher spectra,
\(B_{\mathrm{vis}}\), Cramér--Rao bounds, local objective curvature, and empirical
basins. Verify the score-conditioning inequality and identify equality or
near-equality regimes. Null hypothesis: coherent phase provides no usable pose
advantage after nuisance elimination.

### E2. Nested-rank monotonicity and equality cases

At fixed linearizations, enumerate nested ranks and verify the Loewner differences
\(J_x(r)-J_x(r+1)\succeq0\), rank monotonicity, and equality characterization.
Construct an added current mode orthogonal to the pose tangent as an equality
control and a mode aligned with it as a strict-loss control. Then test whether the
local theorem predicts nonlinear calibration breakdown. Null hypothesis: the
spectral loss does not predict nonlinear behavior.

### E3. Dual \(\rho/B_{\mathrm{vis}}\) prediction

Generate map modes spanning low to high principal-angle overlap. Test whether
\(\rho_i\) predicts map degradation when pose is unknown and whether
\(\sigma_i(B_{\mathrm{vis}})\) predicts pose degradation when map/current vary.
Include priors, limited aperture, strong scattering, and gauge controls. Null
hypothesis: either spectrum is redundant with residual magnitude or singular values
of \(S\) alone.

### E4. Dual continuation versus strong baselines

Compare:

- wrong-geometry coherent SOM;
- established phaseless SOM;
- coherent SOM with fixed rank;
- direct joint full-wave inversion with the same data, forward physics,
  initialization, optimizer, priors, and compute budget;
- receiver/source-extension baseline;
- PRASC-SOM;
- known-pose and known-map oracles.

Sweep initial pose error, SNR, contrast, aperture, frequency ladder, and rank
schedule. Report basin volume, pose error, map error, complex residual, phase
residual, runtime, forward/adjoint counts, and failure rate. Before running, fix a
primary endpoint, a common forward/adjoint-solve budget, pose/map noninferiority
margins derived from wavelength, resolution, and noise, a seed/grid budget, and a
multiplicity-aware confidence procedure. A method-level advantage is accepted only
if the predeclared lower confidence bound on improvement is strictly positive on
the primary endpoint while all noninferiority gates pass. Include two adversarial controls:
the known limited-aperture false basin must test whether gate (36) prevents or fails
to prevent cycle skipping, and non-pose phase corruption (for example clock drift,
mutual coupling, or unmodeled phase noise) must establish a regime in which
PRASC-SOM is expected to lose. The method claim fails if adaptive SOM has no
predeclared advantage at matched budget or matched accuracy.

### E5. Stacked acquisition, symmetry, and gauge

Test single versus separated/directional transmitters, multiple frames,
asymmetric trajectories, and physically consistent multifrequency data. Include
two negative controls: duplicated symmetric acquisitions and an unanchored global
rigid transform. Verify (31), compare greedy design to random design, and report the
smallest quotient singular value before reconstruction. Null hypothesis:
acquisition scores do not predict recovery.

High-fidelity 3D Maxwell and measured-array experiments follow only if E1--E5
support the mechanism.

## 9. Remaining theorem package before a TGRS or TAP submission

The present draft is a defensible theory foundation, but the following items remain
open and are delegated to the accompanying GPT Pro problem set:

1. an efficient-information version of coherent-to-phaseless contraction after
   nuisance elimination, with equality and counterexamples under mismatched noise;
2. weighted/gauge-complete versions of Theorem 1 and Proposition 1, including
   singular \(K_0\), priors, and strictness conditions;
3. a fixed-rank local convergence theorem for the full PRASC-SOM update and an
   event-safe result or counterexample across rank/observability events;
4. a quantitative phase-basin condition that links pose covariance, aperture,
   path geometry, multiple scattering, and the next admissible frequency;
5. sufficient conditions and a tractable design guarantee for stacked
   acquisitions with shared map and frame-specific currents;
6. equation-level integration of the exact Twofold SOM domain fold;
7. a theorem or controlled counterexample showing when the SOM policy can improve
   complexity or basin relative to direct joint inversion;
8. a physical fixed-permittivity multifrequency and 3D Maxwell extension.

For a TGRS submission, the strongest route is an imaging-method paper with
comprehensive synthetic and measured reconstruction evidence. For a TAP
submission, the antenna/array calibration model, phase-center/mutual-coupling
nuisance, and electromagnetic validation must be correspondingly stronger.

## 10. Claim ledger

| Claim | Status in this draft |
|---|---|
| coherent SOM is phase sensitive | prior art / physical fact |
| phaseless SOM retains no current constraint | rejected |
| joint pose plus image estimation is new | rejected |
| unrestricted Green perturbation can be fit by current | often true but scientifically vacuous |
| retained lift plus leakage is well-defined after metric choice | proved, finite dimensional |
| increasing rank cannot improve pose information at a fixed linearization when effective nuisance spaces are nested | proved |
| \(\rho=\sin^2\theta\) on the current-cleaned observable support | proved under stated assumptions |
| phaseless processing cannot increase matched joint Fisher information | proved under regularity |
| medium/high-frequency phase can improve local pose precision | asymptotically supported; constants/basin open |
| state witness \(T_U\) recovers hidden pose | falsified in current experiments |
| naive visible-coordinate solver improves direct inversion | falsified in current experiments |
| rank/frequency-adaptive coherent SOM improves calibration | open; central experiment claim |
| stacked acquisition can restore quotient pose rank | exact criterion given; design guarantee open |

## 11. Limitations

1. The main results are local and finite dimensional.
2. The continuum compact-operator lift, nonclosed ranges, and regularized graph
   limits remain open.
3. The pose-information monotonicity holds at fixed \(A,B,S\), whitening, and gauge;
   it does not compare different frequencies or nonlinear iterates directly.
4. The information-contraction theorem requires a common coherent parent
   experiment or a parameter-independent channel. Different sensor noise laws
   require a separate comparison.
5. The large-\(kr\) phase scaling does not prove a global basin under limited
   aperture or multiple scattering.
6. Current preliminary experiments use a 2D scalar Helmholtz model and
   low-dimensional real contrast.
7. The exact Twofold SOM domain fold is not yet implemented or claimed.
8. No current experiment demonstrates PRASC-SOM reconstruction superiority.
9. Prior-art retrieval was bounded and partly rate-limited; global novelty is not
   asserted.

## 12. Conclusion

The corrected research direction is not “replace unknown Green operators by
equivalent currents.” Unrestricted current equivalence is too flexible to add
identifying information. Nor is the contribution simply “estimate pose and image
together.” The central proposal is to preserve coherent phase and use SOM's
retained-current spectrum as a control variable for self-calibration.

The theory reveals a precise tension. Raising retained rank can recover richer
current content, but at a fixed linearization it monotonically enlarges the
nuisance space that can impersonate pose. The dual spectra
\(B_{\mathrm{vis}}\) and \(\rho\) expose the corresponding loss on the pose and map
sides. Meanwhile, frequency increases local phase precision while shrinking the
unambiguous displacement scale. These relations motivate a dual continuation:
calibrate with low frequency and conservative rank, then unlock medium/high
frequency and higher rank only when the observed information margins permit it.

That mechanism is the paper's plausible scientific value. Its algorithmic advantage
over direct joint inversion remains an explicit falsifiable claim for the next
experimental phase.

## Reproducibility and AI-use statement

The prior autonomous research loop, raw deterministic outputs, figures, negative
results, parent numerical-rank correction, and independent audits are preserved
under the project research and experiment directories. Those artifacts are used as
constraints on this corrected theory; they are not relabeled as validation of the
new method. The autonomous workflow generated and falsified hypotheses, while the
top-level scientific review retained responsibility for type correctness, theorem
claims, novelty boundaries, and interpretation.

## References

1. X. Chen, “Subspace-Based Optimization Method for Solving Inverse-Scattering
   Problems,” *IEEE Transactions on Geoscience and Remote Sensing*, 2010.
   [DOI](https://doi.org/10.1109/TGRS.2009.2025122)
2. Y. Zhong and X. Chen, “Twofold Subspace-Based Optimization Method for Solving
   Inverse Scattering Problems,” *Inverse Problems*, 25, 085003, 2009.
   [DOI](https://doi.org/10.1088/0266-5611/25/8/085003)
3. L. Pan, Y. Zhong, X. Chen, and T. S. Yeo, “Subspace-Based Optimization Method
   for Inverse Scattering Problems Utilizing Phaseless Data,” *IEEE Transactions
   on Geoscience and Remote Sensing*, 2011.
   [DOI](https://doi.org/10.1109/TGRS.2010.2070512)
4. L. Bellomo, S. Pioch, M. Saillard, and K. Belkebir, “An Improved Antenna
   Calibration Methodology for Microwave Diffraction Tomography in Limited-Aspect
   Configurations,” *IEEE Transactions on Antennas and Propagation*, 2014.
   [DOI](https://doi.org/10.1109/TAP.2014.2308534)
5. G. Huang, R. Nammour, and W. Symes, “Full-Waveform Inversion via
   Source-Receiver Extension,” *Geophysics*, 2017.
   [DOI](https://doi.org/10.1190/geo2016-0301.1)
6. L. Métivier and R. Brossier, “Receiver-Extension Strategy for Time-Domain
   Full-Waveform Inversion Using a Relocalization Approach,” *Geophysics*, 2021.
   [DOI](https://doi.org/10.1190/geo2020-0922.1)
7. Y. Li, K. Lee, and Y. Bresler, “Identifiability in Bilinear Inverse Problems
   With Applications to Subspace or Sparsity-Constrained Blind Gain and Phase
   Calibration,” *IEEE Transactions on Information Theory*, 2017.
   [DOI](https://doi.org/10.1109/TIT.2016.2637933)
8. G. R. Karthik and P. K. Ghosh, “A Scalable Deep Learning Model for
   Simultaneous Reconstruction and Transmitter Localization in Inverse
   Scattering,” *PIERS*, 2023.
   [DOI](https://doi.org/10.1109/PIERS59004.2023.10221374)
