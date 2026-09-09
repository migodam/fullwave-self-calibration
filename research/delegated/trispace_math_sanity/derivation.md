# Equivalent-current representation of $dG(X)[\delta X]j$ — derivation

All statements are labeled. Notation follows the report and
`SOM_SLAM_THEORY_CONTEXT.md`: current space $\mathcal H_J$, data space
$\mathcal Y$ (whitened, realified where needed), sensing operator
$G_S:\mathcal H_J\to\mathcal Y$, domain propagator $G_D$, contrast $\chi$,
$M(X)=I-D_\chi G_D(X)$, induced current

$$
j(X)=M(X)^{-1}D_\chi\,e^{\mathrm{inc}}(X),
\qquad
F(\chi,X)=G_S(X)\,j(X).
$$

## 1. Pose tangent decomposition

**[proposition candidate]** With a world-fixed grid ($D_XG_D=0$,
$D_XD_\chi=0$) and pose coordinates $\{X_l\}$,

$$
B_l:=D_{X_l}F
=\underbrace{dG_S(X)[\delta X_l]\,j}_{\text{sensing path}}
+\underbrace{G_S(X)\,M(X)^{-1}D_\chi\,d_Xe^{\mathrm{inc}}(X)[\delta X_l]}
_{\text{interior path}}.
$$

Proof sketch: differentiate $G_S\,M^{-1}D_\chi e^{\mathrm{inc}}$; the
resolvent term vanishes because $D_XG_D=0$. If $G_D$ and $D_\chi$ do depend on
$X$ (moving grid, pose-dependent parameterization), add
$G_S M^{-1}\big[D_\chi\,d_XG_D[\delta X_l]j + d_XD_\chi[\delta X_l]E^{\mathrm{tot}}\big]$
inside the bracket (cf. report Remark 4.54). **[derivation]**

Verified in the harness (`src/helmholtz.py`, `build_AB` plus explicit Green
gradients) at $N=8$, one pose: $\|B_l-(\text{sensing}+\text{interior})\|\le
2.6\times10^{-18}$, and the finite-difference check of $B$ against
$F(\chi,X\pm\delta X)$ gives relative error $1.0\times10^{-9}$. The rotation
column ($l=\theta$) has zero interior part because the transmitter is at the
platform rotation center ($dE^{\mathrm{inc}}/d\theta=0$).

## 2. The sensing-path equivalence question

### 2.1 Exact representability

**[proposition candidate]** Let $v=dG_S(X)[\delta X]\,j\in\mathcal Y$ and
assume $\operatorname{Ran} G_S$ is closed. Then there exists
$\delta j\in\mathcal H_J$ with $G_S(X)\,\delta j=v$ if and only if

$$
v\in\operatorname{Ran} G_S(X).
$$

The minimum-norm solution is unique in $(\ker G_S)^\perp$ and equals
$\delta j_{\mathrm{eq}}=G_S^\dagger v$; the full solution set is
$\delta j_{\mathrm{eq}}+\ker G_S$. Proof: standard Moore–Penrose/closed-range
orthogonal decomposition $\mathcal Y=\operatorname{Ran} G_S\oplus
\ker G_S^*$. **[derivation]**

### 2.2 Residual

**[proposition candidate]** For any $v\in\mathcal Y$ the least-squares fit of
$G_S\delta j$ to $v$ has residual

$$
r(v)=v-G_S G_S^\dagger v
=P_{\ker G_S^*}v
=P_{\operatorname{Ran} G_S^\perp}v,
\qquad
\|r(v)\|=\operatorname{dist}(v,\operatorname{Ran} G_S).
$$

The angle $\theta(v,\operatorname{Ran} G_S)$ satisfies
$\sin^2\theta=\|r\|^2/\|v\|^2$. With data whitening $W$ all formulas use
$\operatorname{Ran}(WG_S)$; realification is an isometry and commutes with
real projectors, so no correction is needed beyond stacking real/imaginary
parts. **[derivation]**

### 2.3 Non-closed range

**[proposition candidate]** If $\operatorname{Ran} G_S$ is not closed and
$v\in\overline{\operatorname{Ran} G_S}\setminus\operatorname{Ran} G_S$, no
exact $\delta j$ exists; zero residual is approached only asymptotically by
$\delta j_n$ with $\|\delta j_n\|\to\infty$. A Tikhonov sequence
$\delta j_\varepsilon=G_S^*(G_SG_S^*+\varepsilon I)^{-1}v$ satisfies
$G_S\delta j_\varepsilon\to v$ while $\|\delta j_\varepsilon\|\to\infty$.

**[derivation]** Minimal model: $Ge_n=e_n/n$ on $\ell^2$, $v_n=1/n$. Then
$v\in\overline{\operatorname{Ran} G}\setminus\operatorname{Ran} G$,
$(\delta j_\varepsilon)_n=(1+n^2\varepsilon)^{-1}$, and
$\|\delta j_\varepsilon\|^2\sim\frac{\pi}{4\sqrt\varepsilon}\to\infty$ while
the residual tends to $0$. This is the equivalent-current analog of the
report's Theorem 4.3 / P1 asymptotic-compensation distinction; the
Moore–Penrose representative does not exist.

## 3. The full pose tangent

**[proposition candidate]** With the decomposition of §1,

$$
B_l\in\operatorname{Ran} G_S
\iff
dG_S[\delta X_l]\,j\in\operatorname{Ran} G_S,
$$

because the interior term is $G_S$ applied to a current. When the condition
holds, the canonical equivalent current is

$$
\delta j_{\mathrm{equiv},l}
=\underbrace{M^{-1}D_\chi\,d_Xe^{\mathrm{inc}}[\delta X_l]}_{\delta j_{X,l}}
+G_S^\dagger\,dG_S[\delta X_l]\,j,
$$

and $G_S\,\delta j_{\mathrm{equiv},l}=B_l$; the residual of $B_l$ equals the
residual of the sensing term alone. **[derivation, verified]** At $N=8$,
rank $4/4$: $\|G_S\delta j_{\mathrm{equiv},l}-B_l\|\le 9\times10^{-16}$.

### 3.1 Finite-m degeneracy

**[proposition candidate]** If $\dim\mathcal Y=m\le S=\dim\mathcal H_J$ and
$G_S$ has full row rank, then $\operatorname{Ran} G_S=\mathcal Y$ and **every**
$dG_S[\delta X]j$, hence every $B_l$, is exactly representable: the residual
is identically zero, the pose-induced data subspace
$W_P:=\operatorname{span}\{B_l\}\subseteq\mathcal Y=\operatorname{Ran} G_S$,
and the pullback $G_S^\dagger(\mathcal Y)=(\ker G_S)^\perp$ fills every
visible current direction — no data-level current-vs-pose distinction exists
at a single snapshot. Consequences:

- TriSpace is non-degenerate only (i) across stacked poses/frequencies where
  the *same* $\delta X$ must explain all blocks, or (ii) inside a retained
  current subspace $\operatorname{Ran} U_T$, or (iii) in the rank-deficient /
  continuum regimes of §2.2–2.3.
- The defined $V_P:=\operatorname{span}\{G_S^\dagger B_l\}$ has
  $\dim V_P\le p$ and carries no information not already available to
  currents; the correct single-snapshot discriminants are the minimum-norm
  representative and the retained-space residual
  $r_{\mathrm{SOM},l}=P_{(\operatorname{Ran} G_S P_T)^\perp}B_l$, not
  $P_{\operatorname{Ran} G_S^\perp}B_l$.

**[derivation]** Linear algebra + the verified rank-$4/4$ run
(residuals $<5\times10^{-18}$ for all three pose coordinates).

### 3.2 Bridge to the $A/B$ layer

**[proposition candidate]** Write $A=G_ST_\chi$ with
$T_\chi=M^{-1}\operatorname{diag}(E^{\mathrm{tot}})$. If
$B\delta X\in\operatorname{Ran} G_S$ (equivalently the sensing residual
vanishes), then

$$
B\delta X\in\operatorname{Ran} A
\iff
G_S^\dagger B\delta X\in\operatorname{Ran} T_\chi+\ker G_S.
$$

If $B\delta X\notin\operatorname{Ran} G_S$, then automatically
$B\delta X\notin\operatorname{Ran} A$. Proof sketch: under $r=0$,
$B\delta X=G_S(G_S^\dagger B\delta X)$; apply the report's Theorem 4.66 range
condition $G_Su\in\operatorname{Ran}(G_ST_\chi)\iff
u\in\operatorname{Ran}T_\chi+\ker G_S$ to $u=G_S^\dagger B\delta X$.
This translates the map–pose confusability question (the $K_{\mathrm{eff}}$
layer) into a current-space membership question **without** identifying
$K_{\mathrm{eff}}$ eigenvectors with SOM modes. **[derivation]**

## 4. Interior $G_D$ vs sensing $G_S$ perturbations

**[derivation]** The structural differences:

1. *Ambient space and residual.* $dG_S[\delta X]j$ lives in data space and
   needs the §2 pullback; the interior path produces $\delta j_X$ directly in
   current space, so its data image is automatically in
   $\operatorname{Ran} G_S$ — the interior path can never violate the range
   condition.
2. *Resolvent coupling.* The interior path is amplified by
   $M^{-1}$ ($\|M^{-1}\|\le m_0^{-1}$, report Theorem 4.52); the sensing path
   has no resolvent factor. Near resonance the interior/sensing sensitivity
   ratio diverges; the two cannot share one robustness constant.
3. *Spectral object.* $G_S$ maps to $\mathbb C^m$ (rank $\le m$, fast decay);
   $G_D$ is a Hilbert–Schmidt domain propagator with log-singular kernel; its
   coordinate derivatives are not uniformly bounded without standoff/self-cell
   handling (report Remark 4.56). $V_S^\pm$ and $V_D^\pm$ therefore have
   different rank/decay — one reason the TSOM intersection is genuinely
   non-commuting rather than a cosmetic basis choice.
4. *Visibility vs realizability.* A pose perturbation of $G_S$ changes which
   currents are visible — it can reveal $\ker G_S$ content (the residual of
   §2.2 is exactly this effect). A pose perturbation of $e^{\mathrm{inc}}$ /
   $G_D$ changes which currents exist. The former is an autofocus/phase-center
   degeneracy; the latter is an excitation change.
5. *Structure preservation.* **[proposition candidate]** For body-fixed rigid
   motions, $dG_S[\delta X]$ preserves any null-row symmetry of $G_S$
   (coincident receivers move identically), hence
   $dG_S[\delta X]j\in\operatorname{Ran} G_S$ even when $G_S$ is
   rank-deficient. Per-channel phase/gain calibration does not respect the
   symmetry and is exactly the perturbation that exits
   $\operatorname{Ran} G_S$. Verified numerically (§2.2 of
   counterexamples.md).
6. *Regime dependence.* **[conjecture]** In the Born far-field model both
   paths reduce to multiplicative phase factors
   ($-ik\,d\cdot\delta p_T$ for Tx motion, $+ik\,s\cdot\delta p_R$ for Rx
   motion, report convention 4.15), so interior and sensing pose coordinates
   are mutually confusable under narrowband phase wrap (cycle skipping) and
   monostatic geometry $d\approx-s$. With $G_D\neq0$ the interior path also
   perturbs the multiple-scattering resolvent, producing contrast- and
   frequency-dependent signatures that pure re-phasing cannot reproduce;
   the separation shrinks to the Born degeneracy as $\chi\to0$. Quantitative
   verification is **[open]**.

## 5. A well-posed $V_P$ and the three residuals

**[proposition candidate]** At a nominal $(\chi,X)$ define, for each pose
coordinate $l$,

$$
\delta j_{\mathrm{equiv},l}=\delta j_{X,l}+G_S^\dagger\,dG_S[\delta X_l]\,j,
\qquad
V_P:=\operatorname{span}\{\delta j_{\mathrm{equiv},l}:l=1..p\},
$$

with $\dim V_P\le p$ ($p=$ pose dimension, $3T$ stacked). $V_P$ is a tangent
object, canonical via the $(\ker G_S)^\perp$ representative, and
**scene-dependent** ($j$, $\delta j_X$ depend on $\chi$) — unlike $V_S^\pm$,
which depends only on $X$. A "precomputed" $V_P$ classification is therefore
not scene-independent in full-wave. Three nested residuals must be kept
separate:

$$
r_S=P_{\operatorname{Ran} G_S^\perp}\,dG_S[\delta X_l]j
\ \text{(pure geometry, current-unrepresentable data)},
$$

$$
r_{\mathrm{SOM},l}=P_{(\operatorname{Ran} G_S P_T)^\perp}\,B_l
\ \text{(representable in full current space, not in the retained basis)},
$$

$$
r_{\mathrm{Schur}}=\min_{h}\|Au-Bh\|
\ \text{(map--pose confusability at the }A/B\text{ level)}.
$$

They answer different questions and must not be summed or interchanged.
**[derivation]**

## 6. A coherent TriSpace SOM self-calibration algorithm

**[conjecture / proposal]** Objective (one function, gauge-anchored):

$$
\min_{\chi,X}\ \|y-F(\chi,X)\|^2_W+\lambda R(\chi)+h^\top J_X h
\quad\text{s.t.}\quad
j=M^{-1}D_\chi e^{\mathrm{inc}}(X),\ \
j\in\operatorname{Ran} U_T(X),
$$

with the global $SE(d)$ gauge fixed by an anchor (fixed pose 0 / known
background / external phase reference), otherwise the objective is constant
along gauge orbits and pose/map are interchangeable.

1. **Init.** $X^{(0)}$ from odometry/prior; known-geometry SOM reconstruction
   of $\chi^{(0)}$.
2. **S/D fold at fixed $X$.** SVD of $G_S(X)$; gap-selected
   $V_S^+(X)$; TSOM retained basis $U_T=[V_S^+,\,B_{S^-D^+}]$ (via CS
   decomposition, not $P_SP_D$ products); constrained contrast-source solve.
3. **Pose block.** Gauss–Newton step
   $\delta X=\arg\min_\delta\|y-F(\chi,X)-B(X)\delta\|^2_W+\delta^\top
   J_X\delta$. Use the residual split: data components in $r_S$ identify pose
   without current competition; components in
   $\operatorname{Ran}(G_S P_T)$ are jointly explained by
   $(\delta c,\delta X)$ via the Schur solve, so pose and retained currents
   must be updated together, not alternately on the same residual.
4. **Map/current block.** Reconstruct $\chi$ from the constrained current
   ($\chi_n=j_n/E^{\mathrm{tot}}_n$ on cells with
   $|E^{\mathrm{tot}}_n|$ above threshold, or a joint $(\chi,c)$ solve with
   $c=P_T j(\chi)$).
5. **Subspace refresh guard.** Recompute $U_T$ only when the pose change
   exceeds a trust radius or a singular value approaches the truncation
   threshold (rank event); otherwise keep the projector fixed so the block
   descent stays on one stratum.
6. **Certificates per iterate.** (i) gauge anchored; (ii)
   $\sigma_{\min}$ of the reduced joint Hessian on the quotient; (iii)
   per-coordinate angles between pose-induced data directions and
   $\operatorname{Ran}(G_S P_T)$; (iv) rank-event flags; (v) the three
   residuals of §5 against the noise floor.

**[conjecture]** With Lipschitz derivatives and fixed-rank strata, the
iteration is a block-coordinate descent on a single non-convex objective, so
subsequences converge to stationary points; rate and basin guarantees are
**[open]**.

**[open]** (a) When $G_S$ has full row rank at finite $m$, $r_S\equiv0$ and
the "third space" collapses to the retained-space and $A/B$ layers — whether
TriSpace then offers more than SOM-regularized joint inversion is unproven;
(b) $V_P$ scene-dependence vs. amortization; (c) rank-event crossing
formulation; (d) exact $V_D^+$ definition per the cited Chen TSOM version;
(e) non-convex local minima and cycle skipping are not addressed.
