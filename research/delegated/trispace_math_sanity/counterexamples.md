# Counterexamples, identifiability and gauge failures

Each entry: what it kills, construction, and verification status. No novelty
claim is implied by any of these; several are finite-dimensional special cases
of statements already in the main report (P1–P10), restated here at the
TriSpace/current level.

## C1. Single-snapshot full-row-rank degeneracy

Kills: any claim that $V_P$ is a proper current subspace at a single snapshot
with $m\le S$ receivers/pixels.

Construction: $G_S:\mathbb C^S\to\mathbb C^m$ generic with $m\le S$ and full
row rank ⇒ $\operatorname{Ran} G_S=\mathbb C^m$; every
$dG_S[\delta X]j$ is representable and the residual is identically $0$, so
the pose-induced data subspace $W_P$ is entirely inside the current-visible
data space — no data-level pose/current distinction exists, and the pullback
$G_S^\dagger(\mathcal Y)=(\ker G_S)^\perp$ fills every visible direction.
Minimal instance: $G_S=[1,0]:\mathbb C^2\to\mathbb C$,
$dG_S[\delta X]j=(\delta a)j_1\in\mathbb C=\operatorname{Ran} G_S$.

Status: **[derivation]** + verified in the harness at $N=8$, one pose,
rank $4/4$: sensing residuals $<5\times10^{-18}$, equivalent current
reproduces $B_l$ to $9\times10^{-16}$. The non-degenerate formulations are
stacked-pose coherence, retained-subspace, rank-deficient, and continuum
versions (C2–C4, C8).

## C2. Structure-preserving rigid motion (negative counterexample)

Kills: the intuition "rank-deficient $G_S$ ⇒ pose perturbations leave
$\operatorname{Ran} G_S$".

Construction: two receivers exactly coincident ⇒
$G_S$ has a null-row symmetry $e_1-e_2\in\ker G_S^*$. Any body-fixed rigid
motion moves both receivers identically, so the derivative rows remain
coincident and $dG_S[\delta X]j$ stays orthogonal to $e_1-e_2$; the residual
vanishes for **all** $\delta X$. Verified: coincident pair, rank $3/4$, residual
ratio $1.3\times10^{-15}$ under $\delta p_x$.

Status: **[proposition candidate]** (proof sketch: rigid motion is an
isometry commuting with the coincidence projector) + numeric verification.
Consequence: the range condition is only *exercised* by perturbations that do
not respect the null-space symmetries of $G_S$ — i.e., per-channel
calibration, independent antenna motion, clock skew (C3).

## C3. Per-channel calibration exits the range (the minimal unrepresentable direction)

Kills: "every pose/calibration perturbation has an equivalent current".

Construction: rank-deficient $G_S$ as in C2. A per-channel phase calibration
$dG_S=i\delta\varphi\,e_1(G_Sj)^\top$ restricted to row 1 gives
$v=i\delta\varphi\,(G_Sj)_1 e_1$. Since $e_1\notin\operatorname{Ran} G_S$
(receiver 2 would have to read zero at the same point), the residual is
$P_{\operatorname{Ran} G_S^\perp}v\neq0$. Verified numerically: residual ratio
$1/\sqrt2$ (single channel; matches projection onto
$(e_1-e_2)/\sqrt2$) and $\approx0.89$ for an all-channel phase vector; the
least-squares current fit leaves exactly that residual.

Status: **[derivation]** + verified. Minimal abstract form:
$G_S=\begin{pmatrix}1&0\\0&0\end{pmatrix}$ (an invisible mode
$j_0=(0,1)^\top\in\ker G_S$),
$dG_S=\begin{pmatrix}0&0\\0&1\end{pmatrix}$ (a perturbation that couples the
invisible mode to receiver 2) ⇒
$v=dG_S j_0=(0,1)^\top\notin\operatorname{Ran} G_S=\operatorname{span}\{e_1\}$:
residual ratio $1$. This is the exact current-level picture of "an invisible
mode turned visible by pose, yet not explainable by any current at the nominal
$G_S$".

## C4. Continuum non-closed range: asymptotic-only compensation

Kills: unconditional existence of a minimum-norm equivalent current in the
continuum; silent use of $G_S^\dagger$ without closed range.

Construction: $Ge_n=e_n/n$ on $\ell^2$ (compact, dense range), $v_n=1/n$.
$v\in\overline{\operatorname{Ran} G}\setminus\operatorname{Ran} G$; no exact
$\delta j$ exists, truncations $\delta j^{(N)}=1_{[1..N]}$ give
$G\delta j^{(N)}\to v$ with $\|\delta j^{(N)}\|\to\infty$; Tikhonov
$\delta j_\varepsilon$ has $\|\delta j_\varepsilon\|^2\sim\frac{\pi}{4\sqrt
\varepsilon}\to\infty$ with residual $\to0$.

Status: **[derivation]** (standard compact-operator pathology, same mechanism
as report Theorem 4.3 for $\operatorname{Ran} B$). Any TriSpace statement
about $V_P$ in the continuum must therefore carry a norm budget or a
regularization, or be restricted to closed-range / finite-$m$ sensing.

## C5. Global rigid-motion gauge: pose ⇄ map exact counterfeiting

Kills: any TriSpace classification that treats "pose-confounded current" and
"map-generated current" as disjoint from data alone.

Construction: report Theorem 4.9/Corollary 4.10: $F(g\chi,gX)=F(\chi,X)$ ⇒
$A\xi_\chi+B\xi_X=0$. Then $B\xi_X=-A\xi_\chi\in\operatorname{Ran} A\subseteq
\operatorname{Ran} G_S$: the pose direction has an **exact** equivalent
current with **zero** residual, namely
$\delta j_{\mathrm{eq}}=-T_\chi\xi_\chi=-M^{-1}\operatorname{diag}
(E^{\mathrm{tot}})\xi_\chi$, a genuine physical map change. The gauge pair
sits in $V_P\cap\operatorname{Ran} T_\chi$: perfect confusability that only an
anchor/prior resolves. In Born form this is
$-i(q\cdot a)\hat\chi$ vs $+i(q\cdot a)\hat\chi$ cancellation (convention
4.15).

Status: **[theorem candidate]** in the report; restated here at current level.
Numeric analogue already exercised by `family3_gauge_born.py`; a current-level
restatement is part of experiment E4.

## C6. Pose-only stabilizer (radial scene)

Kills: the assumption that every pose tangent direction has nonzero data
effect.

Construction: $\chi(x)=\varphi(\|x-c\|)$; rotation about $c$ gives
$\xi_\chi=0$ (map stabilizer) and $B\xi_X=0$ (report Proposition 4.11). The
pose coordinate is unobservable from data: no current, equivalent or
otherwise, can detect it — a pure pose nullspace, not a map–pose confound.
Unanchored motion graphs add the same global-mode kernel on the pose side.

Status: **[proposition candidate]** in the report; needs the discrete-pixel
caveat (piecewise-constant pixels do not contain $\xi_\chi$ exactly —
representation error, report Proposition 4.12(5)).

## C7. Born narrowband / monostatic phase aliasing

Kills: local identifiability claims extending to global; "interior vs sensing
are always separable".

Construction: under far-field Born (convention 4.15) Tx motion gives phase
$-ik\,d\cdot\delta p_T$ and Rx motion $+ik\,s\cdot\delta p_R$ on the **same**
datum. For monostatic $d\approx-s$ they enter identically; for any geometry
the phases are only known modulo $2\pi$, so large pose errors alias
(cycle-skipping; cf. report Counterexample 4.82). As $\chi\to0$
($G_D$-coupling vanishes) the full-wave interior/sensing distinction
degenerates to this bilinear phase ambiguity; identifiability of the split is
therefore contrast- and bandwidth-dependent.

Status: **[derivation]** for the Born phase factors; the contrast-dependent
full-wave separation is **[conjecture]**, quantitative test **[open]**.

## C8. Retained-subspace leakage and non-commuting folds

Kills: "pose directions live in the retained deterministic subspace" and any
"eight-cell" exact decomposition.

Construction: full-rank $G_S$ (C1) makes $r_S=0$, but with TSOM projector
$P_T$ the residual $r_{\mathrm{SOM},l}=P_{(\operatorname{Ran} G_S P_T)^\perp}
B_l$ is generically nonzero: the pose direction is representable only through
discarded current modes. Related non-commutativity: the dominant current mode
of $G_S^*G_S$ need not be the dominant mode of $G_S^*WG_S$ (report
Counterexample 4.73, $G_S=\operatorname{diag}(2,1)$,
$\operatorname{Ran} B=\operatorname{span}\{(1,1)\}$), so "SOM-first" vs
"pose-elimination-first" re-selection disagree; $P_SP_P$ is not the
intersection projector (report Theorem 4.74). Also
$\|P_S^-(X^\star)P_S^+(\hat X)\|$ measures misclassification leakage and can
jump at rank events (report Theorem 4.48: projectors on different ranks are
distance $1$).

Status: **[theorem candidates]** in the report for the algebra; the leakage
quantities are numeric targets of experiment E5. `family8_som_distinction.py`
already verified the distinctness of $V_S$ vs $V_A$ (right subspaces) but not
the $V_P$ construction.

## C9. Moving-grid / pose-dependent map parameterization

Kills: coordinate-independent meaning of $V_P$ when $\chi$ is stored on a
robot-local grid.

Construction: if $G_D$ or $D_\chi$ depend on $X$ (grid rides on the
platform), the interior bracket gains $D_\chi\,d_XG_D[\delta X_l]j$ and
$d_XD_\chi[\delta X_l]E^{\mathrm{tot}}$ terms (report Remark 4.54), and the
"equivalent current" mixes grid-motion with physical current change; $V_P$
depends on the chosen chart. The harness is world-fixed, so its results do not
automatically transfer.

Status: **[derivation]** (bookkeeping); chart-dependence of the resulting
classification is **[open]**.

## Summary of gauge/identifiability failure modes

1. $SE(d)$ gauge: pose and map interchangeable for free (C5).
2. Pose-only nullspace / stabilizer (C6) + unanchored motion-graph kernel.
3. Full-row-rank single snapshot: everything is an equivalent current (C1).
4. Rank-deficient + symmetry-preserving motion: still everything (C2);
   symmetry-breaking calibration: residual, genuinely pose-only data (C3).
5. Non-closed range: no minimum-norm representative (C4).
6. Born phase wrap / monostatic: global aliasing, contrast-dependent
   interior/sensing split (C7).
7. Truncation/leakage + non-commuting folds (C8).
8. Coordinate/parameterization dependence (C9).
