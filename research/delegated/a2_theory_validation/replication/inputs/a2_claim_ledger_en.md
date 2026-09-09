# English claim / type / assumption ledger extracted from `inputs/A2.md`

Focus requested: **rank, bias, and acquisition** claims with their evidence state
([proved] / [proved recombination] / [conditionally proved] / [falsified] / [open]) and
the type/assumption rules those claims depend on.

Source ranges used: §1 verdict 1–42; §2 type ledger 43–238; §3 theorem package 243–1750;
§4 verdict table 1823–2057; §6 claim ledger 2111–2136; §8 candidate-combination audit
2423–2437; §9 unresolved list 2439–2468.

## A. Type / assumption ledger (from §2) — rules the experiments must encode

### A.1 Three ranks are different objects (from §1, §2.1)

| Quantity | Meaning | Operational rule |
|---|---|---|
| $L_{\mathrm{det}}$ | number of data-determined current components (original SOM cutoff) | Increasing it shrinks the complementary ambiguous space; it is not the free-nuisance rank |
| $r_{\mathrm{num}}$ | numerical state-approximation dimension | Refining it does not automatically reduce the physical model’s own Fisher information |
| $r_{\mathrm{free}}$ | extra free-current nuisance dimension | This is the quantity subject to the exact information-loss/monotone-erosion accounting |

Original SOM’s deterministic cutoff must not be identified with a free coefficient rank;
that identification is explicitly marked **falsified** (claim ledger, §6).

### A.2 Parameterization rules (from §2.2)

1. Independent-current form (2.1) and reduced physical form (2.2) are both legal
   parameterizations of the same physical experiment when the state equation is exact.
2. A block of the form $Kc+A\delta\chi+Bh$ (2.3) may only be used after explicitly declaring
   a model-error/current-correction (conservative envelope) model; otherwise the
   parameterizations are being added to fake a third model.
3. Fixed-dielectric multifrequency data cannot be produced by reusing one unexplained
   unscaled scattering potential at different frequencies.

### A.3 Whitening/realification factor (from §2.3)

The real covariance of $\mathcal E(W\varepsilon)$ for proper circular complex noise is
$I/2$. Use convention (2.4) (factor $\sqrt2$) so the real noise covariance is $I$, or keep
the complex embedding and record Fisher $=2\times$ Gram. The factor does not change rank,
principal angles, or relative $\rho$, but changes absolute precision gates. Improper or
across-acquisition correlated noise requires full real-covariance whitening.

### A.4 Gauge/quotient rules (from §2.4)

- Joint-quotient reduction precedes map/pose blocking; efficient information is invariant to
  the right-inverse choice.
- Absolute pose is not a well-defined target on an unanchored quotient; an external anchor or
  explicit gauge fixing is required.
- Unanchored additional acquisitions do not remove the global rigid gauge; source-model
  symmetry stabilizers are a different object.

### A.5 Metric/pseudoinverse rules (from §2.5)

Weighted pseudoinverse (2.5), current projector (2.6), coefficient lift (2.7) as transcribed
in the experiment-blueprint extraction. Only $H_S$ is lifted; a full-$B$-then-$H_D$ lift is
not the same operation. Data/state quantities must not be conflated with priors or parameter
metrics.

## B. G0 minimal models and physical counterexample (from §3.1) — statuses

| Item | Statement | Status in A2.md |
|---|---|---|
| G0-A | Complete, exact current coordinates create no information advantage over the reduced physical form (same predictions, feasible points, objectives, Fisher information) | [proved]; reparameterization/preconditioning only |
| G0-B | Truncation can strictly reduce pose risk iff $c_*^2<\sigma^2/\epsilon^2$ for $y_1=h+c+\varepsilon_1$, $y_2=\epsilon c+\varepsilon_2$; with a prior $C^2<\sigma^2/\epsilon^2$ it is truth-free | [proved] as a bias–variance/regularization statement; **not** information creation and not SOM-exclusive |
| G0-C | Small residual + good truncated Fisher does not imply correct calibration (zero-noise counterexample with truncation bias $T$ and arbitrarily small residual $\epsilon T$) | counterexample; implication [falsified] |
| G0-D | Zero envelope visibility ($B_{\mathrm{vis}}=0$ after admitting a free complex current) does not prove physical-model nonidentifiability | [proved]; narrows the old nuisance-saturation explanation; old matrix computations are not overturned |

## C. Rank / bias / acquisition claims with evidence states (focus ledger)

From §6 claim ledger (2111–2136), §3 status labels (533, 1145), and §4 verdict table
(2039–2057). Evidence-state terminology: `[proved]`, `[proved recombination]`,
`[conditionally proved]`, `[falsified]`, `[open]`; a few are `[numerically supported by old
loop]`.

| Claim | Evidence state | Originality/scope state |
|---|---|---|
| Matched coherent-to-phaseless joint Fisher contraction | [proved] | [textbook] |
| Nuisance-eliminated efficient-information contraction, including equality condition | [proved] | [proved recombination] |
| General coherent experiments are strictly better than phaseless | negated by counterexample | [falsified] |
| $\rho$ and the normalized pose Schur spectrum share canonical correlations | [proved] | [textbook / proved recombination] |
| Original $\rho$ equals $\sigma(B_{\mathrm{vis}})$ | negated | [falsified] (raw $B_{\mathrm{vis}}$ carries absolute scale; $\rho$ is known-pose-normalized; $B=0$ example) |
| Nested free-current space makes $J_x$ Loewner-nonincreasing | [proved] | physical application of a generic projection result |
| Original SOM’s $L_{\mathrm{det}}$ is the free-nuisance rank | negated | [falsified] |
| $\rho(r)$ is monotone | negated (Eq. 3.19: $3/4\to1\to0$) | [falsified] |
| Calibration-neutral current admission ($V^TQZ=0$) | [proved] | [proved recombination]; SOM control use is a candidate, not new linear algebra (state label line 1145) |
| Rank–acquisition information budget $J_{\mathrm{final}}-J_{\mathrm{original}}=\mathcal I_{\mathrm{acq}}-\mathcal L_{\mathrm{rank}}$; $J_{\mathrm{final}}\succeq J_{\mathrm{original}}\iff\mathcal I_{\mathrm{acq}}\succeq\mathcal L_{\mathrm{rank}}$ | [conditionally proved] (under declared model, nested rank, invertible shared-map innovation assumptions) | retrieval-bounded novel candidate; not global novelty |
| Truncation can lower pose risk | [conditionally proved] (G0-B condition) | concretization of bias–variance principle |
| Envelope zero visibility equals physical nonidentifiability | negated (G0-D) | [falsified] |
| $T_U$ / old state penalty restored hidden pose | [falsified by old loop] | limited to mechanisms and ranges actually tested |
| Unrestricted lift residual approaches machine precision | [numerically supported by old loop] | vacuity control, not a positive contribution |
| Multisource illumination removes the global map–pose gauge | negated | [falsified interpretation] |
| PRASC-SOM obtains a larger nonlinear basin | [open] | cannot be written into results |
| Safeguarded stationarity on a fixed stratum | [conditionally proved] | generic optimization theory; not core novelty |
| Lower total cost than tuned direct joint inversion | [open] | must count all auxiliary costs |

### §4.6 formal verdict table (equivalent focus items)

| Claim | Verdict |
|---|---|
| Complete SOM reparameterization automatically improves information or basin | **Negated** |
| Increasing numerical state rank automatically erodes physical pose Fisher | **Negated** |
| Information loss of free-current expansion is exactly computable | **Proved** |
| Calibration-neutral current additions can be selected | **Proved**, given model and current linearization |
| New acquisitions can be certified to compensate rank loss | **Proved**, given assumptions |
| Truncation lowers risk when bias is small enough | **Conditionally proved** |
| PRASC-SOM is faster, more accurate, or larger-basin than tuned direct inversion | **Not established; performance-type method claim currently fails** |
| These safety strategies require SOM | **False** |

## D. Additional assumption-boundary entries (needed to read the ledgers correctly)

- Theorem 1 (phaseless contraction) is [proved] but its mathematical property is statistical
  experiment comparison, **not** an inverse-scattering-specific innovation (line 533).
- Calibration-neutral admission [proved recombination]: the admissible-space statements are
  (3.20)–(3.21), with the complex-linear restriction (3.22); “neutral” does not imply
  sufficient current expressivity, which is exactly why the Theorem-7 bias gate is required.
- Rank–acquisition budget is conditional on: fixed whitened/realified/gauge-reduced tangent
  model, nested current enlargement, full declared shared-map nuisance chart with
  $G=A^TA\succ0$, and acquisition evaluated *after* enlargement. It is a local
  necessary-and-sufficient statement, not a basin theorem.
- Theorem 7 bias-aware risk (3.23) assumes fixed linear experiment, declared model-error set
  $D_rz$, full column rank $B_v$, and the estimator $B_v^\dagger y$; adaptive rank selection
  voids direct post-selection use of that covariance.
- Rank certificates without oracle: exact free-nuisance rank certificates are impossible
  from matrix error bounds alone (counterexample $N(t)=te_1$, Eq. 3.39). Only a
  numerically resolvable rank or a conditional diagnostic may be reported unless structural
  constant rank is certified.

## E. §8/§9 science-status entries relevant to rank/bias/acquisition paper claims

From the five candidate-combination audit (§8.3) and unresolved-problem list (§9):

| Candidate combination / blocker | Status |
|---|---|
| SOM rank as calibration control | viable only after separating $L_{\mathrm{det}},r_{\mathrm{num}},r_{\mathrm{free}}$; SOM-specific policy remains a candidate |
| Nested current nuisance → Loewner loss | theorem holds; generic projection mathematics is not itself new |
| Dual $\rho$ / $B_{\mathrm{vis}}$ gating | usable application combination; cannot alone certify pose correctness or nonlinear basin |
| Low-frequency/low-rank dual continuation | both elements are prior art; candidate novelty is only the certified joint gate |
| Observability-driven acquisition adaptation | idea is prior art; the exact shared-map innovation + rank-loss budget is the narrower theoretical object |
| Physically nonempty safe $(r,k)$ interval | open; needs explicit Helmholtz/Maxwell family and constants (blocks electromagnetic-specific strong theory) |
| Matched-budget computational gain | open; requires E4 fixed protocol plus a direct+same-controller ablation |
| Nonlinear adaptive-gate coverage | open; needs uniform bounded-error regions or independently calibrated coverage bounds |

## F. Verification caveat

The claim ledger in A2.md ends with an explicit note: the old experiments’ numerical and
failure explanations come from uploaded review files; **A2.md did not rerun or independently
reproduce them** (line 2134). The statuses above are A2.md’s documented statuses, not fresh
experimental verification.
