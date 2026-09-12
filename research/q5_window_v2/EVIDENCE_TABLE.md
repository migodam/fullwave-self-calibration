# Theorem / assumptions / proof / evidence ledger

This ledger overrides any pre-existing unsourced V2 branch narrative. The immutable A5/V1 parent corrections retain priority for inherited claims. Detailed mathematical derivations are in MANUSCRIPT_SECTIONS.md; accessible closest-prior formulas are in PRIOR_ART.md.

| ID | Proposition / exact scope | Proof or executed evidence | Verdict |
|---|---|---|---|
| T1 | Coordinatewise finite-error fiber criterion; rectangular task range; specified error sets | Midrange necessity/sufficiency proof, Section 2 | Proved general fact; not novelty |
| T2 | Finite range chord for projected electric-dipole tensor, class D | Determinant identity; 90 exact rational design boxes covering nominal kR in [.6,1.5] | Proved for specified relative-error budget |
| T3 | Electric-only material ambiguity, known direction and same range, class D | Exact gain compensation, rational annulus enclosure; materials 2 and 2.3 | Exact restricted counterexample; not class C |
| T4 | Class-M amplitude monotonicity | 5500 exact rational closed boxes; analytic series tails; outward exports | Computer-assisted proof |
| T5 | Class-M E/M ratio finite inverse | Positive whole-domain real-ratio derivative; deterministic quotient error bound | Conditional finite proof, known geometry/resolved modes |
| T6 | Reference readout with bounded leakage/drift | Exact rational budget, 16 repeats/channel, probability >=.997 | Sufficient conditional specification; apparatus unverified |
| T7 | Class-M size error allocation | 550 material boxes times an entire ka band of relative halfwidth 50 ppm | Computer-assisted bound; other constitutive assumptions exact |
| T8 | Class-C continuum forward coercivity | Static nonpositive Fourier multiplier, dynamic kernel bounds, two-ball Schur estimates, exact rational inequalities | Uniform inverse norm <=500/43 |
| T9 | Class-C forward residual-to-data transfer | Receiver distance >=.503 m; receiver norm <6; gain bound | Error <88 times actual continuous residual; residual unknown |
| T10 | Sharp finite amplitude boundary | Whole-domain strict concavity; other-channel overlap; explicit finite competing pairs | Necessary/sufficient for amplitude experiment only |
| T11 | Material-window splitting at fixed raw amplitude error | Monotone exact tangential radial factor and two-sided rational cutoff enclosure | Two distinct amplitude task cutoffs, about 26.7393 and 6.56835 |
| A1 | Standard numerical multisphere VSW solver | Interface, plane-wave reconstruction, quadrature and order diagnostics | Runnable numerical implementation, no solver novelty |
| E1 | Independent DDA-to-VSW recovery, new registered 12 cases | All 48 fits completed; 3 starts each; raw observations and interrupted checkpoint preserved | Reference 7/12, other three methods 2/12; development only |
| E2 | Attribution oracle diagnostics | Protocol frozen after E1 inspection and before oracle execution | Post-hoc: true gain 7/12; true shift 2/12; both 7/12 |
| G1 | Complete class-C finite material separation / physical competing worlds | No proved interval field oracle or actual continuum residual enclosure | UNRESOLVED; window_cover refuses certification |
| G2 | Physical antenna/reference implementation meets the budget | No apparatus, full aperture/probe/drift certificate or calibrated power conversion | UNRESOLVED |
| G3 | TAP-level algorithmic novelty / final advantage | Closest-prior overlap plus unread direct neighbors; no locked final campaign | UNRESOLVED; not submission-ready |

## Numerical evidence and uncertainty conventions

Class M: sigma is the complex standard deviation, not each real quadrature standard deviation. The inherited real amplitude-reference event and the separately specified measurable complex-reference event are different models. In the final hardware specification, within-channel independent repeats reduce variance by 16; a union bound over three radial complex-noise events does not require independence between channels. Common biases do not average away.

Class C: each scene uses sigma=0.01 times the RMS of its own noiseless gain-weighted field. Methods share that scene's identical observations. This supports paired recovery diagnostics, but is not a common-covariance two-world statistical proof. Reference sigma is .01 complex. Fixed/random acquisitions add one complex scalar; their physical cost is not proved equal to the reference path's cost.

All local Jacobians, solve residuals, truncation changes and cross-solver differences are diagnostics unless a separate continuum proof is provided. None is converted to a model-error radius here.

## Failure preservation

- V1's failed two-cycle selector is not modified, retrained or promoted.
- The coarse single-box class-D range certificate failed; its result remains alongside the passing refined cover.
- The new recovery run was interrupted and resumed from exact stored observations; both the incomplete checkpoint and final result remain.
- All five reference failures remain failures; sensor-only residual is never the task criterion.
- DDA grid refinement is not assumed monotonically accurate.
- The old branch report was inspected as an unverified input. Its finite-secant statement was subsequently proved and audited locally; its unrelated unverified content is not automatically accepted.

## Implementation budget and replication

The reference sufficient specification uses 48 complex modal/reference samples, not necessarily 48 raw antenna positions: synthesizing one modal sample can require an aperture scan. Tested numerical recovery uses 312 coarse and 526 finer DDA cells, a 60-coefficient two-sphere L=3 inverse, 144 base complex observations and one optional scalar. Three starts per method, at most 80 optimizer function evaluations per start. Function-evaluation counters are not identical to the number of Maxwell solves because numerical Jacobians and coefficient caching are involved.

Observed main-recovery runtime was approximately 52.83 seconds across execution/resume sessions, maximum process RSS about 283 MiB. These are environment-specific logs, not hardware throughput guarantees or per-method exclusive memory measurements. Reproduction on a CPU needs the pinned NumPy/SciPy versions; there is no GPU requirement. The original raw arrays are in the delivered evidence archive; lightweight repository summaries carry their cryptographic hashes and all 48 scene/method outcomes.
