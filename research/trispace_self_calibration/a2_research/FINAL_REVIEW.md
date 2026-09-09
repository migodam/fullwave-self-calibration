# Final scientific review and critic disposition

Date: 2026-09-06. Final authority: parent Codex scientific review. Worker reviews
are evidence inputs, not acceptance decisions. Original reports/results remain.

## Decision

Deliverable: a complete **research draft and reproducible bounded study**.
Publication acceptance: **NOT READY for TAP or TGRS**. The conditional theory
survives the checks, and exploratory direct coherent calibration is useful in
the tested image family. A useful SOM-specific reduced solver, high-frequency
admission/cost advantage, covered phase-branch selection, and measured unknown
antenna recovery remain unmet. No additional 3D work was authorized in this
cycle's locked scope; it is an alternative future validation route, not a
required completed deliverable.

## Independent critic: item-by-item response

Source: `research/delegated/a2_critic/CRITIC.md` (F1--F10).

| Item | Parent judgment | Implemented change / remaining gap |
|---|---|---|
| F1: no reduced iterations | Accept, severe | Abstract and Section 7.3 now lead with zero accepted reduced steps in 2,400 runs; no method-performance acceptance claim |
| F2: cost/schedule invalidates comparison | Accept | Explain 158.67-unit full-band evaluation and skinny-product overcounting; retain failed frozen endpoint; fresh cost protocol and test scenes required |
| F3: bounds too loose | Accept | Publish state/derivative bound ratios by frequency; high-frequency admission remains open, not relabeled as success |
| F4: novelty/prior overlap | Accept overlap warning, qualify reviewer wording | No first/SOM-exclusive claim; existing calibration and source-extension overlap explicit. Source-receiver extension is adjacent, not evidence of identical hardware-pose estimation. Passive VIE-bound priority search remains incomplete |
| F5: vacuous gates | Accept | Analysis narrative says VACUOUS for coincident-output noninferiority; original numerical booleans remain traceable but count as no scientific evidence |
| F6: efficient physical Fisher not tested | Closed for bounded test | Parent square-root elimination of nine material directions in the ten existing physical Fisher matrices; all PSD differences, minimum eigenvalue 0.00713--0.01064; no new-scene claim |
| F7: measured geometry not validated | Accept absence of antenna truth; reject original fit interpretation | Parent discovered convention/passivity and test-gain leakage errors in worker pilot. Initial fit withdrawn. Corrected post-hoc model pilot documented below; object center is not antenna calibration |
| F8: scaled identity residual unclear | Accept | State backward-error threshold 100 alongside maximum 0.274; not 27.4% relative error |
| F9: certificate title overstates | Accept | Title becomes model-conditional state error bounds; floating-point vs interval and parameter/model confidence distinctions explicit |
| F10: policy completion vs success | Accept | Report 240/240 final-data policy completions for four methods vs 0/240 fixed rank at each budget, while all recovery successes are zero |

The reviewer's statement that all E4 estimates remain at initialization is true
at budget 200, not at 800 for direct/intensity. The reviewer's use of
"epsilon non-identifiability" is not supported by band-dependent fits of an
incorrect model; neither that inference nor its old measured numbers is adopted.
An oracle known-pose optimization result is a reference, not a guaranteed upper
performance bound; the paper figure removes that original plotting label.

## Measured-data audit and correction

Accepted data: canonical 2001 Fresnel `dielTM_dec8f.exp`, 14,112 x 7 finite rows,
all 36 x 49 x 8 index combinations. SHA-256:
`476cc9d1cfc98797545ab4adf69302dc5aeb45848a24cf8d7b3d222940cc79eb`.
Primary descriptor defines source/receiver angular indices, radii, GHz column,
and total/incident real and imaginary fields. The source is a horn, not an
ideal point source. Unpublished exact object-center direction is an unknown
object parameter, not automatically a gauge with fixed sensor coordinates.

Rejected initial pilot: `research/delegated/a2_measured/model.py`, results_model,
and its synthetic object-location recovery are NOT accepted antenna evidence.

- With outgoing H1, consistent convention is exp(-i omega t), plane phase
  exp(+i k u.r), passive Im(epsilon)>0. The worker mixed H1, incoming plane
  phase, and negative imaginary permittivity, then used the opposite flux sign.
- Independent parent flux checks give outgoing H1 positive power; integrated
  lossless flux -2.2e-15, positive-loss -2.792, negative-loss +4.908 (a common
  positive physical factor omitted). Old incident expansion has maximum error
  1.844 against the specified direction; corrected expansion error 7.6e-16.
- `run_fits` re-profiled complex gains on the test views. Its "held-out"
  numbers are not out-of-sample predictions with a fixed fitted model.

Corrected implementation: `measured_corrected_pilot.py`. It leaves the original
worker files intact, fixes the wave convention, fits real epsilon only, uses
whole-model conjugation to match the primary descriptor's exp(+i omega t)
convention (p. 1570), cross-checks it against training incident data, and freezes training gains
before test prediction. Four fixed starts give the same optimum (7--10 optimizer
iterations each). Training/test normalized squared scattering residuals are
0.0216179/0.0216197. Object center (1.3224,26.0638) mm, epsilon 3.4270. These are
post-hoc descriptive fits, not confidence bounds or known geometry recovery.
The incident point-source approximation itself fits poorly, especially at high
frequency; the better matching sign does not validate the actual horn model.
The time convention itself is explicitly stated in the primary descriptor;
the error was the worker's mixed implementation, not absent metadata. The exact cylinder
boundary/scattering series does not make its plane illumination exact.

No measured antenna-offset truth, clock truth, noise variance law, or blind
calibration comparison exists in this dataset. These gaps remain. The corrected
fit contradicts the worker's blanket high-frequency-failure narrative; the
remaining modest parameter discrepancy cannot diagnose its unique cause.

## Evidence state after revision

**Accepted conditional mathematics:** typed nuisance vs numerical rank; efficient
information contraction; normalized dual spectrum; rank loss and complex-safe
admission; fixed-linear bias risk; shared-map innovation; passive discretization
bound and residual/tangent propagation; conditional descent (not all hypotheses
verified in the constrained finite-budget solver).

**Bounded executed evidence:** E1/E2/E3/E5 algebraic checks and physical tangents;
passive-bound enclosure/refusal; residual-enrichment reference probes; E4 record
integrity and cost failure; 48 larger-map test runs and six mismatch controls;
public data ingestion and the corrected cylinder model-adequacy pilot.

**Unmet:** physical E3 bias-risk MC non-vacuity (0/20 evaluable); E3 near-parallel
physical fixture; useful high-frequency reduced steps; classical-SOM comparator
and block-Krylov/ROM competition; source-stabilizer controls; nonlinear active
acquisition; covered phase-branch gate; known-offset measured validation;
general moving-array/3D/automotive evidence; global priority of the passive bound.

## Implementation notes affecting reproduction

The frozen E4 modules are unchanged after `frozen_parent.json`; do not repair
their unused legacy `prepare_scene_impl` command path in place. The canonical
runner is `study.py final`, which uses its own validated scene construction.
The old helper has unresolved names and is not a supported entry point.

The physical E1 replication has an unused small-lambda shortcut returning 0.5
for a noncentrality-score variance whose correct limit is 0.25. Every recorded
physical row used quadrature (lambda >= 2.7e-4), so the published matrices and
parent Schur checks are unaffected. Its metadata phrase about a nonzero
mean-shift intensity Fisher exactly at zero mean is also not adopted. The
shortcut is a known unexercised limitation of that historical script, not a
validated branch. Never use it for new zero-amplitude experiments without a
separately versioned correction and zero-limit tests.
