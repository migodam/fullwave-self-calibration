# A4 literature, originality firewall, and narrative decision

**Status:** targeted nearest-prior audit, not an exhaustive systematic review. Fourteen primary articles/manuscripts were accessed across the A4 work; the depth differs. Full-text retrieval is not equivalent to auditing every proof. Original SOM/Twofold/phaseless full papers and missing A2 audits were not obtained in this execution. No absence-of-search-results argument establishes novelty.

## Evidence tiers

F = all short-paper pages visually inspected; K = primary full text available and key equations/sections inspected, not every proof replicated; P = primary text retrieved but the relevant proof-level comparison remains partial; M = metadata, supplied A3 discussion, or limited book excerpts only.

| Record | Tier | What was actually inspected | Effect on A4 |
|---|---|---|---|
| Comon--Deruaz 1996 [1] | F | Four PDF pages, local ambiguity equations, far/near cases and remedies | Generic self-calibration, projection, maneuvers and references are prior |
| Weiss--Friedlander 1989 [2] | F | Nine PDF pages, likelihood, scoring algorithm, local accuracy and simulations | Joint geometry/source optimization and local information are not new |
| Rockah--Schultheiss I [3] | P | Primary 14-page text; complete assumption-to-proof mapping unfinished | Far-field array-identifiability priority not claimed |
| Rockah--Schultheiss II [4] | P | Primary 12-page text; full estimator/proof reproduction unfinished | Near-field calibration is not itself new |
| Epp--Janz 2013 [5] | K | Field model and normalized magnetic-field gradient equations (9)-(15), ambiguity and sensitivity discussion | Amplitude-independent EM ranging is prior; strongest physical neighbor |
| Torres et al. 2021 [6] | K | Electromagnetic localization model, aperture/CRB framework, assumptions | EM-specific Fisher information is prior |
| Bellomo et al. 2014 [7] | K | Multipole incident-field calibration, effective-height and scattering-parameter equations, inversion context | Antenna-aware diffraction tomography already exists |
| Idriss--Raj 2025 [8] | K | Primary multiframe model and complex calibration factors, MFSOM/experimental sections | Calibration plus SOM cannot be a core novelty claim |
| Kaipio et al. 2019 [9] | K | Approximation-error construction, mean/covariance and Born-surrogate inverse problem | Sampled approximation error must be a baseline |
| Calvetti--Somersalo 2026 [10] | K | Error subspaces, priorsketching/BAE connections and signal suppression | Covariance-to-projection algebra cannot carry priority |
| Beutel et al. [11] | K | T-matrix/multipole formulation and layered-sphere implementation context | Modal diagonalization is established physics; software not executed here |
| Li--Lee--Bresler [12] | K | Bilinear blind calibration formulation and sample-complexity assumptions | Generic gain identifiability is not the EM theorem |
| Thyagarajan et al. [13] | K | Gain invariants and closure quantities | Gain-invariant dimension is not new geometry information |
| Huang et al. [14] | P | Primary source-receiver extension formulation and motivation | Extended FWI is a strong antecedent, not identical to receiver-position recovery |
| Chen/Zhong--Chen/Pan/book [15]-[18] | M | Supplied A3 treatment and limited Library search | No new SOM priority or family-wide refutation claimed |

## Paper-by-paper lessons and closest-prior comparisons

### [1] Comon--Deruaz
Their physical problem jointly includes sensor geometry and gain/phase deviations in an unknown source field. Local equations reveal ambiguity dimensions and how additional knowledge removes them. The conceptual sequence is failure, identifiable combinations, then a remedy; not an optimizer success plot. For A4, nuisance projection, near/far distinctions and maneuver/reference ideas belong to prior art. A possible increment is the exact vector-Maxwell radial/tangential spectrum after radial-material profiling under a different acquisition, not a renamed array-rank argument. Reviewers will require this distinction in the introduction, not hidden in an appendix.

### [2] Weiss--Friedlander
Unknown source locations and array shape enter one maximum-likelihood calibration problem, with a practical iterative estimation method and local accuracy analysis. The lesson is to connect assumptions to an executable estimator and numerical behavior. A4 cannot claim joint optimization or local Fisher profiling. Its proposed physical increment must survive comparison with their source-field model, and a small local covariance cannot be sold as global branch coverage. Our unresolved-cell mechanism is supporting machinery, not a retroactive guarantee for their or our optimizer.

### [3] Rockah--Schultheiss I
The far-field calibration setting is a foundational identifiability and accuracy neighbor. A4 must specify its externally anchored scatterer rather than imply absolute joint localization of arbitrary objects and sensors. The present audit has not completed a theorem-by-theorem correspondence for every array configuration, so no universal acquisition count is attributed to this paper and no precedence claim is cleared. Its role is to block broad originality claims until that comparison is completed.

### [4] Rockah--Schultheiss II
Near-field sources and estimator implementation directly precede the proposed use of near-field diversity. A4's potential distinction is not spherical phase curvature alone: it uses a longitudinal/transverse Maxwell field ratio under three controlled electric modes and a scalar material response. Different source knowledge, gain structure and noise normalization matter. The exact proof-level overlap remains an explicit literature gap. A reviewer can reasonably reject a broad near-field self-calibration novelty claim on this basis alone.

### [5] Epp--Janz
The paper recovers information about an arbitrary time-varying electric dipole from electromagnetic field spectral data, including normalized magnetic-field derivatives. Equations (9)-(15) eliminate unknown dipole factors and recover range; the paper also discusses ambiguous or inconsistent observations. This is a serious antecedent to amplitude-independent ranging, not a peripheral citation. A4 uses a three-illumination electric response tensor of an unknown radial scatterer, rather than spatial magnetic-field gradients. The exact nuisance-profiled spectrum and fixed-SNR optimum are independently derived here, but their research novelty still requires a broader source-localization review. The tempting claim 'near-field ratios locate without knowing material' is too broad.

### [6] Torres et al.
An electromagnetic received-field model is used to derive near-field localization bounds, making aperture and polarization assumptions explicit. This is a model-to-information-to-design narrative. A4 does not invent electromagnetic CRBs. Its gain/material assumptions and retained reactive/longitudinal terms must be compared explicitly rather than assuming all prior near-field models omit them. We retain the restricted comparison and do not claim that A4 dominates their estimator or covers their entire unknown-orientation setting.

### [7] Bellomo et al.
Antenna modeling and incident-field calibration are integrated into quantitative limited-aspect microwave diffraction tomography. Key multipole and effective-height/scattering-parameter equations make the calibration physically testable, with nonlinear inversion and experimental data. The lesson for TAP is that an ideal illumination assumption must eventually become a measured field and an error budget. A4 currently lacks that bridge for pure electric-l=1 modes. Our numerical calibration of a tensor in an ideal homogeneous background does not replace their hardware evidence or establish phase-center ground truth.

### [8] Idriss--Raj
Complex calibration factors are estimated within quantitative radar imaging and a multiframe SOM framework. Consequently, 'self-calibration plus SOM' is not a defensible new contribution. A4's question is narrower attribution and reliability: which changes can be assigned to geometry rather than nuisance or model error? Neither the existence of factors nor favorable measured reconstruction automatically answers that. The prior should be acknowledged in the main text; our SOM history belongs in a short provenance note, not a competing algorithm novelty claim.

### [9] Kaipio et al.
Approximation error between a more accurate scattering model and a Born surrogate is sampled and incorporated statistically. This motivates the full sampled-error and conditional-linear baselines in A4. Error mean and covariance are part of the inference model, not an electromagnetic error certificate by themselves. The narrative lesson is to define the model reduction/error family before evaluating reconstruction, and to retain contrasts where compensation fails. A4's more elaborate risk selector did not beat the direct fine solve, so no mature replacement for this framework is claimed.

### [10] Calvetti--Somersalo
This preprint connects prior-informed subspace ideas with Bayesian approximation error and explains how nuisance/error directions interact with signal. Projection or covariance inflation cannot be recast as a new theory merely by calling it phase attribution. The useful A4 question is whether a specific electromagnetic condition tells us which physical intervention is necessary. The work is cited as a preprint, not represented as a peer-reviewed TAP result. A4's numerical selector remains an application of these mature tools without a demonstrated incremental advantage.

### [11] Beutel et al.
T-matrix representations express scattering through incoming/outgoing modal coefficients, including sphere and multiple-scattering constructions. A4 uses the established rotational invariance of a radial target to reduce material dependence to one electric-l=1 coefficient. This is a physical starting point, not a discovered diagonalization theorem. The independent A4 implementation uses analytic Mie coefficients and an in-house dense DDA solver; it does not claim to have executed treams or to have compared against its published performance.

### [12] Li--Lee--Bresler
Blind gain/phase calibration is formulated as a structured bilinear identification problem with sample-complexity conditions. A4 should not mistake enough data dimensions for recovery of physical geometry after allowing a distributed material family. Its tensor positive condition specifies how Maxwell field structure produces the surviving directions. Generic bilinear guarantees remain prior methodology, and their assumptions cannot be imported into the nonlinear scattering problem without a reduction that preserves all nuisance freedoms.

### [13] Thyagarajan et al.
Interferometric gain-invariant quantities and closure structures are treated through an Abelian gauge formulation. This is a direct originality firewall against gain-graph or closure repackaging. A4 needs neither a new name for the cycle space nor a claim that a nonzero invariant count proves geometry identifiability. Any gain-invariant acquisition must still survive material profiling and model error. Closure material is therefore supporting context or removed from the short manuscript.

### [14] Huang et al.
Source-receiver extension enlarges a full-wave inversion search space to mitigate optimization difficulty, with penalties enforcing physical consistency. It is relevant to branch/cycle-skipping discussion but is not automatically an unknown physical receiver relocation algorithm. A4 cannot claim that rejecting uncovered cells solves general FWI nonconvexity. The full proof and comparative implementation were not reproduced; this is a recorded overlap boundary, not a statement that the prior has been defeated.

## Three narratives, one retained candidate

| Narrative | Strongest content now | Decisive weakness | Disposition |
|---|---|---|---|
| A: attribution/controller | Correct sandwich risk; nine implemented controls; explicit action interface | Same endpoint as fine solve with higher cost; no complete six-action closed loop | Supporting negative result; no main contribution |
| B: EM acquisition/identifiability | Exact nuisance-profiled spectrum; fixed-SNR range window; exact alias and dither proof | Restricted radial target/modal acquisition; near-field dipole prior overlap not cleared | Only retained main-contribution candidate |
| C: branch reliability | Analytic cell bound and retained unresolved state | Generic set-membership machinery; false acceptance under invalid discrepancy family; no win over strong acquisition baselines | Supporting mechanism |

**Final positioning:** B, not an asserted B+A method paper. A defensible manuscript is currently a restricted electromagnetic measurement-design/limits working paper. Originality and experimental sufficiency remain open; changing the paper class does not automatically make it submission-ready.

## Core/supporting/appendix/remove

Core: modal target/acquisition assumptions, exact spectrum, two limiting degeneracies, global sign ambiguity and dither consequence. Supporting: bounded coverage and electronics-reference distinction. Appendix: all counterexamples, nine-method failed controller comparison, DDA model-class stress, full cost and correction audit. Remove from core: SOM/ROM history, rank manifolds, rho, Woodbury/Schur algebra, generic closure counts, standalone multifidelity story, and unlabelled Fresnel antenna-recovery claims.

## Unclosed search and reading requirements

No claim of a completed TAP/AWPL/TMTT/OJAP-wide systematic search, all radar-autofocus variants, antenna phase-center theory, every modern conditional-error approach, or every dipole localization paper. Original A2 audits/questions and several requested original papers were unavailable. These are not reasons to invent priority: the strongest B novelty remains provisional. Bibliographic details for [18] must be completed or that item removed before submission.

See REFERENCES.md for exact records and links.
