# Far-field calibration nearest-prior audit

2026-09-08. Parent-owned evidence QA, not a global novelty certificate.

## Retrieval and limits

The installed ScholarQA CLI was actually called with two focused queries about
gain/geometry identifiability and near/far-field arrays. All four snippet/paper
operations returned HTTP 429; zero candidates were obtained. The complete error
bundle is `results/farfield_scholarqa.json`. This is retrieval failure, not a
negative prior-art result. No repeated quota-consuming search was launched.

A separate ScholarQA batch record check (`results/farfield_verify.json`)
resolved DOI 10.3390/s26154954 to the Weiss et al. 2026 article (Semantic Scholar
paper ID d1b5851e6b1e8851f3cdf19d07061adadb0cf016). The Comon proceedings URL
did not resolve through that service. Its primary PDF was independently
verified below, but it is not counted as a successful ScholarQA record check.
Bibliographic resolution does not supply verified full-text method evidence.

During the final numerical collection, direct access to the already resolved
Weiss article was attempted via NCBI's public OA endpoint and Europe PMC's
full-text XML endpoint. The browsing tool returned internal/non-retryable URL
errors for both; no article text was obtained. This adds no Tier A evidence and
does not imply that the article lacks open-access full text. No additional
rate-limited broad ScholarQA search or GLM call was made for these failures.

An already identified primary conference paper was independently retrieved from
the official EURASIP proceedings and all four pages were visually read, because
its extracted mathematical text has corrupted font encodings. Direct access
to the previously identified 2026 Sensors paper again returned HTTP 429; its
earlier abstract-level evidence is not upgraded to full-text verification.

## Verified primary paper

P. Comon and L. Deruaz, *Array Self Calibration: Identifiability Issues*,
EUSIPCO 1996, AP.5. Official primary text:
https://www.eurasip.org/Proceedings/Eusipco/1996/paper/ap_5.pdf

**Tier A evidence.** Sections 2–3 use a narrowband source/covariance model with
sensor gain and phase errors, first-order parameter perturbations, and a
second-order range expansion. Equations (9)–(15) explicitly resolve a coupled
bilinear tangent with orthogonal-complement projectors and remaining
ambiguities. Sections 3.1–3.5 distinguish array geometries and near/far fields.
Section 4 discusses maneuvers and calibrated information to remove residual
ambiguities. Hence general near/far-field calibration identifiability,
projection-based separation, and movement/reference-based ambiguity repair
are established prior ideas, not ours.

**Parent comparison / inference.** Our conditional result instead fixes a
world-illuminated extended Maxwell scatterer, permits its full multiple
scattering, differentiates the outgoing complex vector field with a controlled
C1 remainder, and bounds the gain-profiled receiver-translation sensitivity
as receiver distance increases under explicitly separated noise scalings.
This precise statement was not found in these four inspected pages. That is
only a local comparison: the specialization may still be a straightforward
consequence of older field expansions, and other relevant papers remain
uninspected. It cannot carry a strong TAP novelty claim by itself.

## Claim ledger and decision

| Claim | Evidence | Decision |
|---|---|---|
| Geometry and sensor phase can be confounded | Comon–Deruaz §§2–3, Tier A | Cite as prior; no novelty claim |
| Different field regimes and array geometry change identifiability | §§3.1–3.5, Tier A | Cite as prior |
| Maneuvers/references can remove some ambiguities | §4, Tier A | Cite as prior, with their assumptions |
| Our full-wave distance upper bound is globally new | No complete search | Not supported; do not claim |
| Full-wave nuisance-profiled acquisition outperforms alternatives | Four project development scenes only | Hypothesis for frozen testing, not established |

Closest backward leads, not yet full-text evaluated: Rockah–Schultheiss 1987
array shape calibration Parts I (far-field) and II (near-field), and
Weiss–Friedlander 1991 eigenstructure calibration. Their appearance in the
reference list establishes leads only; no theorem is attributed to them here
without retrieval. The present manuscript must distinguish a useful validated
engineering synthesis from a new scientific mechanism.
