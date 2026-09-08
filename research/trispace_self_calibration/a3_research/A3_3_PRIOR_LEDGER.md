# Focused verification of supplied A3_3 leads

2026-09-08. ScholarQA CLI Evidence QA; two queries, no systematic-review claim.
Queries: `Idriss Raj modified subspace optimization calibration` and
`Plumlee Bayesian calibration inexact computer models orthogonal discrepancy`.
The collect operation returned status 1: three searches failed HTTP429; the
Plumlee paper search returned five candidates. No repeated paid/provider calls.
An unrelated vision proceedings result was excluded; two recent calibration
papers remain unscreened leads. Retrieved titles alone establish no result.

## Confirmed nearest SOM-calibration source

Idriss and Raj, *Data-Driven Calibration Technique for Quantitative Radar
Imaging*, 2025, arXiv2503.07316v2. ScholarQA verify resolved
`3b55b7fecb4a66fe5d44db2cb294015e17e15f0e`, with no unresolved IDs.
Primary source: https://arxiv.org/html/2503.07316v2
Parent inspected Sections II–IV-A, including equations and numerical discussion.

| Claim | Evidence / tier | Interpretation boundary |
| --- | --- | --- |
| Calibration is complex and indexed by frequency and transmitter | Eq.4, Eq.8 and following indexing clarification; A | It is not merely an amplitude-only method. |
| Measured-data quantitative imaging improves | Table I reports NSE 0.201 without calibration and 0.038 with calibration; IV-A permits separate complex Tx factors; A | Reported by the paper, not independently reproduced here. |
| Gain optimization is iterative | III uses CGD; Algorithm1 updates the simulated field before the gain subproblem; A | Fixed-field separability is a possible implementation adaptation, not automatically a total-runtime advantage. |
| Fitted factors identify electronic rather than geometric error | No inspected passage establishes this | Our tangent-confounding argument does not invalidate their imaging result or diagnose their actual hardware error. |

The supplied A3_3's effective-calibration interpretation is therefore reasonable
as an inference, not a quotation of the authors' conclusion. Do not decompose
a finite complex gain additively into physical components; additive components
are local tangent bookkeeping, not uniquely recoverable quantities.

## Plumlee lead: identity found, theory not yet verified from primary text

*Bayesian Calibration of Inexact Computer Models*, M. Plumlee, 2017;
DOI10.1080/01621459.2016.1211016;
S2 `6b9d374b0660c57e99374455d6c4e5fcd2a9862a`. Subsequent ScholarQA DOI
verification resolved the same title, author and year with no unresolved IDs.
Discovery returned metadata only; publisher withheld abstract. Primary DOI
opening failed in the browser transport. Evidence is Tier C, not Tier A/B.
Do not yet use this record to certify the precise orthogonality construction
or to justify a new EM algorithm. The supplied text remains a research lead.

The existing DISCREPANCY_PRIOR_LEDGER.md already confirms substantial overlap
with approximation-error methods. Neither the new name nor this incomplete
search establishes global novelty.
