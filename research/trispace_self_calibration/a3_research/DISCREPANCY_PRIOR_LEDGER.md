# Discrepancy method prior-evidence checkpoint

2026-09-08. Focused evidence QA, not a completed review or novelty clearance.

Question: which approximation-error and covariance-weighting mechanisms already
cover coarse/fine electromagnetic inverse models, and what incremental result
would be needed beyond those mechanisms?

The ScholarQA discovery request `electromagnetic inverse scattering Bayesian
approximation error` returned HTTP 429 and zero candidates (documented in
DISCREPANCY_WEIGHTING_INTERPRETATION.md). No absence-of-prior inference follows.
A subsequent targeted identifier verification succeeded:

| Record | Evidence | Disposition |
|---|---|---|
| J. Kaipio and E. Somersalo, *Statistical inverse problems: discretization, model reduction and inverse crimes*, 2007, DOI 10.1016/j.cam.2005.09.027; Semantic Scholar 9c9b363e30a248c059c62518f92183eb1f462c81 | Tier C: authors/title/year/DOI resolved; abstract withheld | Relevant foundation to inspect, not yet equation-level or abstract-level support for a specific method claim |
| Candidate DOI 10.1016/j.cam.2006.01.037 | Resolved to Casian and Kodama's unrelated Toda-lattice paper | Wrong candidate rejected; do not cite as inverse-problem methodology |

Primary publisher retrieval of the verified inverse-problem paper failed:
DOI redirect to Elsevier returned a tool fetch error, and the ScienceDirect
article S0377042705007296 returned HTTP 403. Full text was not read. The record
is therefore not promoted beyond Tier C and is not added as a substantive
method citation to the manuscript yet.

The local Woodbury/penalized-error, realification and fixed-design covariance
identities are explicitly standard mathematics. Empirical benefit of our frozen
pilot mode is not itself proof of methodological originality. The pending audit
must compare purpose (calibration versus imaging), mechanism (deterministic
pilot mode versus prior-sampled/error-adaptive models), and evaluation (independent
Maxwell data with uncertain geometry/electronics). Those are questions to examine,
not asserted distinctions from unread papers.

## Direct methodological overlap now verified

Calvetti and Somersalo, *Spotlight, priorsketching and Bayesian approximation
error paradigms*, arXiv:2604.26254v1 (2026), Semantic Scholar
dd72596adbacd1a0c39e1be4fd58a8e4d5e7df32, resolved by ScholarQA verification.
Identified through one 20-record forward-citation page of the Kaipio paper;
browser API retrieval failed, direct Semantic Scholar endpoint succeeded.
No exhaustive citation screening is claimed.

Tier A: primary HTML sections 2.1–3.4 and the start of 4.1 inspected. Equations
10, 16, 24 and 34 establish covariance weighting, truncated projection, their
limiting relationship and low-rank inverse application. Section 2.2 discusses
loss of useful signal when target and nuisance ranges nearly align.
Source: https://arxiv.org/html/2604.26254v1

Disposition: generic directional weighting/projection is not our original
contribution. A low-pilot electromagnetic implementation is a candidate
adaptation; incremental novelty and useful operating conditions remain unproved.
Next controlled question is target-preserving model-fidelity/acquisition choice,
not another name for this algebra. Existing positive fits remain valid evidence
of the tested adaptation, not evidence of priority.

## Inverse-scattering antecedent

Kaipio, Huttunen, Luostari, Lähivaara and Monk (2019), *A Bayesian approach
to improving the Born approximation for inverse scattering with high-contrast
materials*, DOI 10.1088/1361-6420/ab15f3, arXiv:1901.00909v2,
Semantic Scholar 9102acc38d0ce56d5e1dac3fbe365b7fdda0d0e1.
Title-match retrieval and separate ScholarQA record verification succeeded.

Tier A: primary introduction and Sections 2–4.2 inspected at
https://arxiv.org/html/1901.00909 . It treats 2D scalar far-field material
inversion using Born-model error statistics, including error/material
cross-correlation. Training uses 3000 forward samples; computational costs are
separated into offline and online work. The reported acquisition geometry is
prescribed, not the joint geometry/electronics problem posed here.

Disposition: approximation-error compensation in inverse scattering is prior
art. Our full-wave, low-pilot geometry/electronics adaptation differs in setup,
but those differences alone do not prove a nontrivial new contribution or
superiority. A strong sampled-error baseline and honestly charged costs remain
necessary before such claims.
