# Frozen shape/material stress test of discrepancy weighting

2026-09-08. Registered after the two-scene weighting and correct-model controls,
before generating or inspecting the new cases. This is a factorial development
stress test, not a random-population confirmatory benchmark.

## Falsifiable question

Does the observed rank-one discrepancy-weighting benefit survive a change in
material contrast and boundary shape, rather than depending on the favorable
94–96% error alignment of the original epsilon=2.5 ellipsoid?

Four cases, all new noise/electronics/translation seeds:

| Seed | Known support | Real epsilon |
|---|---|---:|
| 8401 | Ellipsoid | 1.8 |
| 8402 | Ellipsoid | 3.5 |
| 8403 | Box | 1.8 |
| 8404 | Box | 3.5 |

Use existing support dimensions (semiaxes/half-lengths 0.16, 0.11, 0.075 m),
fixed loss 0.03, known plane-wave directions/polarizations, 12 training and 17
held-out receivers. True translation norm is 90 mm, common delay length 0.06,
and illumination electronics follow the existing scene RNG construction.
Only epsilon is replaced after that construction; it is not passed to fitting.
All cases remain known-support, one-material-parameter inversion, not arbitrary
3D imaging or unknown support recovery. No phase-center/coupling claim.

ADDA N64 LDR supplies independent-code data at k=(3,6,9,18). Inversion is N32
CM+RR; discrepancy modes use N16 versus N32. New k18 N96 ADDA fields provide a
separate N64-to-N96 reference sensitivity check on training/held-out receivers.
Always report that difference; values below 1% indicate only a bounded reference
sensitivity check, not proof of monotone continuum convergence. No case may be
discarded or regenerated because that comparison or the inverse outcome is bad.

Noise is proper complex, 30 dB relative to each case's low-band field, with one
fixed absolute sigma across frequencies/methods. One low observation block,
one independent high block and one low-band noisy electronic reference (sigma
0.02, seed+10000) are generated once and saved. Both reference conditions run.

## Methods and budget

For each case/reference condition, fit one low-only pilot from epsilon=2,
all other parameters zero. Continue the SAME low+high observations from that
pilot with: raw weights, isotropic C=I+||delta||^2 I/n, and rank-one
C=I+delta delta^T. The latter two have equal trace inflation and strength one.
delta is the real-whitened coarse/fine predicted high-block difference at the
pilot, including electronics; no ADDA truth or high residual enters it.
All methods keep existing bounds, parameter scales, 35 evaluations and tolerances.
32 fits total: four cases × two references × (one pilot + three continuations).
Rank two is omitted here: it was a declared alternative, not the primary contrast.
No directional or strength tuning, rescue starts or selective additional fits.

Before collection, test each shape adapter's 13-column derivative on N8 and
equivalence to the original ellipsoid adapter; test high-block whitening layout,
shared data and frozen-input provenance. Run serially, save each endpoint before
evaluation, keep all failures and source/input digests, and refuse incompatible
resume. Do not mutate earlier modules or their globals. Cost records separate
data generation, reference refinement, pilot, mode, continuation and evaluation;
no speed claim is made from these development timings.

## Evaluation and decisions

Report each case, not merely a pooled mean: material and receiver-position
errors, held-out low-band sensor phase/field and structural field, raw high-band
fit, reference sensitivity, optimizer status, mode cost and total cost components.
Primary contrast is rank one versus isotropic; raw and low-only delimit simple
alternatives. Mixed outcomes remain mixed, not a scalar score selected afterward.

Post-hoc truth-based error alignment may explain behavior but must be labelled
oracle diagnosis and never used by the estimator. If alignment/benefit fails on
the box or stronger contrast, restrict the claim or change the algorithm in a
new registered development round. These four support/material cases cannot
justify general statistical or novelty claims even if all improve.
