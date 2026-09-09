# Registered strong sampled-error likelihood baselines

2026-09-08. Registered after the shape/material stress outcomes and before
sampling or computing these baselines. These are adaptations of established
approximation-error methodology, not proposed original algorithms. Primary
equation-level antecedents and inspection scope are in DISCREPANCY_PRIOR_LEDGER.md.

## Purpose and frozen setup

Test whether the single low-pilot discrepancy mode offers a useful accuracy/cost
tradeoff against stronger sampled-error models. Use the four existing stress
cases and both reference conditions with unchanged observations, low pilots,
physical material/geometry/electronic freedom, bounds and held-out metrics.
This remains development, with no population significance or novelty claim.

For each known support (ellipsoid and box), generate 24 independent parameter
draws using RNG seeds 8801 and 8802 respectively, never a test truth or fitted
endpoint. Draw independent uniform coordinates:

- real epsilon in [1.3,4.5], loss fixed0.03;
- each receiver translation in [-0.12,0.12] m;
- delay length in [-0.1,0.1] m;
- four log amplitudes in [-0.15,0.15], four phases in [-0.3,0.3].

At k18 compute ADDA N64 LDR and inverse N32 CM+RR predictions for each draw,
including its electronics. Error samples are accurate-minus-inverse, unlike
the unsigned covariance proxy from a coarse/fine pilot. Save every parameter,
field difference, executable/source hash, true residual check and cost. The
prior is a declared additional modeling assumption, not learned from test truth.
It is not an assertion that the experimental scenes were iid draws from it.

## Two baselines

Realify unwhitened high-block errors e_i and scaled parameters u_i=theta_i/SCALE.
Use sample means and covariance with denominator23. Construct:

1. `sampled_eem`: fixed error mean and covariance, independent of theta.
2. `sampled_conditional`: error mean m_e+L(u-m_u), with
   L=C_eu(C_uu+lambda I)^-1 and residual covariance
   C_ee-L C_ue. Set lambda=1e-6 trace(C_uu)/13 in advance, solely for numerical
   stabilization. Require the residual covariance to be PSD within numerical
   tolerance; do not silently repair a substantive negative eigenvalue.

For each experimental sigma, add measurement covariance (sigma^2/2)I and use
its symmetric inverse square root on the high residual after subtracting the
appropriate error mean. The conditional high Jacobian includes L/SCALE with
the same residual sign convention. Low and reference residuals retain their
original weights. Both objectives are error-adjusted likelihood fits with bounds;
no additional material MAP regularizer is introduced. Consequently this is not
a claim of reproducing a prior paper's entire posterior estimator.

Start from each saved low pilot; same35 evaluation cap and tolerances. Sixteen
new inverse fits: four cases × two references × two error models. Compare all
with prior low-only/raw/isotropic/rank1 outcomes. Do not retune after results.

## Checks, costs and failure handling

Before training or fitting: check covariance/conditional algebra, whitening and
the 13-column corrected residual Jacobian on a small physical model with a
synthetic, frozen error model. Check source and immutable input hashes. Sample
generation and inverse runs are serial CPU, no parallel timing. Save training
samples individually and accepted estimates before evaluation; retain failures,
stages, interruption histories and incompatible-resume guards.

Report offline training cost per shape, first-use total and reuse cost separately.
Do not divide training cost by a favorable invented number of future inversions.
New data generation is not free; low-pilot cost and online correction overhead
also count. Historical cost sums are not fresh timing measurements. Low-band
held-out sensor and structural metrics use the actual unchanged inverse model;
the high correction is not silently applied to untrained held-out channels.

A baseline win is retained even if it removes the single-mode advantage. A
single-mode win here is conditional on this training prior and sample count,
not a refutation of approximation-error methods as a class. A later sample-
size sensitivity test must be registered separately, not stopped when a desired
ordering appears. No trained error covariance is a verified continuum enclosure.
