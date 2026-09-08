# Bounded SOM-seed adjudication and research decision

2026-09-08. Parent scientific interpretation after all 240 prespecified final
fits completed. This is not a theorem that every SOM method is useless.

## Question actually tested

Does a Twofold spectral seed, followed by the same task enrichment and chart
reuse as a Krylov seed, provide either a ten-percentage-point joint-success
advantage or a 20% conditional median time saving in the specified small
full-wave numerical benchmark? This was fixed before opening the forty scenes.
All six methods share data, nine material coefficients, physical state equations,
frequency strategy and an eight-second decision deadline. Atomic overruns are
charged, not silently accepted as on-time updates. Protocol and sources were
hashed; no outcome-driven tuning or enlargement was performed.

## Complete results

| Method | Joint successes / 40 | Mean total consumed seconds | Accepted reduced updates |
|---|---:|---:|---:|
| Direct GN | 37 | 1.543 | 0 |
| Direct adjoint | 37 | 2.377 | 0 |
| Generic task | 37 | 3.483 | 326 |
| Sensing task | 37 | 3.302 | 330 |
| Twofold task | 37 | 3.186 | 331 |
| Krylov task | 37 | 2.821 | 334 |

No exception rows occurred, so the review-identified exception/audit-stage
ambiguity did not affect these results. Some methods reached the decision
deadline; the largest observed positive overrun was 0.0472 s. These are
nonpreemptive decision budgets, not exact identical consumed time. Runtime
comparisons assume sufficiently stationary execution conditions. One live
process memory sample was approximately 181 MiB RSS; it is not a per-method
peak-memory benchmark.

Primary paired success: neither method has an exclusive success. The
prespecified conservative one-sided 95% upper bound on Twofold minus Krylov
success probability is **0.08810**, below 0.10. This excludes the declared
ten-point benefit in the specified population, subject to scene-sampling
assumptions. It does not prove exact equality of success probabilities.

Primary speed: 37 pairs are jointly successful and objective-matched to the
frozen relative 1e-4 criterion. The conditional median consumed-time ratio
Twofold/Krylov is **1.18616**, with exact order-statistic 95% interval
**[1.17750, 1.19526]**. The exact sign-test p-value against a median ratio at
most 0.8 is **7.28e-12**. Thus the prespecified conditional 20% saving is
excluded; the estimate instead indicates about 18.6% greater median time.
This is a conditional median, not a mean speed claim over failures or other
models. Secondary methods are descriptive, not substituted into the primary
test after viewing outcomes.

## What is and is not decided

- **Decided for this protocol:** the two prespecified advantage routes are
  excluded for this Twofold-initialized physical ROM relative to matched Krylov.
- **Not decided universally:** classical free-current SOM, alternative spectral
  designs, larger/high-contrast materials, unknown gains, different frequency
  continuation or hardware families. This solver comparison intentionally uses
  the same discrete physical model for data and inversion.
- **Mathematical distinction retained:** numerical basis rank is not independent
  current nuisance dimension; prior rank-hiding arguments cannot prove an
  advantage of this physically constrained numerical solver.

## Research-direction decision

Do not keep SOM exclusivity in the title or core contribution. Retain SOM and
Twofold as replaceable numerical seeds and report this bounded negative result.
Do not spend the remaining project only searching for a favorable SOM case.
The main line becomes phase-error attribution and reliable full-wave
self-calibration, with graph/branch/fidelity mechanisms required to demonstrate
actual added value. The latter scientific contribution is **not** established
merely by this negative adjudication. TAP readiness remains open.

Evidence: ROM_FINAL_PROTOCOL.md; results/rom_final.json;
results/rom_final_analysis.json; results/rom_final_source_manifest.json;
results/rom_final_integrity.json; and research/delegated/a3_current_critic/
PROTOCOL_REVIEW.md. Original A2 and all development failures are preserved.
