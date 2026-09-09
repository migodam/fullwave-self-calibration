# A2 Rank/Bias/Acquisition Controls: Quantitative Audit and Replication Supplement

Status: written quantitative audit of the already-executed synthetic E1/E2/E3/E5 modules and the optional physical-tangent run in this directory. It adds numeric tables, seed-level uncertainty, baseline comparisons, and an explicit pass/fail/not-yet-run ledger. It introduces no new algorithm, model, or hardware claim. It does not execute E4, does not claim nonlinear superiority, and is not a publication acceptance statement.

## 1. Execution and test suite

- Interpreter: venv Python 3.13.13; numpy 2.5.2; scipy 1.18.1; matplotlib 3.11.1 (env_report.txt).
- Self-runnable assertion tests (pytest not installed): 50/50 passed, all exit codes 0 (results/test_suite_report.txt, results/test_suite_results.json): E1 9/9, E2 11/11, E3 15/15, E5 14/14, physical smoke 1/1.
- Seed/sample ledger (recorded in each JSON's settings): E1 seeds 101-110 x 2000 observations; E2 deterministic exact controls plus 12 fixed tangents (seeds 201-212) giving 36 pair checks; E3 seeds 301-310 x 1000 draws for the Theorem-7 tangent plus a 500-draw adaptive-rank caveat; E5 seeds 401-412 giving 24 rank-budget records plus 20 random innovation-sweep fixtures; physical optional N=8 full aperture, 12 seeds each for phase1 and phase2, declared C1=4 and C2=8 current columns, no oracle rank.
- Runtime note: the only recorded wall time is the physical phase-0 forward (0.0191 s). The synthetic tangent checks are sub-second and were not separately wall-timed in the artifacts.

## 2. Preregistration ledger (protocol vs executed evidence)

| Protocol item | Status | Evidence / note |
|---|---|---|
| Common 1: matched phaseless from the same noisy total-field observations | pass | E1 uses the exact induced noncentral-chi-square/Rice intensity law from the same complex parent observations; no incident-intensity subtraction. |
| Common 2: fixed conventions (time, noise, units, metric, gauge anchors) | pass | e^{-i omega t} physical core; unit-noise whitening; sqrt2 realification; joint gauge anchors tested in E5. |
| Common 3: shared map/trajectory, no independent per-frequency maps | pass | E5 shared-map vs independent-map algebra; physical phase2 uses 4-frequency shared frames. |
| Common 4: separate tuning/test seeds; no oracle rank enters an algorithm | pass | Seed ranges recorded per experiment; physical runner uses declared cutoffs and 'no oracle rank'. |
| Common 5: truth only for evaluation; same nonoracle initial map | pass / N-A | Synthetic tangent checks have no nonlinear inversion; physical runner uses a declared fixed point (alpha=0.2, x=0). |
| Common 6: physical reduced state separate from free-current envelope | pass | Zero envelope score is never labeled physical nonidentifiability (E1 nuisance_phase, E3 b_zero, physical phase1 zero B_v). |
| Common 7: record rank tolerance, absolute singular values, nuisance codimension, uncertainty bounds | pass | E2 saturation/threshold controls; physical declared_C and singular-value lists; CI statistics in results/ci_stats.json. |
| E1 scalar matched models + mismatch control | pass | Four matched models obey contraction; direct-intensity mismatch control reverses ordering by design. |
| E1 physical full-wave scene (fifth minimal model) | not run | E1 artifacts contain scalar models only; the physical tangent run does not compute a coherent-vs-phaseless Fisher comparison. |
| E2 exact nested-rank controls | pass | Loss identity residual 0 (machine precision); rank-drop predictions match in all cases. |
| E2 12 fixed physical linearizations, seeds 201-212 | pass (with note) | Executed in the physical optional phase1 as fixed declared-point tangents, not as nonlinear iterates; T5b 24/24. |
| E2 numerical-rank / saturation / gap-closure controls | pass | Saturation J_x_fro 7.0e-16 (scaled 0.0971); 1e-14/1 singular-value ratio never counted as rank; projector gap closure is a discontinuity, not a rank change. |
| E3 confounding model and eps sweep | pass | Crossover at eps=1; truncation better for eps<1, full model better for eps>1. |
| E3 Theorem-7 tangent MC, seeds 301-310 x 1000 | pass | Variance/bias MC within 3 SE on 10/10 seeds; sample covariance matches inverse information. |
| E3 physical extension full and limited aperture, with/without priors | partial | Full-aperture pose-prior checks pass; limited-aperture tangents are not present in the artifacts. |
| E3 controls: weak-absolute, nearly-parallel ranges, small-residual-large-bias, gauge nulls, prior-only | partial | b_zero (unit retention, zero absolute pose info), small-residual-large-bias, adaptive-rank caveat, and prior checks are present; a separately labeled nearly-parallel map/pose-range fixture is not. |
| E4 primary matched-budget method comparison | not run | Owned by the parent integration workflow; not part of this bounded package. |
| E5 minimal pair, duplicate frame, shared-vs-independent map | pass | J_stack=2 vs sum_J=0 for complementary frames; duplicate frame gives ~1e-15. |
| E5 rank/acquisition budget identity, seeds 401-412 | pass | 24/24 identity records within 1e-12; safety condition I>=L holds 0/24 for random seeded enlargements and is correctly flagged as degradation. |
| E5 physical 8-candidate acquisition-policy comparison | not run | The synthetic fixed-old-model greedy/pair/exhaustive/random comparison DID run (section 6); the physical tangent run does not enumerate 8 physical candidate actions. |
| E5 gauge null and anchor restore | pass | Global rigid gauge (1,-1) persists without anchor; anchor row removes it (min eig 0 -> 1). |
| E5 source-stabilizer (isotropic vs directional) | not run | Recorded as not_yet_run; requires a dedicated physical source-stabilizer tangent experiment. |
| Publication evidence boundary (3D Maxwell / measured array; TGRS imaging gate) | not run | Outside this bounded validation package. |

## 3. E1 matched coherent vs phaseless information

| model | J_coh analytic | J_ph analytic | J_coh empirical (95% CI) | J_ph empirical (95% CI) | max bin-weighted RMSE | contraction pass |
|---|---|---|---|---|---|---|
| phase_only | 2.0 | 0.0 | 1.975 (1.933, 2.017) | 0.0 (0.0, 0.0) | 0.2257 | True |
| nuisance_phase | 0.0 | 0.0 | 0.0 (0.0, 0.0) | 0.0 (0.0, 0.0) | 0.2180 | True |
| scale_family | 1.0 | 1.0 | 1.0003 (0.943, 1.058) | 1.0003 (0.943, 1.058) | 4.83e-11 | True |
| total_field_ref | 2.0 | 0.5258 | 2.003 (1.943, 2.062) | 0.5285 (0.5141, 0.5429) | 0.2032 | True |
| mismatch_intensity_noise | 2.0 | n/a | not MC-run | not MC-run | n/a | False (by design) |

Reading: every empirical-minus-analytic 95% CI contains zero for the matched models; the strict/equality pattern (0<J_ph<J_coh for a known reference; J_ph=0 for phase-only; equality for the scale family) matches the efficient-score measurability prediction. The mismatch control (independent additive intensity noise) has J_mismatch_analytic=1960 >> J_coh_analytic=2, reversing the ordering by design and warning against treating a direct-intensity sensor model as evidence against data processing.

## 4. E2 nested rank, neutral admission, numerical rank

- Exact controls: loss_residual_fro=0 for equality, strict-loss, complete-hiding, and loss-without-rank-drop; predicted vs actual rank drops 0/0, 1/1, 2/2, 0/0.
- Nonmonotone retention rho sequence [0.75, 1.0, 1.55e-15]: map retention is not monotone in the retained subspace.
- Neutral admission: real safe kernel leaves J_x unchanged (safe_loss_fro 2.1e-15; unsafe direction lowers min eig by 1.387). Complex safe subspace S=ker(F) intersect ker(F J_c) has jc-invariance residual 6.3e-15 and is lossless. Utility projection: leading eigenvector of P_S T P_S attains trace 10.6996 vs max over 300 random directions 10.6996.
- Fixed random tangents (n=36): backward-scaled residual mean 0.1056 (95% CI 0.0867-0.1244), max 0.2734 (declared pass threshold 100); T5a loss residual max 4.09e-15; T5b rank-drop match 36/36.
- Numerical-rank controls: singular values [1, 1e-8, 1e-14, 0] give absolute rank 2 at tol 4e-12 and relative rank 2 at tol 1e-12; the 1e-14/1 ratio is roundoff and is never counted as physical rank. Saturation control: J_x_fro=7.0e-16 with backward-scaled residual 0.0971 and max projected singular value 5.8e-16, so full-current rank saturates the nuisance. Projector gap closure at t=0 has projector Fro distance 1.4142 while rank(H(t)) stays 2 for all t: a discontinuity/non-uniqueness event, not an H rank change.

## 5. E3 dual spectra and bias-information tradeoff

- Dual spectra: exact pi/4 and pi/2 fixtures reproduce predicted map/pose spectra and multiplicities (residuals <= 2.2e-16); random nuisance fixture max spectral residual 8.9e-16; rho-undefined directions are excluded rather than assigned 0 or 1; scaling control keeps rho invariant while J_x scales exactly as s^2 (0.008139 / 0.8139 / 81.39 for s=0.1/1/10).
- Pose prior: K_e <= K_eL <= K0 (min eigenvalues 1.2e-16 and 7.0e-18, all within the -1e-10 tolerance); K_eL monotone in Lambda_x; variational form max residual 1.3e-15; augmented principal angles match the closed form (4.4e-16); Lambda=0 recovers K_e (6.1e-16).
- Theorem 7 analytic risk decomposition: variance trace 3.1997, worst bias^2 6.2183, total 9.418. Deterministic bound 5.066 vs worst analytic error 2.494 (satisfied). Monte Carlo (10 seeds x 1000 draws): variance relative error mean 0.0213 (95% CI 0.0122-0.0304, max 0.0411); bias relative error mean 0.0271 (95% CI 0.0148-0.0394, max 0.0535); sample covariance vs J_x^{-1} Frobenius relative mean 0.0576 (95% CI 0.0435-0.0716, max 0.1026); all seeds within 3 SE.
- Confounding sweep (full risk = 1 + eps^-2, truncated risk = 2, crossover eps=1): eps 0.2 -> 26 vs 2 (truncate); 0.5 -> 5 vs 2 (truncate); 1.0 -> 2 vs 2 (tie); 2.0 -> 1.25 vs 2 (full); 5.0 -> 1.04 vs 2 (full).
- Gate failures: b_zero control has map retention rho=1 (maximal) while pose information J_x=0, so a relative-only map gate passes despite zero recoverable pose. Small-residual-large-bias control has fitted residual 0.001 against pose error 1.0, so a residual-only gate passes while the pose is wrong by the full c*=1. Adaptive-rank caveat: selector and estimator share noise; pooled covariance of the selected estimator differs from the selected-model inverse information (Fro difference 2.26).

## 6. E5 shared-map acquisition, budget, policies, gauge

- Shared-map algebra: complementary frames J1=0, J2=0, J_stack=2 > sum_J=0; duplicate frame J_stack_dup ~1e-15. Random innovation formula matches direct stack (max 3.3e-15; 20-fixture sweep worst 1.8e-14). Singular G fallback: naive pinv formula wrong by Fro 1.061 while variational/direct PSD forms match (min eig J_new=8.577). Random 3-frame J_stack - sum_J min eig 0.2552 (Loewner), with an exact-equality nontrivial common-c case.
- Rank/acquisition budget: exact identity J_final - J_original = I_acq - L_rank holds in 24/24 seeded records (max Fro residual 2.34e-14). The first-order safety condition I_acq >= L_rank holds in 0/24 of those random rank enlargements, and in each record min_eig(J_final - J_original) is negative and equals min_eig(I_acq - L_rank) (e.g. -3.755), so the budget correctly flags those enlargements as degrading J rather than certifying safety. Crafted cases: compensating enlargement I-L ~ 0 (neutral); non-compensating half-strength acquisition has min_eig(I-L)=-9.785 with matching J degradation; explicit neutral control L=0 and I>=0.
- Policies (synthetic fixed old model, 8 candidates, K=3): logdet greedy 8.0418, pair lookahead 8.0418, exhaustive 8.1263, random best-of-12 8.1263; greedy optimality gap 0.0845 (greedy is not optimal). Rank-1 Schur-update charges: greedy 21, pair 62, exhaustive 168, random best-of-12 36; full logdet evaluations: 21/34/56/12.
- Non-submodularity: complementary pair violates submodularity with gaps 3.0445 (lambda=0.1) and 7.6014 (lambda=0.001).
- Gauge: global rigid (1,-1) joint null persists for L=2 and L=4 (min eig 1e-15 and 0); an absolute anchor row removes it (min eig 1, null dim 0). Source-stabilizer control: not_yet_run.

## 7. Optional physical-tangent corroboration (secondary)

- Phase 0: physics import OK; A_r 576x9, B_r 576x3; wall 0.0191 s.
- Phase 1 (12 seeds): T5b rank-drop match 24/24; max T5a backward-scaled residual 0.01006 (per-seed max mean 0.005076, max 0.01006; threshold 100); max Theorem-3 residual 3.17e-11 (threshold 1e-6); pose-prior sandwich min eigenvalues -1.83e-17 and -9.75e-19 (tolerance -1e-8). Theorem 7 is NOT evaluable on this set: with declared C1 the residual visible pose matrix B_v is numerically zero on 12/12 seeds (min eig J_x ~ 2e-31), so its full-column-rank assumption fails and the run correctly skipped 12/12 Theorem-7 checks rather than fabricating one.
- Phase 2 (12 seeds, 4-frequency shared frames): max innovation relative residual 4.48e-16; max budget-identity relative residual 6.64e-17; stack>=sum in 12/12 seeds with min eig J_stack - sum J_l 7.51e-5 to 2.11e-4; recorded rank-budget Loewner I>=L count 12/12 for the declared C1->C2 enlargement. All 12 rank-budget steps used the singular fallback path (near-singular physical Gram matrices), consistent with the E5 singular-fallback control.

## 8. Planned baseline comparisons

| Baseline | Result |
|---|---|
| Physical vs free-current envelope | Both satisfy the same algebraic identities. T5a backward-scaled residual: synthetic max 0.2734 vs physical max 0.01006. T5b: 36/36 synthetic, 24/24 physical. Theorem-3 residual: synthetic 8.9e-16 vs physical 3.17e-11. Innovation identity: synthetic 5.6e-15 vs physical 4.48e-16. Key contrast: the synthetic envelopes have visible pose information and Theorem 7 evaluates (10/10 seeds); the physical declared-C1 tangents have zero visible pose information (0/12 full-rank B_v), so Theorem 7 is not evaluable there. |
| Relative-only vs absolute/bias-aware gates | Relative-only fails b_zero (rho=1 with J_x=0). Residual-only fails small-residual-large-bias (residual 0.001 vs pose error 1.0). The bias-aware full-vs-truncated comparison crosses over exactly at eps=1. |
| Shared-map vs independent-map profiling | Complementary frames: J_stack=2 vs sum_J=0. Random 3-frame: min eig(J_stack - sum)=0.2552. Duplicate frame: J_stack ~1e-15. Shared-map compensation criterion separates visible (rank D=1) from hidden (rank D=0) cases. |
| Greedy vs pair lookahead | Both greedy and pair lookahead reach logdet 8.0418 vs exhaustive 8.1263; greedy gap 0.0845. Pair lookahead does not recover optimality on this fixture and costs more rank-1 updates (62 vs 21). Random best-of-12 happens to match exhaustive here (8.1263). |

## 9. Negative controls and threshold sweeps

- E1 direct-intensity mismatch: J=1960 vs J_coh=2, reverses ordering by design.
- E2 small-singular-value trap: [1, 1e-8, 1e-14, 0] -> rank 2 under both absolute (4e-12) and relative (1e-12) cutoffs.
- E2 projector gap closure: Fro 1.4142 jump at t=0 with rank(H) constant 2.
- E3 confounding eps sweep (0.2, 0.5, 1, 2, 5): crossover at 1.
- E3 adaptive-rank post-selection caveat: pooled-vs-majority inverse-info Fro difference 2.26.
- E5 non-submodularity lambda sweep (0.1, 0.001): gaps 3.04 and 7.60.
- E5 saturated new frame: zero cleaned tangent gives zero innovation (6.8e-15).
- E5 singular fallback: naive pinv invalid (Fro diff 1.061); PSD variational/direct form correct.

## 10. Uncertainty and replication summary

Seed-level uncertainty is recorded in results/ci_stats.json and summarized above: E1 empirical Fisher means with 95% CIs across 10 seeds; E2 36 tangent records; E3 10-seed Theorem-7 relative errors with CIs and 3-SE gates; E5 24 rank-budget records; physical 12-seed phase1/phase2 aggregates. No matched-model identity violation beyond its declared backward-error/MC bound was observed. The only directional findings are expected counterexamples (mismatch control, relative-only/residual-only gate failures, non-submodularity, greedy suboptimality, random rank enlargements with I<L).

## 11. Discrepancies and gaps (preserved, not hidden)

1. E1 physical full-wave coherent-vs-phaseless Fisher scene: not run.
2. E3 limited-aperture physical extension and a separately labeled nearly-parallel map/pose-range control: not present.
3. E5 physical 8-candidate acquisition-policy comparison and source-stabilizer control: not run (source-stabilizer recorded as not_yet_run).
4. E5 synthetic random rank enlargements: 0/24 satisfy I>=L and all degrade J; this is a correct budget prediction of degradation, not an identity failure.
5. Physical Theorem 7: skipped 12/12 because declared C1 saturates the visible pose direction (B_v numerically zero).
6. E4 and the publication evidence boundary are outside this package and were not executed.

## 12. Limitations and non-claims

This is a small synthetic/conditional validation package: low-dimensional algebraic tangents, scalar physical forward with N=8 full aperture, and no nonlinear inversion. It does not claim E4 completion, nonlinear superiority over a matched full-wave baseline, SOM exclusivity, universal basin behavior, hardware success, or publication acceptance. Passing these conditional controls is necessary evidence for the rank/bias/acquisition machinery, not sufficient evidence for the full A2 research plan.

## 13. Conclusion

On the executed synthetic examples, the conditional controls hold: E1 information contraction holds in all four matched scalar models with the mismatch control reversing as intended; E2 loss/rank identities hold to machine precision with the numerical-rank traps and gap-closure events recorded; E3 dual spectra, pose-prior extension, Theorem-7 bias-variance decomposition, and confounding crossover all match analytic predictions within Monte Carlo uncertainty; E5 shared-map innovation, rank budget identity, gauge behavior, and non-submodularity match the theorem package, while greedy suboptimality and random rank-enlargement degradation are preserved as counterexamples. The optional physical tangents corroborate T5a/T5b, Theorem 3, prior extension, innovation, and budget identities at N=8 full aperture, but Theorem 7 is not evaluable there because the declared current saturates the pose direction. E4, the physical acquisition-policy comparison, and the source-stabilizer control remain not run.
