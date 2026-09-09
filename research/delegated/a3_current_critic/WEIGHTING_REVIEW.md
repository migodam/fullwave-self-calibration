# Bounded review of discrepancy weighting mathematics

Date: 2026-09-08. Read only DISCREPANCY_WEIGHTING_PROTOCOL.md and DISCREPANCY_WEIGHTING_INTERPRETATION.md. No implementation inspection, experiments, APIs or literature checks. Parent retains final mathematical and scientific judgment.

## Disposition

The weighting algebra is correct under fixed-weight linear-model assumptions, and the protocol correctly calls the coarse/fine difference a proxy rather than a known physical covariance or continuum bound. **The initially identified pilot-conditioning defect is now resolved:** the parent revised the interpretation and this worker reread the changed passage. The current text explicitly uses deterministic-design identities and excludes automatic pilot-conditional and unconditional uncertainty for the reused-data pipeline. There is no reason in these documents to cancel the registered development ablation.

## Confirmed identities and realification

1. For real r and E, differentiation gives z*=(I+E^T E)^-1 E^T r. Substitution gives r^T[I-E(I+E^T E)^-1E^T]r=r^T(I+EE^T)^-1r. Thus the claimed penalized-error interpretation is exact (interpretation lines 9–20). It creates no new measurement information.

2. Let d=u+iv, R(d)=[u;v] and R(id)=[-v;u]. Their real dot product is zero and their squared norms are both ||d||². With delta=sqrt(2)R(d)/sigma, the rank-two columns are R(d)/sigma and R(id)/sigma. Each has squared norm ||delta||²/2, so rank-two, rank-one and isotropic increments all have trace ||delta||². If d=0 they all reduce to identity; the phrase “two nonzero covariance increments” requires d nonzero. This is a harmless degenerate-case wording qualification, not an algorithmic contradiction.

3. In an optional stochastic interpretation, rank one corresponds to one real amplitude along d; rank two corresponds to a proper complex scalar amplitude along d with the same total expected discrepancy energy. Rank one need not represent proper complex noise. Working in realified coordinates makes either construction legitimate as a surrogate penalty; no proper-noise claim is required for the proxy.

4. With proper complex measurement variance E|epsilon|²=sigma², sqrt(2)/sigma realification gives identity covariance. For fixed C positive definite, symmetric W=C^-1/2 satisfies WCW^T=I. Applying W to both residual and Jacobian is the correct fixed-weight least-squares derivative. The intended implementation tests explicitly cover these points (protocol lines 49–51); this review has not executed them.

5. For deterministic positive-definite P, full-column-rank J and zero-mean noise with covariance I, A_P=(J^T P J)^-1J^T P obeys A_P J=I. Error is A_P b+A_P n, its noise covariance is A_P A_P^T, and for fixed b=F eta with ||eta||<=1 the worst expected target squared error is tr(T A_P A_P^T T^T)+||T A_P F||_2². The inverse weighted normal matrix is a covariance only if the actual residual covariance is P^-1 under the stipulated fixed linear model. Nonlinear fitted-model Hessians require additional approximation/regularity qualifications.

## Initially required clarification, now resolved: conditioning on the fitted pilot is not enough

Interpretation lines 49–52 correctly acknowledge data-dependent P but warn specifically about **unconditional** uncertainty. The preceding phrase “actual conditional covariance” (line 39) can still be misread as a valid covariance conditional on the pipeline's frozen pilot.

**Recheck of revised text:** current lines 31–42 explicitly define deterministic design, zero-mean unit-covariance noise and deterministic P. Current lines 50–57 explain that conditioning on reused low-band data changes the noise mean/covariance and grants neither automatic pilot-conditional nor unconditional uncertainty. They also correctly identify independent pilot measurements as extra data and a changed protocol, with nonlinear uncertainty still requiring care. This resolves the reporting/math-assumption defect. The paragraph above and derivation below document the reason for the correction, not an outstanding objection.

The low pilot is estimated from observations that are reused in the final fit. Consequently, conditioning on that pilot or P generally changes the low-block noise mean and covariance. Even if the high-block measurement noise is independent of the low pilot, the entire stacked noise no longer automatically has conditional mean zero and covariance I. Freezing P during optimization suppresses parameter derivatives of P; it does not remove statistical dependence between P and the reused data.

**Minimal wording repair:** state that the displayed covariance and risk are fixed-design identities for P deterministic or independent of the measurement noise being fitted. State explicitly that neither pilot-conditional nor unconditional uncertainty for the current reused-data pipeline follows automatically.

If J and b are also treated as fixed after conditioning, the appropriate conditional noise covariance is A_P Cov(n|P) A_P^T, and the conditional mean error includes A_P E[n|P] in addition to A_P b. If J/b themselves depend on the pilot, they require corresponding conditioning as well. No new uncertainty computation is needed for the registered experiment because it already prioritizes empirical recovery errors and makes no calibrated-interval claim.

## Confirmed scope and remaining boundaries

- Coarse/fine prediction differences computed at the low pilot and including fitted electronics are empirical discrepancy directions. Their outer products do not establish the magnitude, orientation, mean or covariance of the unknown continuum discrepancy. The protocol states this correctly.
- Equal trace fixes total surrogate inflation, not retained target information or noise precision; isotropic and directional controls can legitimately behave differently. Legitimate pose/material sensitivity may be suppressed together with mismatch.
- No ADDA truth or high-band measurement residual enters mode construction. The low observations do influence it through the pilot, which is declared and acceptable for a development estimator. This is not independent validation or a noise-independent weight.
- The error-set matrix F in the risk identity needs separate justification. Substituting the proxy E for F would not turn the algebra into a certified worst-case bound; the document correctly requires a separately justified error set.
- Because the mode and strength are frozen within each fit, the likelihood-normalization/log-determinant term, if one chooses a Gaussian surrogate interpretation, is parameter-constant. Different methods nevertheless have different C, so their minimized weighted losses cannot directly be ranked as a common physical fit metric. The protocol's separate physical high-block and held-out evaluations are the appropriate comparison route.
- Distinct noisy-reference conditions add data and should remain separate. Shared pilot costs and historical timings are transparently accounted for. Two reused development scenes do not supply final confidence or novelty evidence.

## Parent handoff

The pilot-conditioning caveat has been repaired and rechecked. Optionally qualify the zero-discrepancy rank statement. The current documents define a coherent bounded fixed-weight ablation. Implementation accuracy, empirical benefit, closest-prior overlap and final acceptance are outside this mathematical document review.
