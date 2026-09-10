# A4 + A4_2: material observability is not the geometry window

Parent mathematical integration, 10 September 2026. These deductions are not
a novelty claim. Original A4/A4_2 text and A3 frozen evidence remain unchanged.

## Corrections that determine the research direction

1. **Task reversal.** A4 Theorem1 concerns geometry after material/gain profiling.
   Its Proposition3 says the allowed material is unidentifiable. Thus its kR=1
   fixed-relative-SNR optimum is not an optimum for reconstructing interior or
   exterior permittivity. A4_2 must not inherit that interpretation.
2. **No Fisher-versus-SVD escape.** For a real-whitened linear tangent, squared
   singular values of Q_B A are eigenvalues of its profiled information matrix.
   Switching notation does not remove local, model-mismatch or branch limits.
   Nonlinear finite-pair distinguishability is a genuinely different test.
3. **Not a universal intermediate optimum.** Receiver standoff does not alter
   the internal resolvent under fixed illumination and nonperturbing observation.
   Multiple scattering alone does not imply nonmonotonic distance dependence.
   Fixed angular aperture, fixed metric aperture and fixed power/noise are
   different design problems. Normalizing fields without transforming noise
   changes the comparison.
4. **Gauge admissibility.** Born scale symmetry holds only for c such that cX
   remains in the allowed material class and the compensating gain is admissible.
   Arbitrary complex scaling can violate passivity or fixed loss constraints.
5. **Full-wave scale breaking is not proved by a nonzero quadratic term.** The
   changed field must escape the gain orbit *after measurement*. A4 itself is
   an exact full-wave counterexample. See below.
6. **Tensor components are not independent targets.** To assess one unknown
   component A_q while other material components are also unknown, profile
   [B,A_-q], not merely B. Otherwise componentwise windows can overstate joint
   identifiability. Unequal polarization sensitivities do not prove disjoint or
   split standoff windows; nullspaces may persist at every distance.
7. **Far-field remainder.** A4_2's scalar radial expansion omits explicit k a^2
   dependence. At fixed k,a, the directional term yields O(k a^2/d^2), while
   amplitude contributes O(a/d^2). Do not claim a high-frequency-uniform bound
   from O(a/d^2) shorthand. The prior A3 conditional theorem already exists.
8. **Global ambiguity is not estimated from a small candidate grid.** A finite
   ambiguous pair gives a lower bound on the worst-case modulus, not an upper
   stability certificate. If the gain family permits zero, all materials can
   yield zero data; even the proposed global modulus can be trivially maximal.
   Require a signal/gain lower bound or informative reference when appropriate.

## Proposition A: exact projective obstruction

Let measured nonzero blocks satisfy y_b=g_b f_b(alpha), with unrestricted
nonzero complex gains independent across blocks. Two admissible materials
alpha and beta are indistinguishable exactly when each f_b(alpha) is complex-
collinear with f_b(beta). If f_b(alpha)=t_b(alpha)v_b for fixed v_b and nonzero
t_b, every material pair is indistinguishable: replace
g_b by g_b t_b(alpha)/t_b(beta). This is an exact nonlinear argument, not FIM.

A4's radial electric-l=1 acquisition has exactly this factorization. Changing
range or adding independently gained frequencies preserves it. Geometry can
remain observable through v_b(r), while material is not. This restates and
connects A4 Proposition3; it is not newly discovered prior art.

With bounded gains, the same counterexample applies only if the compensating
gains remain admissible. With a shared gain across multiple blocks, collinearity
must hold with one common multiplier, which can be a stronger condition.

## Proposition B: an observable second-order condition

For one shared gain, let f(s)=s v1+s^2 v2+O(s^3), v1 nonzero, with the remainder
differentiable and s>0. Let P_f be complex orthogonal projection onto f. Then

$$
(I-P_f)\,s f'(s)=s^2(I-P_{v_1})v_2+O(s^3).
$$

Proof: s f'(s)-f(s)=s^2 v2+O(s^3); projection kills f(s), and P_f=P_v1+O(s).
Consequently the relative visible log-scale sensitivity satisfies

$$
\frac{\|(I-P_f)s f'(s)\|}{\|f(s)\|}
=s\frac{\|(I-P_{v_1})v_2\|}{\|v_1\|}+O(s^2).
$$

Therefore a measured noncollinear second-order contribution is sufficient
for a nonzero local scale direction at sufficiently small positive s. Merely
v2 nonzero is insufficient. If v2 is collinear, higher orders may or may not
help. Multiple independently gained blocks require the blockwise projection.
This is an analytic local condition, not a global uniqueness theorem or a
guarantee of recovery above a noise/model floor. The resolvent Neumann expansion
requires convergence near s=0; strong-scattering conclusions cannot be obtained
by truncating it outside that region.

## Current executed evidence and limits

Nine A4 sphere/material/range/frequency counterexamples compensate material by
gain to relative error at most1.55e-16. Separately,18 two-sphere Treams cases
produce nonzero profiled scale changes under shared or per-illumination gains;
independent entry gains erase the change. Complete run time7.92s, deterministic
mechanism experiment, not an inversion benchmark.

At R=.6m, increasing contrast scale from.01 to2 raises relative shared-gain
log-scale sensitivity from.001243 to.07121. A25% material-scale change leaves
profiled relative residual.000310 to.01613. This warns that positive tangent
rank can coexist with extremely small finite separation. No specific hardware
noise level was imposed, and no recovery probability follows.

Largest lmax3→4 field change is1.09e-5; derivatives were also compared. This is
one truncation sensitivity check, not an independent solver/continuum proof.
The experiment does not isolate inter-sphere multiple scattering from internal
sphere response: a no-interaction control is required before causal attribution.
The linear scaling control is an algebraic homogeneous reference, not a physical
Born calculation. These limitations determine the next experiment, not an
excuse to claim the positive result early.

## Next substantive experiment

Hold total measurements and gain-sharing fixed; compare one versus multiple
calibrated modal orders, shared versus independent frequency gains, and an
independent electronic reference. Use two separately unknown radial material
regions and known support first. Profile the *other* material region and
geometry when evaluating each target. Measure actual material recovery and
finite competing solutions, not just singular values. Add no-interaction and
physical Born controls before claiming nonlinear scattering is the source of
the benefit. Hardware modal synthesis and anisotropy remain explicit later gates.
