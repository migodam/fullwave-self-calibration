# Implemented algorithm and withdrawn method claim

Completed native checks are incorporated into the specification: retain the full residual derivative, enforce the gain annulus, and propagate reference noise through the concentrated GLS objective without silently adding a log-determinant. Scalar equivalence is not a novel optimizer. See Q5_THEORY_FINAL.md for corrected formulas. The selector-superiority claim remains withdrawn.


The shared forward interface maps (epsilon_1, epsilon_2, receiver shift) to a complex receiver field. For each nonlinear candidate, the complex gain is profiled analytically and constrained to the declared modulus annulus. A noisy complex gain reference enters as an additional weighted residual, not an exact calibration. Three bounded least-squares initializations compete on the complete profiled residual.

The development action selector constructs finite competitors by moving each material parameter by ±0.25 where feasible, optimizing the other material, receiver displacement and gain. It scores 576 scalar field candidates by worst separation over those competitors, then compares with a reference score. Candidate score is not a uniform certificate. All scientific acceptance fields remain unresolved.

Selected, random and fixed single-complex-observation baselines use the same added scalar budget and shared electronics assumption. The selector did not improve independent recovery and is withdrawn as a proposed superior method. Correct reference-first GLS and joint profiling can be mathematically equivalent; no superiority over the equivalent formulation is claimed.

The present implementation is a development diagnostic, not the plan's certified diagnosis–intervention–acceptance algorithm. No lower-bound-based acceptance mechanism has been implemented. See src/recovery_cycle1.py and src/scalar_cycle2.py, with immutable protocol versions under docs/.
