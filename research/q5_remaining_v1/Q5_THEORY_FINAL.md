# Q5 theory audit — not a completed final theorem package

## Completed automatic verification: integrated parent verdict

The native process ended with exit 0, one development round, an 18-turn cap and evaluator decision `abandon`; no native writeup ran. A bounded verification is not a stand-alone paper, but its valid checks remain useful. The worker labels 53 entries as 41 matching, two weaker bounds, seven unmatched and three incomparable. Matching entries include declared inputs: these are not 41 independent proofs. Parent verification matched 31 artifact hashes.

### Admissible fixed-loss Born ambiguity

Write u=epsilon-1 and ell=(0.03,0.05). If u=q ell, the coefficients are g(q+i)ell. Choose q=30, q'=40, g=1 and g'=(30+i)/(40+i). Materials (1.9,2.5) and (2.2,3.0), with |g'|=0.750182155685, satisfy the declared domain and give identical Born data for any fixed B. The parent's independent check confirms this identity. The original q'=45, g=1 has |g'|≈0.666872 and violates the annulus; multiplying both gains by 1.2 repairs it. Fixed loss therefore anchors scale only away from the parallel ambiguity. This is not a universal full-wave or unknown-geometry counterexample.

### What the completed checks do and do not establish

Mie conversion and binary risk agree with the parent's independent calculations. Worker V1 amplitude derivative lower estimates (4.2929719577e-4, 1.2943650641e-4) are weaker than the supplied values by 0.768% and 1.092%; different interval enclosures do not refute the original bound. The tiny floating-point identity residual was not reproduced; its precise magnitude is not a physical invariant. Original random arrays and runtime are not independently reproducible from the available context.

The scalar noisy-reference GLS objective agrees at approximately 1.34e-16 relative. Adding a log-determinant changes that objective. Enforce the gain annulus and retain the residual-dependent derivative term. These are correctness checks, not a novel optimizer or advantage over equivalent GLS.

An abstract additional nuisance column reduced a worker visibility-bound ratio below one. Unless it is shown to arise from admissible Maxwell geometry, this is only an algebraic scope warning, not a physical geometry counterexample.

### Parent corrections overriding raw worker derivations

In worker results/DERIVATIONS.md, epsilon means relative permittivity, not contrast, and finite-precision arithmetic is not exact evaluation. For j=j_1 and D_j=j'+j/z, the correct identities are

$$D_j'=-j'/z-j+j/z^2,\qquad D_j''=-j''/z-j'+2j'/z^2-2j/z^3.$$

The raw document's alternative D_j'=-D_y+2D_j/z and its second-derivative coefficient -2j'/z are wrong. The worker's directly differentiated code uses a different, correct expression for D_j'' at the parent's checked points; textual mistakes alone do not invalidate that implementation.

A midpoint Lipschitz enclosure needs L times the full halfwidth, not an unexplained factor 0.5. Endpoint-infimum and strict continuum claims are not promoted from the raw report. Exact size parameters, outward rounding and enclosure coverage remain proof obligations. The integration audit is results/integration_audit.json. These corrections take precedence over all raw worker summaries.


The fixed-loss Born anchor, restricted lossless modal boundary, and simultaneous lossy two-sphere inverse problem are distinct classes. No theorem is transferred between them without its assumptions.

## Verified implementation identity

With outgoing Riccati-Hankel convention, let the electric Mie coefficient be a=N/(N+iD). Then T=-a and q=N/D=i a/(1-a)=-i T/(1+T). The independent parent check matched the corresponding Treams eigenvalue at ten material/size points to 8.74e-19 absolute error. Helicity-basis diagonal entries mix electric and magnetic responses and are not the appropriate comparison. This is a unit test, not a convergence theorem.

## Conditional implementation-error corollary

The parent separately reconstructed the supplied two-world scalar calculation using SciPy spherical Bessel functions and root finding, without importing worker formulas. Worlds (2,3), g=1 and (2.359740580371099,3.9893194150178264), g=0.8 give whitened distance 0.201236624037519 and equal-prior Gaussian testing error 0.443423189256257. Each material change exceeds twice 0.1, so the corresponding disjoint estimation targets inherit the binary obstruction within this restricted class. This confirms the numerical A5 example, not its novelty or physical modal-readout realization. Evidence: results/parent_boundary.json.

Let h_i=|t_i| have h_i'≥m_i>0 and h_i≤H_i on a declared material interval. Observe Y_i=g t_i+e_i+n_i and Z=|g|+e_a+n_a. Suppose |g|≥a_min, |e_i+n_i|≤B_i and |e_a+n_a|≤B_a<a_min. Divide |Y_i| by Z, clip to the image of h_i, and invert h_i. Reverse triangle inequality gives

$$|\widehat\epsilon_i-\epsilon_i|\le\frac{B_i+H_iB_a}{(a_{\min}-B_a)m_i}.$$

Writing the bounds as model/reference bias plus noise, tolerance δ is certified if

$$b_{\mathrm{model},i}+(H_i+\delta m_i)b_{\mathrm{ref,bias}}\le\delta(a_{\min}-b_{\mathrm{ref,noise}})m_i-b_{\mathrm{noise},i}-H_i b_{\mathrm{ref,noise}}.$$

Proof: ||Y_i|-|g|h_i|≤B_i; division introduces at most h_i B_a in the numerator and denominator at least a_min-B_a. Clipping cannot increase distance to a feasible h_i; inverse Lipschitz continuity supplies 1/m_i. The displayed budget is an algebraic rearrangement. Probability is only that of the explicitly bounded noise event.

This adds an implementation requirement to the inherited restricted-class theorem, not a novelty claim or a guarantee for unknown geometry. Negative budget defeats this certificate, not all estimators. Supplied interval constants need an independently audited enclosure before use. Grid and solver differences are diagnostics, not proven error radii.

## Outstanding decisive mathematics

A conditional numerical budget illustrates the tight second-region requirement. Using provisional conservatively reduced slopes (0.000429, 0.000129), amplitude caps (0.00272, 0.00131), and the inherited A5 noise-event bounds gives material bounds (0.04982, 0.09399). At target 0.1, zero additional reference bias leaves additive modal-model budgets only (1.6076e-5, 5.7848e-7). The latter is roughly 0.044% of the second amplitude cap. These are **conditional arithmetic**, not proven hardware feasibility: the worker's enclosure backend, exact size-parameter construction, exported rounding and noise-event assumptions still require full audit. The supplied and reconstructed interval lower constants differ; differing interval expressions may produce different valid lower bounds, so this difference alone is not a refutation.

Uniform finite material separation for the simultaneous two-region domain remains unresolved. Candidate optimization yields an upper bound on the infimum, never a covering lower bound. No scalar-reference theorem alone repairs arbitrary structural model discrepancy. Exact-real interval implementation, floating-point export rounding and modal readout assumptions remain part of the proof audit.
