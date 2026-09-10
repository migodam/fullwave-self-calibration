# Development cycle 1 — registered before new recovery results

2026-09-10. Not final validation. Seed 2026091101,24 development scenes: first12
same-principle multipole data, last12 independent DDA data. Uniform real epsilon
domains [1.5,4] and [2,5], fixed loss(.03,.05); x shift uniform ±2mm; common gain
modulus uniform [.75,1.25], phase uniform [-pi,pi]. Known centers/radii as A5.
k18,12 receivers radius.6m,four inherited plane-wave illuminations. Structural
field:36 receivers radius.15m; held-out prediction:16 receivers radius.6m.

Data: Treams lmax4 or original DipoleVIE spacing.011m,fill_quadrature4.
Inversion: Treams lmax3. These are deliberate development settings, not endorsed
continuum reference. Convergence is a separate prerequisite before final tests.
Noise proper complex with sigma=.01*field RMS; one noisy complex gain reference
sigma_ref=.01. Paired data across methods. Methods no_reference and joint_reference.
Three deterministic starts (2,3,0mm),(1.6,2.2,-1mm),(3.8,4.8,1mm), max80 function
evaluations each, tolerances1e-8. Gain is analytically profiled onto its allowed
modulus annulus; all three starts retained. Realification Re-all then Im-all.

Report per-material errors, position and gain errors, held-out prediction and
gain-free structural-field error. Scientific success eps≤.1,x≤.5mm,gain≤.05.
Heuristic residual chi-square99% diagnostic uses2m-5 dof (including2 reference
real observations where applicable); not coverage/acceptance certificate. Every
scientific acceptance remains unresolved until a justified error envelope and
finite-domain separation gate exists. Zero coverage may not be called safe success.

Record all failures, source hashes, model costs and solve time. No field renormalization
without noise adjustment. Do not tune on these outcomes and call them fresh test.
Final128 and tuning12 remain ungenerated. No acquisition advantage claimed here.
