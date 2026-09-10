# A4 closure supplement: causal interaction control

Registered before execution on 2026-09-10. Mechanism development, not a frozen
recovery benchmark. Preserve round-one sources and results unchanged.

Use the round-one two-sphere geometry, k=18, R=.6, 12 receiver positions and
four plane-wave/polarization illuminations. Six contrast scales .01,.03,.1,.3,1,2.
For each, compare the interacting lmax=4 Treams cluster with the coherent sum
of two isolated exact sphere responses at their actual centers. The latter
retains each sphere's internal full-wave response but removes inter-sphere
rescattering. Evaluate log-scale central differences at steps 1e-4 and 5e-5.
Use one shared unrestricted nonzero complex gain. Record absolute and relative
projected derivatives, finite 25% competing-scale separation, and field change.

Add a physical first-Born volume integral using k^2 free-space dyadic Green,
incident plane waves, and piecewise uniform sphere contrasts. Use midpoint
voxel quadrature at spacings .007 and .005 m, four subcells per coordinate for
boundary filling. This is a numerical Born reference, not a full-wave solver.
Report quadrature difference and weak-contrast Treams/Born mismatch separately.
No global, recovery, noise-threshold or continuum-certification claim.

Checks: Born scale is exactly gain-hidden; finite differences agree within
relative 1e-3 (record failure instead of changing threshold); all six cases
remain in results. Do not require interacting sensitivity to exceed isolated
sensitivity: this is the tested question, not a test harness assumption.

One CPU process, no model-provider calls and no new dependencies. Cost includes
all new field evaluations, Born quadrature, and diagnostic construction, but
not prior round-one runtime. Record hashes of protocol, source and dependencies.
