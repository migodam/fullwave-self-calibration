# Frozen bounded confirmatory ROM-seed protocol

Frozen 2026-09-08 after `timed_rom_checks.json` passed, before generating any
scene in seeds 9201–9240. This closes a narrowly defined comparison, not every
SOM algorithm or the full TAP acceptance plan. Original reserved seeds
1001–1020 remain untouched. No optional enlargement or early significance stop.

## Population and scope

Forty independently seeded 2D scenes, grid16, nine fixed Gaussian material
coefficients, full aperture, k=(3,6,12), 30 dB proper complex noise relative to
the nominal-background total field. Each has two or three randomly selected
material coefficients increased by Uniform(0.35,0.90) above 0.08. Anchored true
shared pose is zero; nominal translation magnitude is Uniform(0.04,0.08) m and
rotation magnitude Uniform(0.01,0.05) rad with random direction/sign. Separate
random streams generate material, noise and initial pose. Numerical forward
and inverse discretization are the same intentionally: this tests solvers, not
independent-physics validity or general calibration robustness.

Methods: direct GN, direct adjoint L-BFGS-B, generic-task, sensing-task,
Twofold-task, Krylov-task. Same data, material prior, bounds, rank ceiling128,
seed ceiling64 where appropriate, chart guards and physical trial acceptance.
All may take up to 200 iterations; none has the previous RHS cap. Joint
frequency fitting is identical for all. No method-specific tuning.

## Decision budget and timing

Each method has an eight-second **decision deadline**, including model setup.
Only iterates accepted before that deadline count. Record consumed wall time,
atomic-solve overrun, and all available accepted iterates. A numerical primitive
cannot be interrupted safely mid-call, so this is not claimed to enforce
identical consumed wall time. If it finishes after the deadline its output is
not accepted. All its work and time are charged. Offline evaluation is separate.
Methods execute serially with one numerical thread, in cyclically rotated order
across scenes. No concurrent heavy jobs. Retain every failure and timeout.

## Frozen estimands and analysis

Joint success at the deadline: lever-weighted pose error <=0.05 m and absolute
material-coefficient RMSE <=0.05. The lever is 1.5 m, as in development. These
are a restricted numerical benchmark's tolerances, not universal hardware specs.

Primary spectral comparison: **Twofold-task versus Krylov-task**. Sensing and
other method comparisons are secondary/descriptive, not a post-hoc substitute
if the primary fails.

1. Success advantage is the paired scene-level mean success difference. Use a
   conservative one-sided 95% upper bound from the difference of a 97.5%
   Clopper–Pearson upper bound for Twofold-only wins and a 97.5% lower bound for
   Krylov-only wins. If this upper bound is below 0.10, exclude the prespecified
   ten-percentage-point benefit in this benchmark population. Otherwise call it
   unresolved, not falsified. Complementary lower bounds govern positive claims.
2. Speed is evaluated only on pairs where both succeed and their final physical
   objectives agree within relative 1e-4. The estimand is the conditional median
   total consumed-time ratio Twofold/Krylov, NOT an unconditional mean speedup.
   Report the eligible count, ratios and exact order-statistic median interval.
   Test the 20% saving hypothesis (median ratio <=0.8) with the exact sign/binomial
   tail at 5%. State conditional selection and nonpreemptive-overrun limitations.

The two conditions are reported separately. To exclude the union of the two
specified advantage routes, both component nulls must be rejected; no positive
family-wide statement follows. Do not describe an eligible-subset speed result
as one applying to failures or all geometries. Report all endpoints, pose and
material errors, costs and activation regardless of success. No result is
examined for tuning while the fixed forty-scene run is in progress.

Complete final-population generalization to q49, unknown gains, larger phase
errors, other material priors or classical free-current SOM remains outside
this bounded protocol. The experiment may close a narrow claim while those
broader scientific gates remain open.
