# Independent development experiments

## Completed native verification integrated with recovery evidence

The automatic run has ended; it did not generate another recovery campaign. Its 53-entry ledger labels 41 matches (including input checks), two weaker bounds, seven unmatched constants and three incomparable records. Parent verification matched 31 local artifact hashes. The approximately 44.34% scalar risk is independently corroborated.

The repaired fixed-loss Born pair is (1.9,2.5) versus (2.2,3.0), with partner gain modulus 0.750182155685. The old q'=45, g=1 pair violates the annulus. Scalar reference GLS agrees numerically; annulus constraints and complete residual derivatives matter. Parent formula corrections and interval limitations appear in Q5_THEORY_FINAL.md. No new 128-scene test, hardware result or continuum certificate was produced. Counts below are unchanged.

## Development recovery


Two known spheres at (-0.06,0,0) and (0.055,0.02,0) m, radii 0.035 and 0.025 m; k=18; losses (0.03,0.05); four illuminations/polarizations; 12 three-component receivers at radius 0.6 m. Material domain [1.5,4]×[2,5], receiver x translation ±2 mm, gain modulus [0.75,1.25]. Complex noise standard deviation is 1% of field RMS. Complex reference standard deviation is 0.01.

Treams lmax=3 is inverted. Same-principle sanity data use lmax=4. Independent data use DipoleVIE spacing 0.011 m and fill quadrature 4. The latter is neither continuum truth nor a certified discretization bound. Independent solver convergence remains open.

| Development condition | Joint task success | Median material absolute errors |
|---|---:|---:|
| Same model, no reference | 11/12 | 0.0229, 0.0271 |
| Same model, reference | 12/12 | 0.0083, 0.0110 |
| Independent DDA, no reference | 1/12 | 0.3794, 0.4780 |
| Independent DDA, reference | 6/12 | 0.0929, 0.1051 |
| Independent, selected EM scalar | 1/12 | 0.3873, 0.4593 |
| Independent, random EM scalar | 1/12 | 0.3791, 0.4795 |
| Independent, fixed EM scalar | 1/12 | 0.3814, 0.4695 |

The two-action policy succeeds on 3/12, choosing reference in five scenes. All second-cycle scenes reuse development data. No confidence claim on an unseen population is made. First and second cycle elapsed times were 94.24 and 87.67 seconds; these are complete cycle costs, not matched per-method runtime superiority evidence. Memory was not instrumented.

The corrected independent no-reference median geometry error is 0.0440 mm despite large material errors. Its fixed-world structural error median is 14.99%; reference reduces that median to 3.90%. These illustrate wrong attribution and conditional reference benefit, not universal recovery.

Original cycle-one chi-square flags omitted the factor 2 required by proper complex noise, and structural fields were evaluated at shifted rather than fixed physical points. Both metrics are superseded by results/development_audit.json. Even the corrected residual test is heuristic: active constraints, nonlinear fitting and model mismatch invalidate a generic calibrated acceptance interpretation.

Formal 12-scene tuning and 128-scene final evaluation were not generated. Coverage certificates, independent convergence, electronics-breakdown cases, complete mechanism ablations and statistical gates remain outstanding. Stopping the current unsuccessful method after two development cycles follows the registered route-switch rule.
