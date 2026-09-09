# A4 v1 frozen protocol

Frozen after development v2, before inspecting any seed below. This is a local
preregistration (timestamp and source hashes), not an externally registered study.
Primary questions are restricted and may yield negative results. No A3 endpoints
are run, edited or reclassified.

## P1: analytic electromagnetic condition

41 dimensionless radii, seed 91031; exact vector Maxwell electric-l=1 model;
three regular modal illuminations and three calibrated Cartesian components.
Compare analytic profiled singular values against direct realified projection.
Report all results, low/high radius slopes, and degeneracies.

Independent forward check: four support classes (sphere, near-sphere, ellipsoid,
box), two contrasts, three wavenumbers, two grids, two receiver radii: 48 DDA
state solves and 96 receiver cases. Exact Mie comparison only for spheres.
DDA extensions are numerical evidence, not an error enclosure or a theorem.
Diversity: two orientations x three component counts x three illumination counts
x three gain models = 54 algebraic physical-tensor controls.

## P2: coverage and acquisition

Seed 92031. Twenty-four independently sampled scenes, four for each of six
strata: easy anchored hemisphere, antipodal ambiguity, intentionally missing
initial bank, low SNR, nonspherical model discrepancy, and an exact nonlinear
two-world discrepancy control. One noise draw per scene, not repeated noise
reported as independent scenes. Independent validation noise and matched
potential acquisition noises are shared across methods.

Six policies: naive one-root multistart, algebraic two-root multistart,
profiled local-Fisher acquisition, random acquisition, branch-separation
acquisition, and branch acquisition plus budgeted coverage.

The initial domain consists of one or two explicit 3D boxes. Its inclusion of
truth is a physical-prior assumption, not certified by the data. Branch cells
retain all pending regions at budget exhaustion. Maximum 6000 visited cells.
Task error tolerance 15 mm. Training exclusion alpha=0.005, validation alpha=0.005.
The conditional beta=0 model applies only to the exact-modal strata; report
out-of-model strata separately. Double precision is not interval certification.

Primary endpoints: wrong acceptance, coverage missing, wrong selection conditional
on a qualifying represented branch, correct estimates rejected. Also report
geometry error, accepted fraction, unresolved regions, all fitting/design/
coverage/validation time and measurement counts. No claim of broad superiority
from 24 heterogeneous scenes. This is a symmetry-branch experiment, not a test
of all full-wave cycle-skipping phenomena.

## P3: attribution baseline stress

Seed 93031. Eight new scenes, nine methods = 72 fits. Joint 3D receiver position,
one real permittivity, shared complex gain and common delay. Known sphere radius
45 mm, fixed loss part 0.05, k=[4,8,22]/m; independent pilot data and reference
noise. Fine simulator: analytic full-wave electric Mie response. Coarse simulator:
Rayleigh scattering coefficient but the same exact exterior Maxwell dyad.

Methods: low-only; raw all-frequency; isotropic high-frequency downweight;
rank-one high-frequency error downweight; sampled mean/covariance approximation
error; conditional Gaussian-regression approximation error; coarse-to-fine;
full fine; a fixed-pilot surrogate-risk controller choosing raw/isotropic/rank-one/
sampled-EEM/full-fine. EEM is a standard borrowed baseline, not a new method.

32 independent prior training pairs per scene, charged without amortization.
Identical independent pilot and initialization for all methods. All noisy
references are charged: three fitting plus three pilot reference values.
No true parameter is supplied to the controller. Known sigma is part of the
simulated measurement model, not an estimated noise oracle.

Task loss = (Euclidean geometry error / 15 mm)^2 +
(permittivity absolute error / 0.15)^2. Report components separately. Paired
controller/full-fine loss and total cost are the primary contribution gate.
If loss is the same and cost greater, the controller is not a main contribution.

## P4: electronics intervention

Seed 94031, sixteen fresh material-only scenes with externally known geometry.
Three arbitrary per-frequency complex gains. Without reference, material loss
is exactly flat after gain profiling. Add three noisy complex electronics
references and refit material and gains jointly; never substitute true gains.
Report rank change, estimation errors and cost. This does not by itself close
joint geometry/material reliability or demonstrate superiority to prior methods.

## Statistical and cost rules

Keep every run, exception and bound failure. No retries selected by quality.
Do not tune on final scenes. Paired confidence summaries are descriptive;
finite-scene binomial bounds are not guarantees for a different population.
Reference/acquisition counts must be reported separately from CPU time.
Analytic Mie and in-house DDA are independent formulations; no external ADDA
or Treams execution is claimed. No hardware measurement is fabricated.
