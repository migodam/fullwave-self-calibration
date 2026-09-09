# Bounded Flash worker: A2 numerical self-calibration solver and E4 runner

Root: `/Volumes/migodam's-external-brain/Research/Inv_SLAM`. All relative paths
refer to this root. Do not spawn workers or AI Scientist, change configurations,
install infrastructure, write secrets/reasoning, or edit historical experiments.
Use apply_patch. Write code/results/summary only in `research/delegated/a2_solver`.
Read root AGENTS.md, A2 validation protocol, and the files below yourself.

Physics implementation is being finalized in `research/delegated/a2_physics`.
Use that API; do not modify it. The parent controls scientific conclusions.
Read the parent `research/trispace_self_calibration/a2_research/THEORY_AUDIT.md`
and import `certificates.py` from that folder for the passive-medium certificate.

## Task and execution boundary

Implement the E4 solver, cost ledger, finite tuning and final-test runner. RUN
ONLY tuning seeds1-10 and unit/smoke checks now. Do NOT generate or inspect final
seeds1001-1020 until the parent explicitly runs your `--mode final` command.
No final test access from tuning. Do not alter the physical specification.
Use existing experiment venv, CPU with BLAS threads1, no installation.

## Models and methods

All physical methods solve the SAME total-complex Gaussian objective with nine
real shared material coefficients and three shared geometry parameters, r_free=0.
Use coefficient boxes [0.03,2.5] and pose box translations[-.7,.7]m, rotation
[-.7,.7]rad. Fixed optional quadratic prior must be offered identically to all
methods; default zero (box constraints still explicit). Scale pose by lambda_min
and lever arm1.5; material scale1. All data are proper complex, total-field SNR30dB.

Methods:
1. direct: exact physical solver with analytic adjoint gradient and L-BFGS-B.
2. prasc: certified SOM-informed reduced physical states/sensitivities, same
   L-BFGS-B and total physical objective. Fall back to exact solves if certificate
   is unavailable or fails. No freely data-fitted current nuisance.
3. fixed_rank: same reduced numerical family with a fixed rank selected in tuning,
   no adaptive certificate (record its actual physical error for evaluation).
4. phaseless: same geometry/material boxes and matched parent total intensity.
   Exact induced Rice/noncentral-chi-square negative loglikelihood, stable i0e/i1e;
   gradient chain to total field. No additive Gaussian intensity surrogate.
5. direct_control: exact solver with the same diagnostic frequency/control policy
   as prasc; separates generic continuation/control from reduced-state effects.

Core.forward supplies full A/B at high RHS cost; use its analytic adjoint
gradient when available for direct/phaseless. If absent implement an adjoint
wrapper against returned state matrices. Never hide parameter derivative costs.

SOM-informed numerical basis: at current geometry stack the per-pose sensing
S (avoid duplicated illumination rows), take leading right singular vectors.
Optionally append sequential TSOM domain vectors P_(S-)V_(D,+) by RRQR, with
actual rank recorded; do not call it an exact intersection. At a fixed linear
algebra solve, all state coefficients minimize ||M U c-b||, not measurement
fit. Allow data-subspace residual corrections; do not freeze noisy j_det.
Rank ladder [8,16,24,36,64,96,128,192,256], capped at n, monotonically promoted
within a stage. For rank n use direct solve. Keep L_det/r_num/r_free distinct.

Important derivative design from parent: you may approximate the PHYSICAL tangent
equations M t_v=b_v-M_v*j_tilde in the same basis and use the parent residual
certificate. This is an inexact physical Jacobian, NOT a claimed exact derivative
of the moving-basis surrogate. Thus no unaccounted DU term is asserted zero.
Freeze the numerical basis during a line search; rebuild at accepted chart
events, reset quasi-Newton curvature when the chart changes. Include all S_v
terms. If you instead differentiate a moving surrogate, include DU/Dj_det.

Evaluate passive_bound from estimated positive material and known tau only.
Bound state data error and total derivative error. Label any uncertified
heuristic as a diagnostic. A practical initial candidate acceptance is residual
error <=0.05*sqrt(number real measurements) and relative physical gradient
error <=0.1 (use certificates.gradient_error_bound); tighten near stationarity.
If such conservative bounds force direct fallback, REPORT THAT, do not silently
replace them by true-error checks or estimates. Objective acceptance may use
certified intervals, otherwise exact physical checks (charge them).

For frequency policy use the same low-to-high cumulative frequency data schedule
for all methods; distribute common200 work cap across cumulative stages with
cumulative budget fractions [.12,.28,.50,1]. Secondary800-cap schedule identical
fractions. A covariance gate is ONLY a surrogate, not a coverage certificate.
Compute/use additional observability gates only if their cost can be explicitly
charged. Record unavailable diagnostics rather than accessing true parameters.

## Budget and failure

Count forward/adjoint RHS exactly, factorizations and SVD separately, operator
products and reduced solves. Calibrate reduced work by a repeatable N16 timing
microbenchmark on tuning only; save raw timing samples and conservative factors.
Wall time includes basis/certificate/acceptance costs. Do not charge a multi-RHS
solve once per matrix. Catch budget exhaustion before the next full evaluation,
preserve last finite ACCEPTED iterate; report incomplete stages and budget failure.
Keep evaluation-only exact Jacobians/refinement/oracles separate from estimator
time and never feed them into estimator decisions.

## Preregistration generation and tuning

Two coefficient truth ranges: weak uniform[.15,.5], stronger uniform[.7,1.4].
Two apertures full/limited as physics spec. Assign strata cyclically to seeds,
so final20 seeds have5 each. Geometry truth shared error zero; supplied nominal
initialization offset is the unknown correction to recover. Do not use truth in
algorithm. Initial material is the same bounded low-frequency material-only fit
at the nominal supplied pose for all methods, with that initialization's work
charged identically. If init alone would exceed a budget, report and use a
predeclared fixed nonzero nominal alpha=.5 for ALL methods, documenting this
design change before final generation; never initialize zero scatterer.

At each test seed initialize x at12 error combinations from A2, with radii
lambda*[1/8,1/2,1],4 angles and half squared lever-arm metric in translation and
rotation, alternating rotation sign by seed. Noise on inverse observations and
an independent validation replicate generated with N32. Noise sigma must derive
from a declared nominal reference total field, not be secretly tuned per fit.

Before final data generation, calculate reference design covariances at nominal
alpha=.5,x=0 for each aperture, material task Q_res with width
max(2*grid_h,lambda_min/(2*NA)), NA1 full/sin(pi/4) limited. Use A2 tolerance and
noninferiority definitions and publish actual numbers. Do not use regularization
or pseudo-inverse null zeros to declare an unsupported task full-rank. Save
rank/condition diagnostics and flag reference design failure for parent review.

Finite tuning grid for optimizer: maxls[10,20], ftol[1e-8,1e-10], gtol1e-6;
fixed ranks[16,36,64]. Compare on1-10 only with equal candidate count per physical
method, pick joint error/objective behavior without test knowledge. No endless
parameter search. Parent may simplify if CPU constraints require, but record it.

## Artifacts/API

Provide solver.py with callable solve(model,y,sigma,alpha0,x0,method,budget,settings)
and run_e4.py --mode tune|final|summarize --budget 200|800 [--methods ...].
Final mode MUST require an explicit --frozen-config path, no implicit tuning.
Raw one-record-per-run JSONL/CSV includes seed,stratum,init,method,config hash,
all metrics, estimates, failure, covariance/coverage diagnostic flags and counts.
Evaluator uses A2 joint success, fixed resolved map RMSE, pose lever-arm error,
heldout residual (with declared discretization allowance, not silently assumed
noise-only), and no false calibrated labels. Return statuses as evidence, not
SciPy success codes. Write tests, `tuning_summary.md`, timing and a compact summary.
