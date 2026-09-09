# Algorithm A4: implemented modules and admission boundaries

## 1. Scope and variables

Two implemented pipelines are deliberately distinguished. They are not advertised as one validated, autonomous six-action controller.

**Geometry/coverage:** r in R^3, in a declared union of nominal boxes; per-block unknown complex amplitude c_b is profiled exactly. A block is nine complex components for three regular electric-l=1 illuminations and three calibrated receiver components. k and known receiver dither d_b are supplied. Material, common delay and shared electronic amplitude are contained in c_b. This enlarges nuisance freedom and makes geometry independent of material coefficients in the restricted radial class. It does not estimate material.

**Joint quantitative calibration:** theta=(r_x,r_y,r_z,epsilon_real,g_real,g_imag,ell), seven real parameters. Imaginary permittivity is fixed at0.05, sphere radius0.045m, known support/center, k=(4,8,22)rad/m, shared complex gain and common delay length. Three noisy complex electronics references are additional data. Material and geometry are measured jointly in the declared scaled task. There is no unknown-support or arbitrary volumetric material inference.

## 2. Whitening and profiling

For proper complex entry noise of variance sigma^2, use sqrt(2)[Re z;Im z]/sigma for real inference. Never project raw, nonwhitened nuisance columns. An SVD of scaled nuisance columns permits rank-deficient/duplicate material and electronics directions. For geometry B_vis=(I-P_[A,G])B; full joint least squares is equivalent to profiling all other identifiable parameters for a task block. A scalar amplitude is optimized by c=(f*y)/(f*f), with complex conjugation in the inner products.

## 3. Fixed-estimator risk and implemented surrogate

For fixed J,W,T and true unit noise, L=(J^T WJ)^-1 J^TW. The covariance contribution is ||TL||_F^2. For deterministic discrepancy d add ||TLd||^2; for a verified ellipsoid d=E eta,||eta||<=1 the worst squared bias is ||TLE||_2^2. These are different from the stochastic sample-model risk used in the implementation:

R_hat=||TL||_F^2+||TL m_e||^2+||TL E_s||_F^2,

where E_s E_s^T is a sampled covariance. Its Frobenius term is not the ellipsoid operator norm. The result is a surrogate, not a confidence certificate. The independent low-band pilot and independent pilot references avoid reusing the same noise to construct and fit the weight. Nonlinear linearization and representativeness of the prior error family still remain assumptions.

The task map is T=diag(1/0.015,1/0.015,1/0.015,1/0.15,0,0,0), with its zero rows omitted. Report geometry/material errors separately as well.

## 4. Joint estimation algorithm actually run

1. Fit the common coarse-model low-band pilot from independent observations and independent electronics references.
2. Draw32 prior parameters from declared nominal/bounded ranges, not from truth. Evaluate fine-minus-coarse fields, including electronics. Estimate mean/covariance; fit a regularized linear conditional mean for the conditional-error baseline.
3. Freeze transforms. Isotropic and rank-one choices suppress only the high-frequency block, with matched covariance trace. Sampled and conditional-error models also correct the mean; the latter is a linear regression surrogate, not an exact conditional posterior.
4. At the pilot, compute coarse and fine task Jacobians and sandwich risks for raw, isotropic, rank-one, sampled error and full fine. Charge every screen/model evaluation. The risk selector chooses the smallest risk; this version has NO explicit cost penalty.
5. Solve every baseline from the same pilot and the same matched noisy data/reference. Coarse-to-fine pays for both stages. Save all endpoints before assessment by parameter truth.
6. Report final task error and total charged time, not the surrogate objective value as success.

The selector chose full fine in all8 frozen scenes. This version is retired as a main-contribution candidate rather than relabelled successful because its chosen estimates are accurate.

## 5. Coverage with an explicit unresolved state

Let f_b(r)=vec M_k(r+d_b)/||M_k(r+d_b)|| and P_b=f_b f_b*. The profiled residual is ell(r)^2=sum_b||(I-P_b(r))y_b||^2. With m complex entries and a declared discrepancy norm beta, use tau=sqrt(chi2_(2m,1-alpha)/2)+beta in complex-whitened units. The code uses alpha=.005 in training and .005 in validation, a conservative .01 two-gate budget; the training outer-set theorem itself only needs its own noise event.

For cell center c, half widths h and rho=||h||, compute the exact lower radius of its translated box to each origin. The analytic bound ||D P_b||<=sqrt(3)/R yields a residual lower bound

LB=max(0,ell(c)-sqrt(sum_b[||y_b|| min(1,sqrt(3)rho/R_min,b)]^2)).

The implementation subtracts a small floating-point pad. This is not directed-rounding interval arithmetic. Tests compare bounds to sampled cell points; they do not replace a machine-verified proof.

```
queue <- every root cell covering the declared prior domain
inside <- empty
while queue nonempty and examined < budget:
    C <- cell whose farthest point is farthest from current estimate
    if LB(C) > tau: exclude C
    else if every point of C is within task tolerance of estimate:
        retain C in inside
    else if C cannot be refined under resolution limit:
        retain C as unresolved
    else: split its longest edge and retain both children
retain every unprocessed cell as unresolved
accept_conditionally only if estimate is feasible, inside is nonempty,
    unresolved is empty, and independent validation is consistent
otherwise reject or acquire; never silently discard unresolved cells
```

Frozen coverage uses6000 examined cells and15mm geometry tolerance. Exact scalar profiling avoids the invalid use of a local material optimizer's residual as a cell lower bound. Extending this bound to general full-wave material profiles is unimplemented.

## 6. Acquisition and validation

A noiseless modal tensor gives two candidates +/-r from its simple normalized eigenvalue. This algebraic candidate generator is a strong prior-based baseline, not an oracle. A nonzero known receiver dither removes this exact ambiguity within the model. Candidate dither pool: +/-35mm along each Cartesian axis. The branch score maximizes the minimum squared angle between frozen candidate patterns at k=12. The Fisher control maximizes the ideal local minimum singular value at the leading candidate; random selects from the same pool.

All acquisition policies add one nine-complex-value tensor and an independent nine-value validation tensor. Original training/validation each contain27 values at three frequencies. Frozen potential measurements/noises are shared across policies. Absolute noise at the new point uses the baseline k=12 noise scale; new-point SNR is not separately renormalized. Physical acquisition and validation counts are separate from CPU time. Scores assume informative nonzero responses; no guaranteed unknown-SNR optimality is claimed.

Selection/fitting uses training only. Independent validation checks a fixed geometry while reprofiling its scalar nuisance; it is a feasibility gate, not the old fixed-mean finite-bank Gaussian theorem.

## 7. Six-action interface: implemented logic, not fully validated evidence fusion

| Action | Required evidence / trigger | Actual status |
|---|---|---|
| accept | Task risk below declared tolerance, branch coverage, supported fidelity, consistent validation | Conditional geometry cover implemented; joint physical accept not certified |
| reject | Missing coverage/evidence, bad model/prior consistency, unresolved budget, or no justified remedy | Executed in branch tests; can fail if error family is falsely declared valid |
| downweight | Candidate sample-model task risk below both raw risk and tolerance after nuisance profiling | Several weighted baselines executed; six-action interface unit-tested only |
| refine model | Supported fidelity failure or lower predicted task risk than suppression | Joint selector executed, always fine in frozen scenes, no cost win |
| add electronics reference | Remaining material/gain ambiguity and an available independent reference that constrains it | Separate known-geometry material experiment; not unified online diagnosis |
| acquire new data | Multiple admissible geometry regions or inadequate covered task precision | Fixed-pool dither policy executed and charged |

`risk.action(Evidence)` assumes that the evidence fields have already been established. It does not infer them from residuals, and should not be mistaken for an implemented general physics classifier. An unknown discrepancy set must not be encoded as `fidelity_supported=True`.

## 8. Complexity and accounting

Ideal tensor evaluation is O(9) per block. Profiling uses inner products; C geometry cells cost O(C*9*number_of_blocks), apart from optimizations and bookkeeping. The in-house DDA forms a3N-by-3N dense complex system: O(N^2) storage and O(N^3) factorization, three incident right-hand sides plus a material derivative right-hand side. It is a verification engine, not a scalable tomography algorithm.

Joint screening needs finite-difference Jacobians of both fidelities;32 fine/coarse error pairs are charged per scene, not amortized after inspecting favorable outcomes. Counts and wall times are retained. Fine is an analytic Mie coefficient here, not an expensive high-resolution Maxwell solve, so this benchmark cannot establish cost superiority for large numerical models.

No monetary/hardware acquisition time is invented. Illumination-synthesis channels, antenna-switching time, vector-field receiver calibration, motion-stage accuracy and solver certification costs are unpriced and remain part of the hardware feasibility gap.
