# PRASC-SOM: Preregistered Theory-Validation Protocol

Version: 5 September 2026. **Design only. No experiments have been executed for this document.**

This protocol tests the rank/bias/acquisition theory in the accompanying theorem package. Passing it is not equivalent to publication acceptance or to a universal basin theorem. The primary study is deliberately small: low-dimensional material parameters, scalar full-wave propagation, and a fixed budget. It is not a 3D Maxwell or general high-dimensional imaging validation.

## Common conventions

1. Preserve raw complex coherent observations. For a matched phaseless comparison, derive total-field magnitude/intensity from those same noisy total-field observations. Do not subtract incident intensity from total intensity and label the result scattered intensity.
2. Fix the time convention, material dispersion, proper/improper noise convention, units, current metric, and joint gauge anchors before generating data.
3. Use shared material parameters, a shared trajectory/extrinsic-error model, and acquisition-specific induced currents. Do not use independent per-frequency maps.
4. Separate training/tuning seeds, calibration of computational costs and thresholds, and final test seeds. No true parameter, true rank, or post-hoc threshold enters an algorithm.
5. Truth is used only for evaluation and explicitly labeled oracle baselines. All algorithms start from exactly the same nonoracle initial map and geometry.
6. Keep physical reduced-state information separate from free-current envelope information. A zero envelope score is not labeled physical nonidentifiability.
7. Record numerical-rank tolerance, absolute singular values, nuisance codimension, and uncertainty bounds. Saturated nuisance codimension implies exactly zero projected rank in the declared model.

## E1. Matched coherent versus phaseless efficient information

**Hypothesis.** The matched transformed experiment has no larger efficient map or pose information. Equality and strictness are predicted by the efficient-score measurability condition.

**Null for an algorithmic advantage.** The additional coherent information does not produce a measurable finite-sample pose improvement in the chosen regular regime.

**Minimal models.**

- $Y=a e^{ix}+\varepsilon$ with known positive $a$ and proper circular Gaussian noise; compare $Y$ to $|Y|^2$.
- $Y=a e^{i(x+\eta)}+\varepsilon$ with nuisance phase $\eta$; both efficient pose informations are zero.
- Zero-mean circular Gaussian scale family; intensity is sufficient for the scale score.
- Total field $Y=c+a e^{ix}+\varepsilon$ with a known reference $c$; intensity can retain scattered-phase information.
- One small physical full-wave scene with real map and pose parameters after gauge fixing.

**Likelihoods.** Use the exact induced magnitude/intensity law for independent complex Gaussian samples. With correlated noise, use the induced joint likelihood or an explicitly controlled integration method; do not silently replace it with independent additive Gaussian intensity noise. A separate direct-intensity sensor model is a mismatch control, not evidence against data processing.

**Controls.** Equal coherent parent noise realizations; fixed parameter support; identical nuisance family; singular Fisher support treatment; independent sensor-noise mismatch.

**Budget/seeds.** Analytic calculations for the scalar cases; 2,000 simulated observations for likelihood-score Monte Carlo checks using seeds 101-110, where necessary. These are proposed computations, not existing results.

**Metrics/plots.** Joint and efficient Fisher eigenvalues; minimum eigenvalue of the coherent-minus-phaseless difference; Monte Carlo integration uncertainty; equality residuals; error distributions and coverage on identifiable targets.

**Success/failure.** A statistically significant negative efficient-information difference under a genuinely matched regular model falsifies the implementation or proof assumptions. A strict joint difference combined with equal efficient information is an expected control, not failure. Coherent estimator superiority is a separate endpoint from the information theorem.

## E2. Nested rank, neutral admission, and numerical rank

**Hypothesis.** Fixed-linearization nested free-current spaces obey the exact loss formula. Neutral admissions preserve pose information. Map retention need not be monotone.

**Null for predictive usefulness.** The calculated local loss does not identify nonlinear calibration deterioration beyond what a raw sensing singular value already predicts.

**Minimal cases.**

- Equality: add a nuisance direction orthogonal to the existing visible pose tangent.
- Strict loss: add a nuisance direction aligned with one visible pose direction.
- Complete hiding: add all visible pose directions.
- Nonmonotone retention: $A=(1,1,0)^T$, $B=(1,0,1)^T$, $C=0$, then span$(e_1)$, then span$(e_1,e_2+e_3)$.
- Complex-coefficient case: require invariance under multiplication by $i$ for any proposed safe current subspace.
- Roundoff control: explicitly saturated nuisance range with all projected singular values at numerical noise level.

**Data/model.** Small exact algebraic examples followed by 12 fixed physical linearizations from nonoracle iterates of the common scalar model. Hold all physics, whitening, and gauges fixed during each rank sweep.

**Budget/seeds.** Deterministic exact examples; physical linearization seeds 201-212. No repeated nonlinear inversion is required for the algebraic theorem checks.

**Metrics/plots.** Eigenvalues of the predicted and measured information loss, equality residual, visible rank versus nuisance codimension, ordinary versus generalized spectra, free rank versus classical deterministic cutoff, hard projector distance versus soft-filter distance, and projected-current dimension versus safe-admission dimension.

**Success/failure.** Reject the implementation if the loss identity fails beyond its backward-error bound. Never count small relative ratios of roundoff singular values as physical rank. Nonlinear predictive failure does not falsify the fixed-linearization identity, but blocks a nonlinear method interpretation.

## E3. Dual spectra and the bias-information tradeoff

**Hypothesis.** Absolute pose strength and map retention predict their respective linear Gaussian error changes; adding a bias term explains failures missed by relative-only gates.

**Null.** Dual gating is redundant with sensing singular values and a residual-only rule for the predeclared outcomes.

**Minimal model.** Use the two-observation confounding model $y_1=h+c+\varepsilon_1$, $y_2=\epsilon c+\varepsilon_2$ and append an independent map observation. Sweep $\epsilon$ and a predeclared coefficient bound. This supplies both truncation-benefit and truncation-failure regimes.

**Physical extension.** Evaluate fixed current-cleaned map/pose tangents from E2 under full and limited aperture, with and without declared pose/map priors. Keep the task map basis fixed. Do not discard failed modes after examining their generalized eigenvalues.

**Controls.** Unit retention but weak absolute information; strong raw Fisher but nearly parallel map/pose ranges; wrong low-rank model with a small residual and large pose bias; gauge nulls; prior-only recovery distinguished from data recovery.

**Budget/seeds.** 1,000 linear Gaussian draws per small tangent, seeds 301-310. One optional nonlinear validation batch uses the same test grid as E4, without changing E4's primary endpoint.

**Metrics/plots.** Empirical covariance versus the appropriate fixed-tangent inverse information, directional map error versus retention, pose error versus absolute visible singular value, worst-case risk versus the analytic bias-plus-variance formula, coverage versus rank, and feasibility of simultaneous approximation/information constraints.

**Success/failure.** The fixed-linear estimator must match the analytic covariance and bias within Monte Carlo uncertainty. Superior nonlinear prediction requires a predeclared held-out comparison against a sensing-spectrum-only diagnostic. If relative gates pass a high-bias case, that falsifies the relative-only policy, as intended.

## E4. Primary matched-budget method comparison

### E4.1 Fixed problem family

The primary model is two-dimensional scalar Helmholtz scattering. The material is described by nine fixed smooth spatial basis functions with real material coefficients and a declared Ohmic-loss model. The free-space Green operators include the correct frequency factors for a fixed material. Use an inverse grid of 16 by 16 cells and a finer data-generation grid of 32 by 32 cells. Grid refinement is an evaluation control, not a new unknown parameter.

Use three acquisition poses, with the first externally anchored. A shared three-dimensional rigid-array error parameter affects the remaining poses. Each pose has two illuminations and twelve receiving channels. The frequency ratios are fixed as $\{1/4,1/2,3/4,1\}$ relative to the highest frequency. The core study has independent proper complex receiver noise with a predeclared SNR. Two contrast ranges and two aperture configurations are fixed before tuning, rather than selected after successful reconstructions.

This low-dimensional model is intended to establish or reject the mechanism economically. General high-dimensional images and measured Maxwell data are separate publication gates.

### E4.2 Primary endpoint

The single primary endpoint is the difference in **joint successful-recovery probability at a fixed computational budget**, PRASC minus a tuned direct joint full-wave baseline, averaged over the fixed initialization grid and test seeds.

A run succeeds only if all of the following hold:

- anchored pose error meets a predeclared wavelength/resolution tolerance;
- map error on a fixed resolution-limited task representation meets a noise-derived tolerance;
- full-physics prediction is consistent with independent validation measurements;
- no required target is falsely labeled calibrated and no unhandled gauge remains.

Optimizer success/status is not a success criterion. Failure to finish within budget counts as failure.

### E4.3 Physical task metrics and noninferiority margins

Let $\lambda_{\min}$ be the shortest wavelength, $R_{\rm eff}$ a fixed array lever arm, and use the pose metric $\|h\|_x^2=\|\delta t\|^2+R_{\rm eff}^2\delta\theta^2$. Define the resolved map representation by a fixed low-pass operator $Q_{\rm res}$ with width no smaller than $\max\{2\Delta,\lambda_{\min}/(2\,\mathrm{NA})\}$, where NA is the declared geometric aperture factor, not a fitted image property.

Before final test data are generated, use a single declared nominal reference material and the actual design/noise model to compute full-rank reference task covariance matrices $C_x^{\rm ref}$ and $C_\chi^{\rm ref}$. These are design scales, not guarantees of error in the nonlinear test scenes. If the chosen task is not supported even in this reference experiment, revise the preregistration before test generation, not after looking at outcomes.

Use a pose success tolerance $d_x=\min\{\lambda_{\min}/16,\ell_{\rm res}/4\}$ in the lever-arm metric. Use a map success tolerance $d_\chi=2\sqrt{\operatorname{tr}(Q_{\rm res}C_\chi^{\rm ref}Q_{\rm res}^T)/q_T}$, where $q_T$ is the fixed number of task coefficients. Set the pose noninferiority margin to $\Delta_x=\min\{\lambda_{\min}/32,\sqrt{\operatorname{tr}C_x^{\rm ref}/p}\}$ and the map noninferiority margin to $\Delta_\chi=d_\chi/2$.

These numerical constants are preregistered engineering choices tied to phase budget, resolution, and noise. They are not universal physical laws. Report the resulting margins in meters/degrees and permittivity units before running final tests. Also report full-resolution map error as secondary information, without redefining the primary task.

### E4.4 Seeds and initialization grid

Use seeds 1-10 only for tuning and computational-cost calibration. Lock all settings before using final scene/noise seeds 1001-1020.

For each test seed, use twelve initialization errors: three radii $\{\lambda_{\min}/8,\lambda_{\min}/2,\lambda_{\min}\}$ and four uniformly spaced translation directions. Split each error radius equally in the squared pose metric between translation and lever-arm rotation. Alternate the rotation sign deterministically with the seed. This gives 240 paired primary runs per method, with 20 independent seed clusters. Use the same initialized map for all methods, obtained from the same declared low-frequency nonoracle procedure.

### E4.5 Budget and baselines

Use a primary cap of 200 forward/adjoint-equivalent solves per run. Count every full-wave right-hand-side solve, whether used for a forward field, adjoint, parameter derivative, certificate, or acceptance test. Convert reduced solves and extra operator products to equivalent work using a calibration fixed on the tuning set and report raw operator products as well. Charge SVD, projector updates, factorizations, and their reuse to wall time; do not describe them as free. Publish both equivalent-work and wall-time curves.

Mandatory primary baseline: direct joint full-wave inversion with the same physical likelihood, priors, gauge, initial candidates, frequency data, and optimizer/tuning budget. Permit it the same frequency continuation and standard preconditioning. A baseline denied these tools would not isolate a SOM-specific effect.

Secondary baselines: fixed-rank coherent SOM-style approximation, the corresponding tuned phaseless method using the matched parent intensity, and an independently implemented source/receiver-extension comparator if its physical parameterization is matched. Include known-map pose-only and known-pose map-only oracles as diagnostics, not as competitors used for tuning the proposed method.

An optional particularly strong ablation gives the direct solver the proposed acquisition policy. If the benefit transfers fully, report it as a general acquisition-control gain rather than an exclusive SOM gain.

### E4.6 Statistical decision and multiplicity

Preserve pairing at each seed and initialization. Estimate uncertainty with a seed-cluster bootstrap (10,000 resamples of the 20 seeds, retaining all twelve initializations together). Report interval sensitivity because 20 independent clusters may give limited power; do not equate a wide inconclusive interval with equivalence.

For the compound primary claim, require all three one-sided simultaneous intervals, using Bonferroni levels $1-0.05/3$:

1. the lower confidence bound for the success-probability difference is strictly positive;
2. the upper confidence bound for the difference in pose RMSE is below $\Delta_x$;
3. the upper confidence bound for the difference in task-map RMSE is below $\Delta_\chi$.

Apply Holm correction to the separately enumerated secondary baseline comparisons. A failed primary endpoint cannot be rescued by selecting a favorable contrast, frequency, aperture, seed subset, or runtime statistic after the fact. Treat all methods' failed runs consistently; report failure-conditioned errors only as secondary diagnostics alongside the unconditional primary endpoint.

### E4.7 Adversarial controls

**Limited-aperture false basin.** Use a prespecified narrow aperture and initial errors covering multiple high-frequency phase branches. Run the known-map pose-only oracle from the same initial poses. Report cases where an incorrect local minimum has a small residual, high Fisher curvature, and a passed covariance-based frequency surrogate. Such cases test the distinction between a local certificate and a basin claim.

**Non-pose phase corruption.** Add unmodeled clock delay, producing a frequency-dependent phase ramp, or a structured mutual-coupling perturbation. Predeclare the corruption levels on the tuning set and do not estimate them as pose. The pose-only PRASC implementation is expected to lose in some of these cases, potentially to a phaseless estimator. This is failure of a misspecified coherent algorithm, not reversal of matched-experiment information ordering. An appropriately augmented coherent model is a separate control, subject to its own observability/gauge analysis.

**Truncation bias.** Include a scene with a controlled material component outside the low-rank representation. Test whether approximation certificates refuse the low-rank step. A favorable pose Gram accompanied by a wrong physical map is failure.

### E4.8 Plots and decision

Publish paired basin maps, unconditional joint success rates and intervals, pose/map error versus charged work, failure counts by certificate type, rank/frequency trajectories, full-physics residual and held-out predictive residual, envelope versus physical visibility, and the bias bound along accepted steps.

A method claim requires the preregistered compound primary gate to pass. Without it, retain the theorem/diagnostic result and withdraw performance superiority. A scalar nine-parameter success does not establish full-image or 3D Maxwell superiority.

## E5. Shared-map acquisition, source symmetry, and gauge

**Hypothesis.** Shared-map compensation inconsistency creates the innovation predicted by the Schur update. The innovation/rank-loss budget predicts first-order safety of combined acquisition/rank changes.

**Null.** The proposed acquisition score does not improve final reconstruction over a random or fixed acquisition schedule at the same charged acquisition and inversion budget.

**Minimal experiment.** Use $a_1=a_2=1$, $b_1=1$, $b_2=-1$, then compare to duplicated $b_2=b_1$. Verify the exact innovation and the failure of adding independently map-profiled per-frame Grams. Add an explicit nuisance-saturation case with zero cleaned rows.

**Physical model.** Reuse E4 without extra material complexity. Enumerate eight candidate source/receiver/illumination actions. Compare the exact shared-map criterion with an incorrectly independent-map criterion, one-step greedy design, pair look-ahead, exhaustive design over this small candidate set, and random selection.

**Gauge/symmetry controls.** Remove the anchor and apply a joint rigid map/trajectory transform; acquisitions must not be claimed to recover absolute pose. Separately compare one isotropic source against separated/directional illuminations, recording only the source-stabilizer directions removed. Restore the anchor for absolute-pose evaluation.

**Budget/seeds.** Seeds 401-412 for tangent tests; reuse the fixed E4 seeds for nonlinear secondary tests. Exhaustive matrix scoring is allowed for eight candidates but is charged if it requires new forward/adjoint computations. No universal frame-count claim is inferred from this finite candidate set.

**Metrics/plots.** Hidden dimension, smallest task-scaled visible singular value, common-compensator disagreement, exact versus predicted innovation, rank-loss versus acquisition-gain eigenvalues, map survival, and final reconstruction under equal budgets.

**Success/failure.** The exact algebra must match within backward-error bounds. A design that changes rank but not conditioning is not credited with precision improvement. Failure of the physical candidate set to span missing signatures is reported as an acquisition limitation, not repaired by an oracle action. Greedy submodularity is not assumed; the complementary pair example should violate it.

## Publication evidence boundary

For a TAP calibration claim, add a small controlled 3D dyadic-Maxwell or measured-array validation containing known phase-center offsets and at least one clock/coupling mismatch control. For a TGRS imaging-method claim, add a genuinely imaging-relevant task with a larger map representation and matched-budget comparative evidence. These are scientific evidence gates, not statements of official journal acceptance rules. Neither is replaced by the successful completion of elementary matrix checks.

## Additional locked evaluation conventions

Reserve an independent noisy replicate at the same acquisition geometry for final predictive validation. It is not used to tune ranks, thresholds, priors, or frequency schedules. The cost and existence of this replicate are the same for all methods. Predictions are evaluated against the known replicate-noise model; any fitted-data residual degrees-of-freedom approximation is reported separately rather than treated as an exact chi-square law.

Every run returns its last finite, physically admissible parameter iterate even when an algorithm declares failure. Unconditional error statistics use that iterate. A run with no admissible iterate receives the maximum error allowed by the predeclared compact parameter domain; the domain and this rule are locked before test generation. Such a run is always a primary-endpoint failure. No failure is silently excluded from RMSE or confidence intervals.
