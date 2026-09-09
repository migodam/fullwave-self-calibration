# Family 4: frequency stacking and trajectory geometry

Date: 2026-09-03 (SGT; UTC stamp in `results/family4_results.json`).
Experiment: `experiment_pose_confounding_spectral_geometry`.

Family 4 reuses, without rerunning, the corrected Family 1 forward/Jacobian
pairs (`src/helmholtz.py`), the Family 1 pose/scene construction, and the
Family 2 smooth p=24 RBF basis and machine-rank rule.  It does not repeat the
earlier family gates.  All numbers below are for the discrete, finite-
dimension N=16 whitened/realified model.

## Exact command and runtime

```bash
.venv/bin/python src/family4_frequency_trajectory.py
```

Wall runtime: 0.78 s.  Platform:
macOS-26.6.2-arm64-arm-64bit, Python 3.12.13,
numpy 2.5.2, scipy 1.18.1,
matplotlib 3.11.1.  Deterministic dense linear
algebra; no RNG used.

Scenario: N=16, T=6 poses, n_rx=4, m_f=48 per frequency, q=18 pose columns,
two-blob chi0 (amp 0.3/0.5, sigma 0.09/0.07 at (-0.15,-0.12)/(0.18,0.14)),
p=24 smooth basis (4x6 Gaussian RBFs on [-0.3,0.3]^2, sigma_b=0.16, unit
columns).  Frequency stacks use the Family 1 90-degree arc at radius 1.6 for
every frequency.  Rank tolerance:
`tol(M) = max(M.shape) * eps_machine * sigma_1(M)`.

## Tolerances used (Family 4 gates)

| gate | value |
|---|---:|
| A min_eig(D) for each prefix increment | >= -1e-10 * max(1, ||larger||_2) |
| B no-prior rel Frobenius K_IS / K_SLAM | < 1e-10 |
| B no-prior max abs rho difference | < 1e-10 |
| B fixed-prior change | expected clearly nonzero (recorded, e.g. > 1e-3) |
| B jointly-scaled-prior rel Frobenius and rho diff | < 1e-10 |
| C duplicate rho movement | < 1e-10 |
| C distinct rho movement | > 1e-4 in at least 2 of 3 directions |
| D trajectory ordering | `hypothesis_holds` boolean; no forced pass |

## Noise / whitening convention and normalization control

W = None (identity noise): whiten_realify(A,B,None) returns sqrt(2)*[Re; Im] row stacks; raw identity-noise stacks throughout. Frequencies and trajectories are NOT block-energy or SNR normalized.

The frequency-diversity and trajectory magnitudes in this record therefore
reflect the raw identity-noise stack.  A declared
block-normalized/SNR-matched control is **deferred** (deferred and declared: raw identity-noise stacking is the tested metric; a block-normalized/SNR-matched control is not silently conflated with it);
the raw result is not silently promoted to a same-energy conclusion.
Frequency row energy differs with f (see the table below), and trajectory
rows also have different norms because equal path length forces different
radii.

### Frequency block norms (raw identity-noise smooth stacks)

| f | k_b=2*pi*f | ||A_s||_F | ||A_s||_2 | A_F/||A_1||_F | ||B_R||_F | ||B_R||_2 | B_F/||B_1||_F |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.0 | 6.283185 | 3.217209e-02 | 1.825708e-02 | 1.0000 | 2.006294e-01 | 1.029932e-01 | 1.0000 |
| 1.4 | 8.796459 | 2.500301e-02 | 2.011244e-02 | 0.7772 | 2.586319e-01 | 1.255715e-01 | 1.2891 |
| 1.8 | 11.309734 | 1.676217e-02 | 1.173556e-02 | 0.5210 | 1.658858e-01 | 9.085180e-02 | 0.8268 |

## Claim-status table

| check | status | executed comparison | key numbers |
|---|---|---|---|
| A: PSD monotonicity of stacking distinct frequencies | **PASS** | min eig of both K_SLAM and K_eff prefix increments | K_SLAM F1->2 min_eig=1.517e-13; K_eff F1->2 min_eig=1.288e-14; K_SLAM F2->3 min_eig=1.149e-13; K_eff F2->3 min_eig=5.350e-14 |
| B: duplicate no-prior invariance | **PASS** | rel Frobenius + rho | relKIS=1.080e-15, relKSL=2.893e-13, max|rho|=1.185e-12 |
| B: fixed finite prior breaks duplicate invariance (recorded, expected) | **observed** | max rho and rel Frobenius | max|rhoX dup - rhoX single|=3.984254e-02; rel F K_eff=1.048203e-01 |
| B: jointly scaled prior restores invariance | **PASS** | rel Frobenius + rho | relF=1.037e-15, max|rho|=1.221e-15 |
| C: distinct frequency moves confounded directions | **PASS** | rho movement per direction | 3/3 directions > 1e-4 |
| C: duplicate frequency leaves rho fixed | **PASS** | rho movement per direction | max movement = 2.900e-16 < 1e-10 |
| D: retained_mass circle > arc180 > arc90 > straight | **FAIL (counterexample, not forced)** | observed order | arc90, straight, circle360, arc180 |
| D: confusable_mass circle < arc180 < arc90 < straight | **FAIL (counterexample, not forced)** | observed order | arc180, circle360, straight, arc90 |

### A raw rows

| comparison | min eig D | tol | gate |
|---|---:|---:|---|
| K_SLAM F=1->2 | 1.516646e-13 | 1.000e-10 | PASS |
| K_eff F=1->2 | 1.288406e-14 | 1.000e-10 | PASS |
| K_SLAM F=2->3 | 1.148563e-13 | 1.000e-10 | PASS |
| K_eff F=2->3 | 5.350473e-14 | 1.000e-10 | PASS |

### B duplicate controls

No-prior rel K_IS = 1.080e-15,
rel K_SLAM = 2.893e-13,
max|rho_dup - rho_1| = 1.185e-12.
Fixed prior (alpha=1): max|rho_X| = 3.984254e-02,
rel F K_eff = 1.048203e-01.
Joint scaled prior ((1+c^2)*alpha*I): rel F = 1.037e-15,
max|rho| = 1.221e-15.

### C three most confounded directions

| dir | rho_single | rho_distinct | movement_distinct | z-resid distinct | rho_dup | movement_dup | z-resid dup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 5.510027e-11 | 1.212691e-01 | 1.212691e-01 | 3.482372e-01 | 5.510021e-11 | -6.566674e-17 | 7.422945e-06 |
| 1 | 1.414337e-10 | 2.239175e-01 | 2.239175e-01 | 4.731992e-01 | 1.414340e-10 | 2.899515e-16 | 1.189259e-05 |
| 2 | 1.095607e-09 | 4.768305e-01 | 4.768305e-01 | 6.905291e-01 | 1.095607e-09 | 1.328091e-16 | 3.309995e-05 |

### D equal-budget trajectory table

| trajectory | retained_mass | confusable_mass | theta_min deg | log_volume | count rho<1e-6 | rho_min |
|---|---:|---:|---:|---:|---:|---:|
| straight | 8.353517 | 15.646483 | 9.887569e-05 | -204.0601 | 6 | 2.9781e-12 |
| arc90 | 9.038781 | 14.961219 | 2.561321e-05 | -220.4378 | 7 | 1.9984e-13 |
| arc180 | 7.999302 | 16.000698 | 3.193125e-04 | -175.8303 | 6 | 3.1059e-11 |
| circle360 | 8.139513 | 15.860487 | 2.029155e-02 | -133.4532 | 3 | 1.2543e-07 |

hypothesis_holds = **False**.  The requested order was
circle360 > arc180 > arc90 > straight for retained_mass (and the same order
inverted for confusable_mass).  Observed: retained = arc90, straight, circle360, arc180;
confusable = arc180, circle360, straight, arc90.

### D geometry and row-norm record (audit rule 6)

| trajectory | design L | open polyline | closed polyline (if counted) | min range | max range | ||A||_F | ||A||_2 | mean A row norm | ||B||_F | mean B row norm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| straight | 5.026548 | 5.026548 | 10.053096 | 1.677099 | 2.979353 | 2.39010e-02 | 1.38025e-02 | 3.36266e-03 | 1.48196e-01 | 1.90108e-02 |
| arc90 | 5.026548 | 5.005903 | 9.531386 | 3.200000 | 3.200000 | 1.60100e-02 | 8.94906e-03 | 2.30568e-03 | 9.77199e-02 | 1.23833e-02 |
| arc180 | 5.026548 | 4.944272 | 8.144272 | 1.600000 | 1.600000 | 3.18323e-02 | 1.61849e-02 | 4.57159e-03 | 1.72253e-01 | 1.84989e-02 |
| circle360 | 5.026548 | 4.000000 | 4.800000 | 0.800000 | 0.800000 | 6.73768e-02 | 3.23981e-02 | 9.62787e-03 | 4.05789e-01 | 5.33797e-02 |

The continuous design length is L = pi*1.6 for every trajectory.  The sampled
open polyline is shorter for non-straight trajectories; for `circle360`
(endpoint=False) a full circumference is attained only if the closing
last-to-first segment is included, and that sampled chord polygon is recorded
above.  Raw row norms also differ between trajectories, so the ranking below
is a scenario result, not a causal angular-coverage theorem.

## Cannot-establish section

* All claims are finite-dimensional and model-specific: they are statements
  about the discrete N=16 whitened/realified Jacobians, not continuum-limit,
  exact-global-SE(2), recovery, or estimator claims.
* No claim is made that the trajectory ordering is universal.  The observed
  equal-budget ordering violates the stated hypothesis (arc90 has the
  largest retained mass), and the result is recorded as a counterexample,
  not forced to pass.
* Equal total path length forces different radii (straight offset 1.6;
  arc90 R=3.2; arc180 R=1.6; circle R=0.8) and different absolute
  measurement locations, so a mass ordering mixes arc curvature, radius,
  and endpoint effects.  The table does not isolate a single causal factor.
* Frequency stacking invariances hold for the c-scaled duplicate of one
  physical frequency block and the same finite-prior/whitening model; they
  do not prove invariance under arbitrary resampling, noise whitening, or
  model mismatch.
* The tested stack uses the raw identity-noise metric (W=None).  A declared
  block-normalized/SNR-matched control is deferred, so the raw
  frequency-diversity magnitudes (including check C movements) may partly
  reflect per-frequency block energy; no same-SNR claim is made.
* Check A establishes absolute matrix monotonicity only.  Normalized
  retention spectra are allowed to move non-monotonically with F because
  K_IS changes; no retention-spectrum decrease was interpreted as a
  contradiction.
* Check C directions are the reconstructed parameter-space vectors
  u = V_A diag(1/s_A) w (audit rule 3), not the raw 24-vector eigenvectors
  of R_op.
* The equal-budget geometry record mixes arc curvature, radius/standoff, and
  per-pose signal energy; the ordering is a reproducible scenario result or
  counterexample and cannot establish that angular coverage alone creates the
  observed order (audit rules 6-7).
* The 3 confounded directions are those of the single-frequency block;
  movements for other directions or larger stacks were not exhaustively
  checked.

## Artifacts

- results: `results/family4_results.json`
- figures: family4_monotonicity.png, family4_duplicate.png,
  family4_shared_z.png, family4_trajectories.png
- this report: notes/family4_report.md

Self-cell formula used: self-cell v2 corrected 2026-09-03: complete equal-area disk integral (adds the -1/k_b^2 lower-endpoint term omitted by v1).

## Artifacts and digests

```text
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242  src/family4_frequency_trajectory.py
c05c77a1402dab3b291c2e2ee3eda3ad9e2ef886fe05e4a7961d4dadd56c8b97  results/family4_results.json
0a72c2f4d80efec0671fa89e49532af182f33cc81acd03cef6e0ca6e9895a86e  figures/family4_monotonicity.png
7c10ea0dcfe3555f25769c7332ad56c67280ba9f7e522132ad6dd73ec49cbd70  figures/family4_duplicate.png
bac13a8c7c1cf52b98e18011a8333e669d980f09beaef6322ba571485d3f425f  figures/family4_shared_z.png
3ffaf43e161302e8132d87e8155c5437587948bad5600e22c3e01dadab91283e  figures/family4_trajectories.png
```
