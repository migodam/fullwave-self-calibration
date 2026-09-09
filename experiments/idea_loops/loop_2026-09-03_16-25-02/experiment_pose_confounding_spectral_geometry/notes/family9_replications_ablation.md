# Family 9: seed replications and component ablations

Date: 2026-09-03T11:58:36+00:00 UTC.  Experiment:
`experiment_pose_confounding_spectral_geometry`.

Family 9 re-runs the deterministic algebraic-spine (Family 2 c1-c3) and
frequency-diversity (Family 4 check C) claims over seven scenes and adds
descriptive component ablations: pose-DOF nuisance subsets, receiver count,
pose count, and smooth-basis p.  It reuses (read-only) `helmholtz.py`,
`family1_pilot.py`, `family2_algebraic_spine.py`, and
`family4_frequency_trajectory.py`; no existing file is modified.

## Exact command and runtime

```bash
.venv/bin/python src/family9_replications_ablation.py
```

Wall runtime: 1.35 s (part seconds:
A 0.51, B 0.00,
C 0.09, D 0.07,
E 0.00).  Platform:
macOS-26.6.2-arm64-arm-64bit, Python 3.12.13,
numpy 2.5.2, matplotlib
3.11.1.

Source SHA-256 (this script):
`a003d3b3067995d5201a9a892c440880669cc27bf3330980bebc344f8db78803`

Reused-module SHA-256:
`family1_pilot.py` `762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7`, `family2_algebraic_spine.py` `173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235`, `family4_frequency_trajectory.py` `ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242`, `helmholtz.py` `ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585`

## Config

N=16 (N^2=256), standard T=6, n_rx=4
(realified rows m=2*T*n_rx; m=48 at standard T/n_rx), k_b=2*pi (f=1.0),
Family-1 90-degree arc radius 1.6, phi
-45.0..45.0 deg,
theta=atan2(-py,-px); rx offsets [[-0.06, 0.0], [0.06, 0.0], [0.0, -0.06], [0.0, 0.06]]; tx offset
[0.0, 0.0].  Standard two-blob chi0 as in Families 1-8.  Smooth basis
p=24 (4x
6 centres, sigma_b=
0.16, unit columns).  Random scenes:
8 positive blobs per scene,
seeds [1001, 1002, 1003, 1004, 1005], rescaled to the two-blob L2
norm.  Whitening is identity (`whiten_realify(A,B,None)`); alpha=
1.0.  Ranks use the family-2 rule
`tol(M)=max(M.shape)*eps*sigma_1(M)`.

## Claim-status table (seed replications over 7 scenes)

| claim                               | status  | executed comparison                              | key numbers                                          |
| ----------------------------------- | ------- | ------------------------------------------------ | ---------------------------------------------------- |
| rank identity (F2 c2)               | PARTIAL | integer residual == 0 per scene                  | 6/7; ring residual = 1 (machine-rank boundary)       |
| no-prior retention identity (F2 c3) | PASS    | max abs diff < 1e-10*max(1,||A_s||_2^2)          | 7/7; diff ~1e-15 on all scenes                       |
| kernel identity (F2 c1)             | PARTIAL | F2 subspace gate; identity_support/nullity match | F2 gate 7/7; identity_support 6/7; nullity match 6/7 |
| frequency diversity (F4 check C)    | PASS    | rho movement > 1e-4 in >=2 of 3 directions       | 7/7; per-scene moved counts [3, 3, 3, 3, 3, 3, 3]    |

Part B-D-E claims are descriptive ablations with no forced pass:

| part | ablation                                         | status   |
| ---- | ------------------------------------------------ | -------- |
| B    | pose-DOF nuisance subsets + marginal destruction | recorded |
| C    | receiver count n_rx in {1,2,4,8}                 | recorded |
| D    | pose count T in {3,6,12}                         | recorded |
| E    | smooth-basis p in {12,24,36}                     | recorded |

## A. Seed replications (per scene)

Per-scene table (rank-identity residual; retention identity max abs diff and
gate; kernel identity nullities KSL/C and nullity match; frequency-diversity
moved count and family-4 gate):

| scene     | rank res | rank pass | retention diff | retention pass | null KSL/C | null match | moved >1e-4 | C pass |
| --------- | -------- | --------- | -------------- | -------------- | ---------- | ---------- | ----------- | ------ |
| seed_1001 | 0        | PASS      | 1.912e-15      | PASS           | 0/0        | PASS       | 3           | PASS   |
| seed_1002 | 0        | PASS      | 1.332e-15      | PASS           | 0/0        | PASS       | 3           | PASS   |
| seed_1003 | 0        | PASS      | 1.443e-15      | PASS           | 0/0        | PASS       | 3           | PASS   |
| seed_1004 | 0        | PASS      | 2.071e-15      | PASS           | 0/0        | PASS       | 3           | PASS   |
| seed_1005 | 0        | PASS      | 2.817e-15      | PASS           | 0/0        | PASS       | 3           | PASS   |
| two_blob  | 0        | PASS      | 1.554e-15      | PASS           | 0/0        | PASS       | 3           | PASS   |
| ring      | 1        | FAIL      | 1.987e-15      | PASS           | 1/0        | FAIL       | 3           | PASS   |

Pass-rate summary: rank identity 6/
7; retention identity 7/
7; kernel identity F2 gate 7/
7 with identity_support 6/
7 and nullity match 6/
7; frequency diversity 7/
7 (all scenes move 3/3 directions; see table).

Cross-check of the standard two-blob row against stored results:
0
max abs movement difference vs `family4_results.json`
(match < 1e-12:
True).

## A caveat: ring-scene rank boundary

The ring scene (radially very smooth) has K_SLAM's smallest eigenvalue below
the family-2 eig rank tolerance while C=(I-ZZ^T)A_s stays full rank at its own
(larger) SVD tolerance.  The literal machine-rank identity therefore records
residual 1 and nullity match false on the ring scene.  This is recorded as a
machine-rank boundary event, not as evidence that the exact algebraic identity
is false; the exact statements would hold under exact arithmetic for these
subspaces if no near-degeneracy crosses the tolerance.

## Check-C direction method note

The executed directions follow Family 4 check C exactly:
`family4.generalized_eigen_directions(A1, P1)` with P1 = I - Z Z^T for
Range(B1).  The literal brief form `A1.T@A1` is dimension-incompatible with
the verified API (R_op = Q_A^T W Q_A requires a data-space W; A1.T@A1 is
p x p, Q_A is m x r).  Diagnostic raised:
ValueError:
matmul: Input operand 1 has a mismatch in its core dimension 0, with gufunc signature (n?,k),(k,m?)->(n?,m?) (size 48 is different from 24)
Dimensions:
{'A1': (48, 24), 'A1.T@A1': (24, 24), 'Q_A': (48, 24)}.

## B. Pose-DOF nuisance baseline

Standard two-blob scene, f=1.0, full B_R has q=18 columns ordered
[px, py, theta] per pose.  Known-pose reference = zero B (retained mass
24.000000, equal to p=24 up to roundoff).
The known-pose row is stored as a full-width zero B (q=18) because
`family2.range_basis`/`rank_svd` need at least one column.

| subset     | q_sub | retained  | confusable | theta_min deg | rho_min   | K_IS fro   | K_SLAM fro | K_eff fro  |
| ---------- | ----- | --------- | ---------- | ------------- | --------- | ---------- | ---------- | ---------- |
| known_pose | 18    | 24.000000 | 0.000000   | 90.0000       | 1.000e+00 | 4.6217e-04 | 4.6217e-04 | 4.6217e-04 |
| x_only     | 6     | 18.020102 | 5.979898   | 0.0308        | 2.882e-07 | 4.6217e-04 | 2.2127e-04 | 4.6082e-04 |
| y_only     | 6     | 18.028847 | 5.971153   | 0.0098        | 2.936e-08 | 4.6217e-04 | 2.3039e-04 | 4.6148e-04 |
| theta_only | 6     | 21.362289 | 2.637711   | 5.0796        | 7.839e-03 | 4.6217e-04 | 4.6124e-04 | 4.6217e-04 |
| xy         | 12    | 12.572748 | 11.427252  | 0.0013        | 5.197e-10 | 4.6217e-04 | 7.5686e-05 | 4.6014e-04 |
| x_theta    | 12    | 15.382237 | 8.617763   | 0.0035        | 3.771e-09 | 4.6217e-04 | 2.2081e-04 | 4.6082e-04 |
| y_theta    | 12    | 15.378367 | 8.621633   | 0.0038        | 4.357e-09 | 4.6217e-04 | 2.2913e-04 | 4.6148e-04 |
| full       | 18    | 9.407162  | 14.592838  | 0.0004        | 5.510e-11 | 4.6217e-04 | 2.4109e-06 | 4.6014e-04 |

Marginal destruction (loss of retained mass / increase of confusable mass
when a nuisance subset is added):

| transition            | delta retained | delta confusable |
| --------------------- | -------------- | ---------------- |
| add x nuisance        | 5.979898       | 5.979898         |
| add y nuisance        | 5.971153       | 5.971153         |
| add theta nuisance    | 2.637711       | 2.637711         |
| add y to x-only       | 5.447354       | 5.447354         |
| add theta to x-only   | 2.637866       | 2.637866         |
| add x to y-only       | 5.456099       | 5.456099         |
| add theta to y-only   | 2.650481       | 2.650481         |
| add x to theta-only   | 5.980053       | 5.980053         |
| add y to theta-only   | 5.983923       | 5.983923         |
| add y+theta to x-only | 8.612940       | 8.612940         |
| add x+theta to y-only | 8.621685       | 8.621685         |
| add x+y to theta-only | 11.955127      | 11.955127        |
| add theta to xy       | 3.165586       | 3.165586         |
| add y to x-theta      | 5.975075       | 5.975075         |
| add x to y-theta      | 5.971205       | 5.971205         |

Pattern observation: single-DOF retained mass is highest for theta alone
(21.36), then y-only (18.03) and x-only (18.02).  Adding theta to x-only or
y-only costs only ~2.64-2.65 retained mass, while adding x or y to
theta-only costs ~5.98 each, and adding the missing translational DOF to
x/y pairs costs ~5.45-5.46; full nuisance leaves retained = 9.41.  These are
recorded observations, not forced claims.

## C. Receiver-count ablation

| n_rx | retained | confusable | theta_min deg | rho_min    | B/A fro |
| ---- | -------- | ---------- | ------------- | ---------- | ------- |
| 1    | 0.000000 | 12.000000  | 0.0000        | -8.882e-16 | 6.2373  |
| 2    | 6.000000 | 18.000000  | 0.0000        | -4.441e-16 | 6.2376  |
| 4    | 9.407162 | 14.592838  | 0.0004        | 5.510e-11  | 6.2361  |
| 8    | 9.399363 | 14.600637  | 0.0009        | 2.684e-10  | 6.2361  |

Recorded monotonicity flags: retained non-increasing =
False; confusable
non-decreasing = False.

Caveat: r_A itself changes with measurement count: n_rx=1 gives r_A=12 and Range(B_R)=R^12 (A fully B-confounded, retained 0); n_rx=2 gives r_A=24 with Range(B_R) contained in Range(A_s) (retained 6); n_rx>=4 keeps r_A=24 with a transverse intersection (retained 9.41 for n_rx=4 and 9.40 for n_rx=8). The literal monotonicity flags over the mixed-r_A rows are therefore not interpretable as a pure receiver-count effect.

## D. Pose-count ablation

| T  | q  | retained  | confusable | theta_min deg | rho_min    |
| -- | -- | --------- | ---------- | ------------- | ---------- |
| 3  | 9  | 15.000000 | 9.000000   | 0.0000        | -1.332e-15 |
| 6  | 18 | 9.407162  | 14.592838  | 0.0004        | 5.510e-11  |
| 12 | 36 | 6.579120  | 17.420880  | 0.0028        | 2.445e-09  |

Recorded monotonicity flags: retained non-increasing =
True; confusable
non-decreasing = True.

Note: q=3T grows with T and r_A stays 24 (A_s has p=24 columns), so retained mass falls monotonically here (15 -> 9.41 -> 6.58) as more pose nuisance columns occupy more of the data space.

## E. Smooth-basis ablation

| p  | x_centers | y_centers | retained  | confusable | theta_min deg | rho_min    |
| -- | --------- | --------- | --------- | ---------- | ------------- | ---------- |
| 12 | 3         | 4         | 0.157018  | 11.842982  | 0.0510        | 7.912e-07  |
| 24 | 4         | 6         | 9.407162  | 14.592838  | 0.0004        | 5.510e-11  |
| 36 | 4         | 9         | 18.749631 | 17.250369  | 0.0000        | -4.441e-16 |

Recorded monotonicity flags: retained non-increasing =
False; confusable
non-decreasing = True.

Note: B_R_standard is fixed (r_B=18) while r_A = p grows, so retained mass increases with p (0.157 -> 9.41 -> 18.75) partly because more A directions are added outside the fixed B range; confusable mass also increases (11.84 -> 14.59 -> 17.25).

## Figures

* [figures/family9_seed_replications.png](figures/family9_seed_replications.png) -
  (a) Family 2 core-identity residuals per scene; (b) frequency-diversity
  movements for the three most-confounded directions.
* [figures/family9_pose_dof_ablation.png](figures/family9_pose_dof_ablation.png) -
  retained and confusable mass bars for pose-DOF subsets including the
  known-pose reference.
* [figures/family9_rx_T_basis_ablation.png](figures/family9_rx_T_basis_ablation.png) -
  receiver-count, pose-count, and basis-p ablation panels.

## Cannot establish

* These are finite-dimensional N=16 numbers on one deterministic pose
  geometry at f in {1.0,1.4}; no continuum theorem, other-N transfer, other
  geometry, other contrast, or universal bound is asserted.
* The ring scene's rank-identity/nullity mismatch is a machine-rank
  tolerance boundary; it neither disproves the exact algebraic identities nor
  certifies them in floating point.
* Kernel identity on scenes with empty null spaces is a vacuous (though
  consistent) check at the declared eig thresholds; it does not test a
  nontrivial kernel.
* Part B/C/D/E trends are descriptive; the recorded monotonicity flags are
  observations over the tested grid, not guarantees.
* Frequency-diversity "movement" is defined on directions fixed by the
  single-frequency kernel and measured with raw identity-noise stacks; block
  or SNR normalization was neither applied nor silently conflated.
