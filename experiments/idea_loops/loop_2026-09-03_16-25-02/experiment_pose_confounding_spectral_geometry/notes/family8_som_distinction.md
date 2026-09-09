# Family 8: current-space G_S/SOM versus map-tangent A/K_eff distinction

Date: 2026-09-03T11:49:23+00:00 UTC.  Experiment:
`experiment_pose_confounding_spectral_geometry`.

Family 8 is a NEW sanity-check family.  It reuses (and does not modify)
`helmholtz.py`, `family1_pilot.py`, `family2_algebraic_spine.py`, and
`family4_frequency_trajectory.py`.  It makes the current-space sensing map
`G_S` (stacked over poses) and the map-tangent Jacobian `A_c` / pixel kernel
`K_eff` explicit objects in the same N^2 pixel space, and verifies that the
algebraic relationships hold while the two mode subspaces remain distinct.

## Exact command and runtime

```bash
.venv/bin/python src/family8_som_distinction.py
```

Wall runtime: 0.45 s.  Platform:
macOS-26.6.2-arm64-arm-64bit, Python 3.12.13,
numpy 2.5.2, scipy 1.18.1,
matplotlib 3.11.1.

Source SHA-256 (this script):
`13157b2ac3b916f18ec937f9c15dab7880ef98eb1f979ddba7b79ab3fe02b198`

Reused-module SHA-256:
`family1_pilot.py` `762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7`, `family2_algebraic_spine.py` `173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235`, `family4_frequency_trajectory.py` `ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242`, `helmholtz.py` `ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585`

## Config

N=16 (N^2=256 pixel dimension), T=6,
n_rx=4 (complex data rows T*n_rx=24;
realified rows 48), k_b=2*pi (f=1.0), Family-1
90-degree arc of radius 1.6 from phi=-45.0
to 45.0 deg, theta = atan2(-p_y,-p_x); rx offsets
[[-0.06, 0.0], [0.06, 0.0], [0.0, -0.06], [0.0, 0.06]], tx offset [0.0, 0.0]; two-blob chi0 (amp
0.3/0.5, sigma
0.09/0.07, centres
[-0.15, -0.12]/[0.18, 0.14]); smooth p=24 unit-column
RBF basis as scene metadata; alpha=1.0.  Ranks use the
family-2 machine rule `tol(M)=max(M.shape)*eps*sigma_1(M)`.  Whitening is
identity (`whiten_realify(A,B,None)`).

## Claim-status table

| claim | status | executed comparison | key numbers |
| --- | --- | --- | --- |
| full-wave composition identity | PASS | stack rel Fro vs hh.build_AB A_c < 1e-12 | rel Fro 0 (per-pose values 0.0, identical float path) |
| Born composition identity | PASS | stack rel Fro vs hh.born_forward < 1e-12 (born_forward available) | rel Fro 7.05442e-17 |
| Range(A_c) subset Range(G_S_stack) | PASS | stable projector/residual norms < 1e-8, rank_A <= rank_GS | proj diff 3.5260e-15, residual 2.3354e-15, ranks 24 <= 24; literal arccos angle 2.1073e-08 rad is the sqrt(2*eps) roundoff floor (literal rad gate False) |
| V_GS vs V_A right singular subspaces are distinct | observed | at least one principal angle > 1e-6 deg (recorded, not forced) | boolean True; max 8.979380e+01 deg; median 2.970670e+01 deg; row proj diff 0.999994 |
| K_eff confounded-eigenvector overlaps (V_A vs V_GS) | recorded | three smallest-eigenvalue u_conf of K_eff(alpha=1), max squared overlaps and top-5 lists | max squared overlaps with V_GS: 2.7206e-03, 3.2680e-03, 3.1695e-03; with V_A: 3.1354e-09, 2.7018e-09, 5.7556e-08 |
| mode-count contrast (SOM vs map information) | recorded | sigma > rel*sigma_max counts for G_S and A_c; eigenvalue > rel*lambda_max for K_IS | G_S=[13, 8], A_c=[23, 12], A_R=[46, 23], K_IS=[23, 12] at thresholds [1e-06, 0.001] |

## Checks / gates

| id | claim | observed | status |
| --- | --- | --- | --- |
| fullwave_composition | full-wave composition A_comp = G_S lu_solve(M, diag(E_tot)) reproduces A_c from build_AB | 0 | PASS |
| born_composition | Born composition A_born_comp = G_S diag(E_inc) reproduces hh.born_forward A_born | 7.05442e-17 | PASS |
| range_containment | Range(A_c) subset Range(G_S_stack): stable projector and residual containment checks < 1e-8 with rank_A <= rank_GS | proj diff 3.5260e-15; residual 2.3354e-15; ranks 24 <= 24; literal arccos angle 2.1073e-08 rad (roundoff floor, literal gate False) | PASS |
| V_not_identical | V_GS and V_A right singular (pixel/current-mode) subspaces are NOT identical | 8.979380e+01 | boolean=True (recorded, not a forced pass) |
| K_eff_overlaps | squared overlaps of the 3 smallest-eigenvalue K_eff eigenvectors with V_A vs V_GS (recorded) | see overlap table | recorded |
| K_eff_formula_consistency | pixel-space K_eff formula matches family4.K_eff(A_R,B_R,alpha) | 2.86628e-16 | PASS |

The V-subspace distinction is a recorded boolean
(`V_not_identical`), not a forced pass: boolean =
**True** (max principal angle
8.979380e+01 deg).

## 1. Composition identity (single unified construction)

One `hh.build_operators` call supplies M, G_S_list, E_inc_list and
E_tot_list; `lu,piv` are computed once from the returned M (the builder's
public dict does not itself return lu/piv; this mirrors build_AB's internal
factorisation).  For each pose t:

* full wave: `A_comp_t = G_S_list[t] @ lu_solve((lu,piv), diag(E_tot_list[t]))`;
* Born: `A_born_t = G_S_list[t] @ diag(E_inc_list[t])`.

Per-pose relative Frobenius errors (denominator = the corresponding block of
`hh.build_AB` / `hh.born_forward`):

| pose t | full-wave comp rel Fro | Born comp rel Fro | full-vs-Born A rel Fro |
| --- | --- | --- | --- |
| 0 | 0 | 7.0869e-17 | 0.175995 |
| 1 | 0 | 7.03537e-17 | 0.184996 |
| 2 | 0 | 7.04068e-17 | 0.196097 |
| 3 | 0 | 7.04068e-17 | 0.20444 |
| 4 | 0 | 7.03537e-17 | 0.208095 |
| 5 | 0 | 7.0869e-17 | 0.208223 |

Whole-stack numbers:

* full-wave composition stack rel Fro =
  0 (gate < 1e-12:
  **True**).
* Born composition stack rel Fro (born_forward available =
  True) = 7.05442e-17
  (gate < 1e-12: **True**).
* ||A_born||_F / ||A_full||_F = 1.00403
  (recorded Born-vs-full-wave scale contrast).
* machine epsilon reference = 2.220e-16.

## 2. Range containment

Column-space containment Range(A_c) subset Range(G_S_stack) is a statement
about the 24-dimensional data-space ranges, measured on the LEFT singular
subspaces U_A and U_GS.  (The right-subspace V_GS^H V_A object requested for
mode distinction is reported in section 3; it is not by itself a containment
measure.)

* rank(G_S_stack) = 24, rank(A_c) = 24
  (rank_A <= rank_GS = **True**).
* direct projector difference ||P_A - P_GS||_2 =
  3.5260e-15; containment residual
  ||(I - P_GS) U_A||_2 = 2.3354e-15
  (stable gate < 1e-8: **True**).
* recorded arccos-based max canonical angle U_A -> U_GS =
  2.1073e-08 rad.  This value sits at the
  arccos(1-eps) roundoff floor ~sqrt(2*eps) = 2.11e-8 rad even for exactly
  identical full-rank subspaces; literal gate
  (< 1e-8 rad) = **False**.  Because
  both ranks equal the 24-dimensional data dimension, Range(A_c) =
  Range(G_S_stack) = C^24 at rank level, and the direct O(eps) projector and
  residual norms are the meaningful machine-precision containment checks.
* orthonormality errors: U_GS 2.516e-15,
  U_A 3.163e-15 (both O(eps), so the
  cross-product singular-value step, not the bases, produces the angle
  roundoff floor).

## 3. Current-space vs map-tangent mode distinction

Right singular vectors: V_GS (columns of the thin SVD of G_S_stack) are the
current-space/SOM pixel modes ordered by sensing strength; V_A (columns of
the thin SVD of A_c) are the map-tangent data-informative pixel modes.  Both
are orthonormal sets in the same N^2 pixel space.

Right-subspace principal angles V_GS vs V_A (24 angles):
count < 1e-8 deg = 0, median =
2.970670e+01 deg, max = 8.979380e+01 deg,
min = 4.660907e-01 deg.

Row-space projector difference ||P(V_A) - P(V_GS)||_2 =
0.999994, i.e. the two pixel-mode subspaces
are nearly orthogonal, in sharp contrast to the column-range containment of
section 2.

Top-10 largest principal angles (deg): 8.979380e+01,
8.908932e+01, 8.683913e+01, 8.507973e+01, 8.294558e+01, 7.882666e+01, 7.368756e+01, 7.066283e+01, 6.226088e+01, 5.613217e+01

K_eff = A_R^T P A_R (P = I - B_R(B_R^T B_R + alpha I)^-1B_R^T, alpha =
1.0).  Three smallest-eigenvalue (most confounded)
eigenvectors u_conf (ascending eigenvalues
-2.0769e-19, -1.5875e-19, -1.5326e-19);
max squared overlaps and total masses:

| u_conf | eig (asc) | max ovl V_GS | max ovl V_A | mass V_GS | mass V_A |
| --- | --- | --- | --- | --- | --- |
| 0 | -2.0769e-19 | 2.7206e-03 | 3.1354e-09 | 5.7891e-03 | 3.1415e-09 |
| 1 | -1.5875e-19 | 3.2680e-03 | 2.7018e-09 | 1.7238e-02 | 2.7405e-09 |
| 2 | -1.5326e-19 | 3.1695e-03 | 5.7556e-08 | 1.2240e-02 | 5.7575e-08 |

Top-5 squared overlaps with V_GS (u_conf index, V_GS column index, value):

| u_conf | V_GS col | squared overlap |
| --- | --- | --- |
| 1 | 18 | 3.2680e-03 |
| 2 | 23 | 3.1695e-03 |
| 0 | 23 | 2.7206e-03 |
| 1 | 2 | 1.5532e-03 |
| 0 | 19 | 1.5444e-03 |

Top-5 squared overlaps with V_A (u_conf index, V_A column index, value):

| u_conf | V_A col | squared overlap |
| --- | --- | --- |
| 2 | 23 | 5.7556e-08 |
| 0 | 23 | 3.1354e-09 |
| 1 | 23 | 2.7018e-09 |
| 1 | 22 | 3.7060e-11 |
| 2 | 21 | 1.4473e-11 |

Pull-through of each u_conf through the representative single-pose current
map c_t = lu_solve(M, diag(E_tot_t)) @ u_conf at pose t =
2: cosine
with the top-5 V_GS current modes (columns ordered by descending sigma_GS):

| u_conf | cosines with top-5 V_GS modes | V_GS top-5 indices |
| --- | --- | --- |
| 0 | 3.582e-05,8.873e-05,1.613e-04,2.951e-04,4.078e-04 | [0, 1, 2, 3, 4] |
| 1 | 7.618e-05,2.395e-04,4.267e-04,1.053e-03,1.583e-03 | [0, 1, 2, 3, 4] |
| 2 | 9.069e-05,2.565e-04,4.858e-04,1.081e-03,1.603e-03 | [0, 1, 2, 3, 4] |

## 4. Mode-count contrast

At relative thresholds 1e-06 and 1e-03 (singular values relative to
sigma_max for G_S and A_c; eigenvalues relative to lambda_max for K_IS):

| quantity | count > 1e-06*max | count > 1e-03*max | max scale |
| --- | --- | --- | --- |
| sigma(G_S_stack) | 13 | 8 | 5.5864e-01 |
| sigma(A_c) | 23 | 12 | 2.6668e-02 |
| sigma(A_R) | 46 | 23 | 2.7082e-02 |
| eig(K_IS = A_R^T A_R) | 23 | 12 | 7.3343e-04 |

Threshold caveat: G_S/A_c counts use sigma > threshold*sigma_max; K_IS counts use eigenvalue > threshold*lambda_max.  Because K_IS = A_R^T A_R, its eigenvalues scale like sigma_A_R^2, so the same literal relative threshold on eigenvalues is NOT commensurable with the relative singular-value thresholds.  A_R sigma counts are recorded as the commensurable companion.  The three counts are therefore
recorded side by side, and the "SOM mode count" (sigma(G_S)) is not silently
identified with the "map-information mode count" (sigma(A_c), eig(K_IS));
each is reported with its own scale and rule.

## Figures

* [figures/family8_composition_identity.png](figures/family8_composition_identity.png) -
  per-pose full-wave and Born composition relative Frobenius errors (log
  scale) with machine-epsilon and 1e-12 gate references.
* [figures/family8_mode_distinction.png](figures/family8_mode_distinction.png) -
  (a) principal-angle spectrum between right singular subspaces V_GS and V_A;
  (b) grouped max squared overlaps of the three most-confounded K_eff
  eigenvectors with V_A vs V_GS.

## Cannot establish

* These are finite-dimensional N=16, T=6, f=1.0 numbers only; the observed
  containment, mode counts, angles, and overlaps are not proved to transfer
  to other N, pose sets, frequencies, or contrast levels.
* The distinction is demonstrated numerically, not proved as a continuum
  theorem; no SOM algorithm quality claim is made and no universal bound on
  the V_GS/V_A separation or on K_eff confoundedness is asserted.
* `K_eff` confounded directions are the smallest-eigenvalue eigenvectors of
  one finite-prior (alpha=1) pixel kernel.  Their overlaps with V_GS/V_A and
  the single representative pull-through are recorded observations, not
  certification that these directions are intrinsically non-estimable.
* Range containment on the tested operators holds up to the recorded
  numerical precision; this does not imply that the current-space sensing
  map and the full-wave map Jacobian have identical singular structure,
  mode counts, or conditioning (they do not, as section 3 shows).
