# Round-2 synthesis (read-only)

Date: 2026-09-03. Experiment:
`experiment_pose_confounding_spectral_geometry`. This synthesis consolidates
the round-2 family reports into the updated manuscript draft
(`manuscript/main_round2.tex`). It is a read-only synthesis: no experiment was
rerun, and every number below was copied from the cited reports and their
JSON result files.

## New families, roles, and artifacts

| family | role | report | result JSON | figures |
| --- | --- | --- | --- | --- |
| Family 6 | declared per-frequency colored-noise whitening and stacking checks | `notes/family6_report.md` | `results/family6_colored_noise.json` | `figures/family6_noise_whitening_sanity.png`, `family6_snr_sweep.png`, `family6_correlation_effect.png`, `family6_retention_spectra_colored.png` |
| Family 4c | same-standoff, per-pose-energy-normalized trajectory control | `notes/family4c_trajectory_control.md` | `results/family4c_trajectory_control.json` | `figures/family4c_trajectory_control.png` |
| Family 7 | resolution / rank-event / stable-transversality diagnostics on finite grids | `notes/family7_refinement_rank.md` | `results/family7_refinement_rank.json` | `figures/family7_resolution_convergence.png`, `family7_frequency_rank_sweep.png`, `family7_contrast_sweep.png` |
| Family 8 | current-space `G_S`/SOM vs map-tangent `A`/`K_eff` distinction | `notes/family8_som_distinction.md` | `results/family8_som_distinction.json` | `figures/family8_composition_identity.png`, `family8_mode_distinction.png` |
| Family 9 | seed replications and pose/receiver/T/basis ablations | `notes/family9_replications_ablation.md` | `results/family9_replications_ablation.json` | `figures/family9_seed_replications.png`, `family9_pose_dof_ablation.png`, `family9_rx_T_basis_ablation.png` |
| Family 5b (precursor) | loose structural Lipschitz bound, cited only as corrected | `notes/family5b_lipschitz_bound.md` | `results/family5b_lipschitz_bound.json` | `figures/family5b_lipschitz_bound.png` |
| Family 5b2 | alpha-aware tight structural operator-Lipschitz bound | `notes/family5b2_alpha_tight.md` | `results/family5b2_alpha_tight.json` | `figures/family5b2_alpha_tight.png` |

All round-2 families were executed with `.venv/bin/python`; each round-2
result JSON embeds its exact command, generated-UTC timestamp, configuration,
and script/reused-module SHA-256 digests. Wall runtimes: Family 4c 0.147 s,
Family 6 1.15 s, Family 7 18.54 s, Family 8 0.45 s, Family 9 1.35 s,
Family 5b 40.87 s, Family 5b2 41.19 s.

## Key numbers copied from the reports

### Family 6 (colored-noise / per-frequency SNR whitening)

- Analytic whitening `rel Fro ||W_f C_f W_f^H - I||/||I||`: `7.370e-16`,
  `5.550e-16`, `5.769e-16` (f = 1.0 / 1.4 / 1.8); max ~ `7e-16`.
- Monte-Carlo whitening (seed 2718): rel Fro ~ `0.107`, recorded as a
  statistical `O(sqrt(2/n))` check, not machine whiteness.
- Algebraic spine under colored noise: rank residual `0`
  (rA=24, rB=18, rAB=42); no-prior retention identity max abs diff
  `1.776e-15`; count rho<1 = 18; rho_min ~ `1.973e-10`.
- PSD prefix increments at snr2=100: K_SLAM `1.220e-7`, K_eff `7.424e-8`
  (PASS).
- Duplicate invariance with same C_f (c=2, alpha=1): no-prior rel Fro K_IS
  `1.176e-15`, K_SLAM `8.837e-14`, max |rho diff| `6.216e-13`; fixed finite
  prior observed: max |rho diff| `0.3398`, rel Fro K_eff `2.546`.
- SNR movements at snr2=100 (u0/u1/u2): `0.127704` / `0.235057` / `0.568961`
  (non-monotone peak near the balanced-SNR region).
- SNR movements at snr2=1e-2 (u0/u1/u2): `5.816e-5` / `9.563e-5` /
  `5.533e-4`; at snr2=1e10: `0.010205` / `0.009420` / `0.026538` (finite
  positive high-SNR plateau).
- Receiver-correlation sweep (ell = 0.001..0.5): recorded, no forced pass;
  theta_min ranges ~`4.253e-4` to `1.363e-3` deg; identity-distance column is
  mixed-scale and is not used for geometric statements.

### Family 4c (same-standoff, per-pose-energy-normalized trajectory control)

- Normalization verification: max `|1 - ||A_block||_F| = 2.220e-16` across
  all five trajectories (gate `<= 1e-12`; PASS).
- Raw vs normalized orderings do not flip: retained desc and confusable asc
  are both `arc90 > straight > arc180 > circle360 > arc270`.
- Four-core hypothesis remains false under normalization
  (`hypothesis_holds_for_four_core = false`, no forced pass):
  observed arc90 > straight > arc180 > circle360.
- Arc270 is the lowest retained-mass trajectory in the tested five-path set.
- Normalized retained mass: arc90 `9.408661`, straight `8.269820`,
  arc180 `7.997790`, circle360 `7.701223`, arc270 `7.535343`.
- Per-pose A-block normalization is a declared gain control; it is not a
  physical noise covariance, and equal standoff does not isolate angular
  coverage alone (continuous/polyline lengths still change).

### Family 7 (resolution, rank events, stable-transversality diagnostics)

- Ranks over N = 16/24/32/40/48: rA=24, rB=18, rAB=42, rK=24; rank identity
  true for every computed row.
- theta_min deg: monotone decreasing `4.253070e-4` (N=16) to
  `4.188332e-4` (N=48); span `6.474e-6` deg.
- Near-confounded direction count: `[4, 4, 4, 4, 4]` (constant on tested
  grid).
- Frequency sweep (21 points, f = 0.6..2.6): rank events `0`; 21 literal
  near-crossing flags, all on the lowest ascending numerical near-null K_eff
  pairs ([0,1] or [1,2]); not certified crossings or transversality failures.
- s=0 Born gates pass exactly: B/A Fro ratio `0.0`; rel Fro
  `K_SLAM - K_IS` = `0.0` (both `< 1e-12`).
- Born/full discrepancy grows `1.894e-3` (s=0.01) to `0.1975` (s=1).
- Contrast rank transition s=0 -> 0.01 is the expected machine-rank exit
  (rB 0 -> 18, rAB 24 -> 42), recorded as observed.
- These are finite-grid diagnostics only; no continuum/transversality proof.

### Family 8 (current-space `G_S`/SOM vs map-tangent `A`/`K_eff`)

- Full-wave composition stack rel Fro: `0` (identical float path; PASS).
- Born composition stack rel Fro: `7.05442e-17` (PASS).
- Range(A_c) subset Range(G_S_stack): projector difference `3.5260e-15`,
  containment residual `2.3354e-15`, ranks 24 <= 24 (PASS at machine
  precision). The literal arccos angle `2.1073e-8` rad is the sqrt(2 eps)
  roundoff floor for identical full-rank subspaces (literal rad gate false).
- Right singular subspaces V_GS vs V_A are strongly distinct: 0/24 principal
  angles below `1e-8` deg; median `29.70670` deg; max `89.79380` deg; row
  projector diff `0.999994`.
- K_eff confounded directions (3 smallest-eigenvalue eigenvectors, alpha=1):
  max squared overlaps with V_GS `2.7206e-3` / `3.2680e-3` / `3.1695e-3`;
  with V_A `3.1354e-9` / `2.7018e-9` / `5.7556e-8` (nearly orthogonal to
  both mode families on this scenario).
- Mode counts at rel thresholds 1e-6/1e-3: sigma(G_S) `[13, 8]`;
  sigma(A_c) `[23, 12]`; sigma(A_R) `[46, 23]`; eig(K_IS) `[23, 12]`;
  threshold semantics recorded separately.

### Family 9 (seed replications and component ablations)

- Pass rates over 7 scenes: rank identity `6/7` (ring scene residual 1,
  machine-rank boundary); retention identity `7/7`; kernel identity F2 gate
  `7/7` with identity support `6/7` and nullity match `6/7`; frequency
  diversity `7/7` (3/3 directions moved > 1e-4 on every scene).
- Ring scene boundary: smallest K_SLAM eigenvalue below the family-2 eig rank
  tolerance while C stays full rank at its own SVD tolerance; preserved as a
  machine-rank boundary event, not evidence against the exact identity.
- Pose-DOF retained mass: known pose `24.000`; x-only `18.020`; y-only
  `18.029`; theta-only `21.362`; xy `12.573`; full `9.407`.
- Receiver-count retained mass (n_rx = 1/2/4/8): `0 / 6.000 / 9.407 / 9.399`;
  mixed-rA rows, not a pure receiver-count effect.
- Pose-count retained mass (T = 3/6/12): `15.000 / 9.407 / 6.579` (monotone
  on the tested grid).
- Basis retained mass (p = 12/24/36): `0.157 / 9.407 / 18.750` (rA grows for
  fixed rB).
- Cross-check vs stored Family 4 result: max movement difference `0`
  (< 1e-12).

### Family 5b2 (alpha-aware tight operator-Lipschitz bound)

- Tight pointwise structural constant `L_struct,tight(X0) = 1.124e-2` versus
  the loose Family 5b value `1.408e7` (cited only as corrected); tightness
  ratio `1.252285e9`.
- Exact alpha-aware factor: `||C^-1||_2 = 1/(sigma_min(B)^2 + alpha)` for
  C = B^T B + alpha I; alpha = 1 keeps the affine C_aff invertible for all
  eps.
- Affine tangent certificate: 2000 unit directions (seed 7171), all
  eps in {1e-3, 3e-3, 1e-2, 3e-2, 1e-1}; every Frobenius and eigenvalue-drop
  check true, zero violations.
- Full nonlinear model: 500 unit directions (seed 9191), worst Frobenius
  quotient `3.228e-4` below L_struct,tight(X0) `1.124e-2`; no sampled
  violations; sampled ingredient envelope `1.125e-2`.
- Full-model result is a derived structural bound with a sampled envelope,
  NOT a formal nonlinear interval certificate.
- FD Jacobian validation: global rel Fro J_A `5.991e-7`, J_B `5.374e-7`;
  O(h^2) ratios ~4.

## Consolidated claim-status table (as reported)

| family | claim | status | key evidence |
| --- | --- | --- | --- |
| 6 | analytic whitening W C W^H = I | PASS | rel Fro `7.370e-16` / `5.550e-16` / `5.769e-16` |
| 6 | Monte-Carlo whitening | recorded | seed 2718 rel Fro ~0.107, statistical |
| 6 | colored kernel/rank identity | PASS | residual 0 |
| 6 | colored no-prior retention identity | PASS | max diff `1.776e-15`; count rho<1 = 18 |
| 6 | PSD prefix increments (snr2=100) | PASS | K_SLAM `1.220e-7`; K_eff `7.424e-8` |
| 6 | duplicate no-prior / fixed prior | PASS / observed | max rho `6.216e-13` / `0.3398` |
| 6 | SNR sweep trend | PASS (tested grid) | peak at snr2=100; positive plateau at 1e10 |
| 6 | SNR non-monotonicity | observed | movements 0.128/0.235/0.569 peak |
| 6 | receiver-correlation sweep | recorded, no forced pass | ell 0.001..0.5 |
| 4c | per-pose A-block normalization | PASS | max dev `2.220e-16` |
| 4c | raw vs normalized ordering flip | no flip (observed) | same orders in both columns |
| 4c | four-core trajectory hypothesis | FAIL (counterexample preserved) | arc90 > straight > arc180 > circle360 |
| 4c | arc270 lowest retained | observed | `7.5353` normalized |
| 4c | normalization is gain control | recorded | not physical noise |
| 7 | constant ranks and rank identity | PASS (computed rows) | rA=24, rB=18, rAB=42 |
| 7 | theta_min decreasing | observed | `4.253070e-4` -> `4.188332e-4` deg |
| 7 | near-confounded count constant | observed | [4,4,4,4,4] |
| 7 | frequency rank events | observed/flagged | 0 events; 21 literal near-crossing flags |
| 7 | s=0 Born gates | PASS | 0.0 / 0.0 rel Fro |
| 7 | Born/full discrepancy | observed | `1.894e-3` -> `0.1975` |
| 7 | finite-grid scope | recorded | diagnostics only, no continuum proof |
| 8 | full-wave composition identity | PASS | stack rel Fro 0 |
| 8 | Born composition identity | PASS | rel Fro `7.05442e-17` |
| 8 | Range(A_c) subset Range(G_S) | PASS | proj diff `3.5260e-15`; residual `2.3354e-15` |
| 8 | V_GS vs V_A distinct | observed | 0/24 angles < 1e-8 deg; median 29.71 deg; max 89.79 deg |
| 8 | K_eff overlaps | recorded | max `3.2680e-3` (V_GS) vs `5.7556e-8` (V_A) |
| 8 | mode counts differ | recorded | G_S [13,8]; A_c [23,12]; eig(K_IS) [23,12] |
| 9 | rank identity over 7 scenes | PARTIAL | 6/7; ring residual 1 (boundary) |
| 9 | no-prior retention identity | PASS | 7/7 |
| 9 | kernel identity | PASS (gate) / PARTIAL (support) | gate 7/7; support/nullity 6/7 |
| 9 | frequency diversity | PASS | 7/7; 3/3 directions every scene |
| 9 | ring rank-boundary event | preserved boundary | K_SLAM smallest eigenvalue below eig tolerance |
| 9 | pose/receiver/T/basis ablations | recorded | trends only, no universal claims |
| 5b2 | exact alpha-aware ||C^-1|| | PASS (derived) | `1/(sigma_min(B)^2 + alpha)` |
| 5b2 | tight L_struct | PASS (derived) | `1.124e-2` vs loose `1.408e7` (1.252285e9x) |
| 5b2 | affine tangent certificate | PASS (affine model only) | 2000 dirs; eps <= 0.1; zero violations |
| 5b2 | full nonlinear sampled quotient | observed | worst `3.228e-4`; no sampled violation |
| 5b2 | full nonlinear interval certification | NOT MET | derived structural bound + sampled envelope only |

## Cannot establish (round 2)

- Continuum transfer, stable transversality, or convergence of the retention
  spectrum under discretization: Family 7 is a finite-grid diagnostic on
  N = 16..48, not a proof.
- Physical or measured colored-noise realism: the Family 6 covariance is a
  declared synthetic model (`sigma_f^2 (I_T kron R_rx)` with Matérn-1/2
  receiver factor); no estimator or physical-noise claim.
- Machine-precision whiteness from the Monte-Carlo check: only finite-sample
  decorrelation at `O(sqrt(2/n))` is established.
- A universal trajectory-ordering theorem: Family 4c normalization is a gain
  control, same standoff still changes length/budget, and the straight path
  still varies standoff.
- A certified eigenvalue crossing or transversality failure from the Family 7
  near-crossing flags: all 21 flags sit in the lowest numerical near-null
  pairs of K_eff.
- A continuum theorem for the `G_S`/SOM vs map-tangent distinction: Family 8
  is finite-dimensional evidence for one scenario; no pullback/composition
  theorem `A = G_S T_chi` and no SOM quality claim.
- An exact-algebraic failure from the Family 9 ring rank mismatch: it is a
  machine-rank tolerance boundary event.
- Universal component-ablation laws from Family 9 parts B-E: the trends are
  descriptive for the tested grid only.
- A formal interval certificate of the full nonlinear operator-Lipschitz
  constant: Family 5b2 certifies the affine tangent model only, conditional
  on numerically computed FD Jacobians; the full-model statement is a derived
  structural bound with a sampled envelope.
- Global, online-SLAM, production, or novelty acceptance: no continuum/global/
  online claim is made by any round-2 family.

## Compile status

- Command: `tectonic manuscript/main_round2.tex`
- Result: success (exit status 0); after the first successful compile, two
  small manuscript corrections (a recorded-number typo and a reproducibility
  figure-list addition) were made and the file was recompiled successfully.
- PDF: `manuscript/main_round2.pdf` (2.06 MiB; written by tectonic after the
  standard aux/rerun passes).
- Notes: only overfull/underfull TeX layout warnings were reported; no
  LaTeX errors. `manuscript/main.tex` was not modified; the round-2 draft is
  `manuscript/main_round2.tex`.
- This file is a read-only synthesis; no experiment was rerun.
