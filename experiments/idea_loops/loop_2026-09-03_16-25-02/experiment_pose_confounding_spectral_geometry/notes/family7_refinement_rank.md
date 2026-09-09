# Family 7: refinement / stable-transversality diagnostics and rank events

Date: 2026-09-03T11:38:33+00:00 UTC.  Experiment:
`experiment_pose_confounding_spectral_geometry`.

Family 7 is a NEW finite-grid diagnostic family.  It reuses (and does not
modify) `helmholtz.py`, `family1_pilot.py`, `family2_algebraic_spine.py`, and
`family4_frequency_trajectory.py`.  Its scope is the discrete
whitened/realified identity-noise model: for every resolution N the same
continuous two-blob chi0 is sampled on that N-grid and the same smooth p=24
unit-column basis is rebuilt on that grid.

## Exact command and runtime

```bash
.venv/bin/python src/family7_refinement_rank.py
```

Wall runtime: 18.54 s.  Platform:
macOS-26.6.2-arm64-arm-64bit, Python 3.12.13,
numpy 2.5.2, scipy 1.18.1,
matplotlib 3.11.1.

Source SHA-256 (this script):
`6576bf5249fab8aaf0deca185b56af8324c76be4eeab499b8455942ffaa29ab0`

Reused-module SHA-256:
`family1_pilot.py` `762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7`, `family2_algebraic_spine.py` `173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235`, `family4_frequency_trajectory.py` `ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242`, `helmholtz.py` `ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585`

## Config

Scenario matches Family 4/6: T=6, n_rx=4,
rx_offsets=[[-0.06, 0.0], [0.06, 0.0], [0.0, -0.06], [0.0, 0.06]], tx_offset=[0.0, 0.0], 90-degree arc
at radius 1.6 from phi=-45 to 45 deg, two-blob chi0 amp 0.3/0.5 sigma
0.09/0.07 centres (-0.15,-0.12)/(0.18,0.14), smooth p=24 basis with
x_centers_n=4, y_centers_n=6, x_span=y_span=[-0.3,0.3], sigma_b=0.16,
unit_columns=True.  Ranks use the stated family-2 machine rule
`tol(M)=max(M.shape)*eps*sigma_1(M)`.  Noise whitening is identity
(`whiten_realify(A,B,None)`).

Part 1 grid: N in [16, 24, 32, 40] always, N=48 guarded
(time guard 180.0 s);
N48 status = `computed`.  Part 2: f in
linspace(0.6,
2.6,21)
at N=16.  Part 3: s in [0.0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0] at N=16,
f=1.0, chi=s*chi0.

## Claim-status table

| claim | status | executed comparison | key numbers |
| --- | --- | --- | --- |
| theta_min_deg trend across tested N | observed | finite differences on N in [16, 24, 32, 40, 48] | decreasing_observed; diffs [-4.0806386778362595e-06, -1.406644138250029e-06, -6.38572840663701e-07, -3.479317355827732e-07]; span 6.473787e-06 deg |
| near-confounded direction count (cos2>1-1e-8) | observed | per-N counts on the tested grids | [4, 4, 4, 4, 4]; constant over grid = True |
| rank identity r(K_IS)-r(K_SLAM) == r(A)+r(B)-r([A,B]) | executed | machine-rank identity for every computed N row | N16:True, N24:True, N32:True, N40:True, N48:True |
| frequency rank events | observed/flagged | rank(B) or rank([A_s,B]) change vs previous f | 0 event(s): [] |
| frequency near crossings | observed/flagged | smallest K_eff gap < 1e-6*max(1,lambda_max) | 21 flag(s) at f = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6] (lowest ascending pair; details in Part 2 table) |
| s=0 B exactly zero | gate | ||B_R||_F/||A_R||_F < 1e-12 | ratio 0.0; pass=True |
| s=0 K_SLAM == K_IS | gate | rel Fro < 1e-12 | rel Fro 0.0; pass=True |
| contrast rank / near-degenerate transitions | observed | rank changes and rho_min<1e-6 as s increases | rank events [{'s': 0.01, 'changed': ['rank_B', 'rank_AB'], 'prev_rank_B': 0, 'new_rank_B': 18, 'prev_rank_AB': 24, 'new_rank_AB': 42}]; near-degenerate onset [{'s': 0.01, 'rho_min': 8.339151591485461e-11, 'theta_min_deg': 0.0005232190726495859, 'rule': 'rho_min < 1e-6 (dimensionless retention rho)'}] |

## Part 1: resolution / stable-transversality table

| N | r_A | r_B | r_AB | r_KIS | r_KSL | rank identity | theta_min deg | confusable mass | retained mass | log_volume | count rho<1e-6 | rho_min | state res | sigma_min(M)/||M|| | build s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 16 | 24 | 18 | 42 | 24 | 24 | True | 4.253070e-04 | 14.5928 | 9.40716 | -167.7461 | 5 | 5.5101e-11 | 9.34e-16 | 0.718892 | 0.057 |
| 24 | 24 | 18 | 42 | 24 | 24 | True | 4.212263e-04 | 14.5809 | 9.41915 | -167.8156 | 5 | 5.4049e-11 | 1.28e-15 | 0.718670 | 0.276 |
| 32 | 24 | 18 | 42 | 24 | 24 | True | 4.198197e-04 | 14.5768 | 9.42324 | -167.8387 | 5 | 5.3688e-11 | 1.92e-15 | 0.718583 | 1.076 |
| 40 | 24 | 18 | 42 | 24 | 24 | True | 4.191811e-04 | 14.5749 | 9.42512 | -167.8492 | 5 | 5.3525e-11 | 2.21e-15 | 0.718541 | 3.094 |
| 48 | 24 | 18 | 42 | 24 | 24 | True | 4.188332e-04 | 14.5739 | 9.42614 | -167.8548 | 5 | 5.3436e-11 | 2.99e-15 | 0.718518 | 7.922 |

`state res` is the largest relative M J = chi*E_inc residual returned by
`hh.build_operators`; `sigma_min(M)/||M||` is the normalised M conditioning
indicator.  Interpretation (finite-grid observations only):

* theta_min deg: [0.00042530697982207, 0.00042122634114423376, 0.00041981969700598373, 0.00041918112416532003, 0.00041883319242973726]
* finite differences between successive computed N:
  [-4.0806386778362595e-06, -1.406644138250029e-06, -6.38572840663701e-07, -3.479317355827732e-07]
* direction: **decreasing_observed** (span
  6.473787e-06 deg, max absolute step
  4.080639e-06 deg, last step
  -3.479317e-07 deg)
* near-confounded counts (cos2 > 1-1e-8):
  [4, 4, 4, 4, 4]; constant over the tested grid =
  **True**

Top-10 cos2 values per N are in the JSON (`top10_cos2`), together with full
cos2_desc, rho spectra, rank tolerance, and build times.

## Part 2: frequency rank-event sweep (N=16, smooth p=24, T=6)

| f | rank(B) | sigma_min(B) | theta_min deg | rho_min | count rho<1e-6 | retained mass | confusable mass | smallest K_eff gap | gap pair | near-crossing | rank event |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.600 | 18 | 3.6285e-05 | 5.615469e-05 | 9.6056e-13 | 4 | 10.64199 | 13.35801 | 8.933e-17 | [0,1] | True | False |
| 0.700 | 18 | 2.6648e-04 | 9.896780e-05 | 2.9836e-12 | 5 | 9.96110 | 14.03890 | 3.302e-16 | [0,1] | True | False |
| 0.800 | 18 | 3.9583e-04 | 1.470500e-04 | 6.5870e-12 | 6 | 9.92894 | 14.07106 | 7.953e-16 | [0,1] | True | False |
| 0.900 | 18 | 4.6635e-04 | 2.101983e-04 | 1.3459e-11 | 6 | 9.67312 | 14.32688 | 1.953e-15 | [0,1] | True | False |
| 1.000 | 18 | 6.7489e-04 | 4.253070e-04 | 5.5101e-11 | 5 | 9.40716 | 14.59284 | 6.935e-15 | [0,1] | True | False |
| 1.100 | 18 | 3.7242e-04 | 5.632034e-04 | 9.6624e-11 | 6 | 8.74966 | 15.25034 | 1.219e-14 | [1,2] | True | False |
| 1.200 | 18 | 2.0375e-05 | 6.186790e-04 | 1.1660e-10 | 5 | 8.57338 | 15.42662 | 1.512e-14 | [0,1] | True | False |
| 1.300 | 18 | 4.2010e-04 | 6.727379e-04 | 1.3786e-10 | 6 | 8.37101 | 15.62899 | 5.025e-14 | [0,1] | True | False |
| 1.400 | 18 | 5.7938e-04 | 6.051455e-04 | 1.1155e-10 | 5 | 8.58151 | 15.41849 | 1.499e-13 | [0,1] | True | False |
| 1.500 | 18 | 7.9719e-04 | 8.398868e-04 | 2.1488e-10 | 5 | 8.48147 | 15.51853 | 5.292e-14 | [0,1] | True | False |
| 1.600 | 18 | 7.9939e-04 | 3.757060e-04 | 4.2998e-11 | 5 | 8.12943 | 15.87057 | 1.385e-13 | [1,2] | True | False |
| 1.700 | 18 | 3.7235e-04 | 1.499307e-03 | 6.8476e-10 | 6 | 7.77958 | 16.22042 | 1.004e-13 | [1,2] | True | False |
| 1.800 | 18 | 1.5363e-04 | 1.323244e-03 | 5.3338e-10 | 5 | 7.64326 | 16.35674 | 7.037e-14 | [0,1] | True | False |
| 1.900 | 18 | 3.9888e-04 | 8.562392e-04 | 2.2333e-10 | 6 | 7.56932 | 16.43068 | 2.468e-13 | [0,1] | True | False |
| 2.000 | 18 | 2.0924e-04 | 1.127502e-03 | 3.8725e-10 | 6 | 7.10142 | 16.89858 | 8.747e-13 | [0,1] | True | False |
| 2.100 | 18 | 2.1953e-04 | 1.630713e-03 | 8.1005e-10 | 6 | 6.99177 | 17.00823 | 1.611e-13 | [0,1] | True | False |
| 2.200 | 18 | 7.3531e-05 | 1.434864e-03 | 6.2716e-10 | 6 | 6.84730 | 17.15270 | 5.631e-13 | [0,1] | True | False |
| 2.300 | 18 | 3.1757e-05 | 3.004975e-03 | 2.7507e-09 | 5 | 6.63141 | 17.36859 | 1.777e-12 | [0,1] | True | False |
| 2.400 | 18 | 2.1663e-05 | 3.460966e-03 | 3.6488e-09 | 6 | 6.99111 | 17.00889 | 4.041e-13 | [0,1] | True | False |
| 2.500 | 18 | 6.8498e-05 | 3.780128e-03 | 4.3528e-09 | 6 | 6.71374 | 17.28626 | 1.897e-13 | [0,1] | True | False |
| 2.600 | 18 | 2.4283e-05 | 4.007016e-03 | 4.8910e-09 | 6 | 7.07854 | 16.92146 | 4.697e-13 | [0,1] | True | False |

Rank events: **0** (details in JSON and below).
Near crossings (smallest adjacent K_eff(alpha=1) gap < 1e-6*max(1,lambda_max)):
**21**.

Literal-flag caveat (recorded, not forced): every near-crossing flag sits at
one of the lowest ascending pairs (recorded pairs are [0,1] or [1,2] over the
sweep), where numerical near-zero K_eff eigenvalues are separated by gaps up
to ~1.8e-12, many orders below the 1e-6 literal threshold.  The requested
inequality is therefore satisfied at every sweep point, and the flag is
dominated by the repeated numerical near-null tail of K_eff.  It is recorded
as the literal criterion and is **not** interpreted as a certified eigenvalue
crossing or a transversality failure.

Recorded rank events:

| f | changed | prev rank(B) | new rank(B) | prev rank([A,B]) | new rank([A,B]) |
| --- | --- | --- | --- | --- | --- |

Recorded near crossings:

| f | smallest gap | threshold | gap pair |
| --- | --- | --- | --- |
| 0.600 | 8.933e-17 | 1.000e-06 | [0,1] |
| 0.700 | 3.302e-16 | 1.000e-06 | [0,1] |
| 0.800 | 7.953e-16 | 1.000e-06 | [0,1] |
| 0.900 | 1.953e-15 | 1.000e-06 | [0,1] |
| 1.000 | 6.935e-15 | 1.000e-06 | [0,1] |
| 1.100 | 1.219e-14 | 1.000e-06 | [1,2] |
| 1.200 | 1.512e-14 | 1.000e-06 | [0,1] |
| 1.300 | 5.025e-14 | 1.000e-06 | [0,1] |
| 1.400 | 1.499e-13 | 1.000e-06 | [0,1] |
| 1.500 | 5.292e-14 | 1.000e-06 | [0,1] |
| 1.600 | 1.385e-13 | 1.000e-06 | [1,2] |
| 1.700 | 1.004e-13 | 1.000e-06 | [1,2] |
| 1.800 | 7.037e-14 | 1.000e-06 | [0,1] |
| 1.900 | 2.468e-13 | 1.000e-06 | [0,1] |
| 2.000 | 8.747e-13 | 1.000e-06 | [0,1] |
| 2.100 | 1.611e-13 | 1.000e-06 | [0,1] |
| 2.200 | 5.631e-13 | 1.000e-06 | [0,1] |
| 2.300 | 1.777e-12 | 1.000e-06 | [0,1] |
| 2.400 | 4.041e-13 | 1.000e-06 | [0,1] |
| 2.500 | 1.897e-13 | 1.000e-06 | [0,1] |
| 2.600 | 4.697e-13 | 1.000e-06 | [0,1] |

## Part 3: contrast sweep Born -> full wave (N=16, f=1.0)

| s | B/A Fro ratio | rank(A) | rank(B) | rank([A,B]) | theta_min deg | retained mass | rho_min | Born-vs-full A rel Fro | KSL vs KIS rel Fro | rank event |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.000e+00 | 24 | 0 | 24 | 9.000000e+01 | 24.00000 | 1.000e+00 | 7.054e-17 | 0.000e+00 | False |
| 0.01 | 2.720e-02 | 24 | 18 | 42 | 5.232191e-04 | 9.15633 | 8.339e-11 | 1.894e-03 | 9.967e-01 | True |
| 0.05 | 1.363e-01 | 24 | 18 | 42 | 5.206346e-04 | 9.18256 | 8.257e-11 | 9.488e-03 | 9.967e-01 | False |
| 0.1 | 2.735e-01 | 24 | 18 | 42 | 5.175523e-04 | 9.22456 | 8.159e-11 | 1.902e-02 | 9.966e-01 | False |
| 0.2 | 5.503e-01 | 24 | 18 | 42 | 5.119802e-04 | 9.30990 | 7.985e-11 | 3.823e-02 | 9.966e-01 | False |
| 0.4 | 1.114e+00 | 24 | 18 | 42 | 5.033950e-04 | 9.37212 | 7.719e-11 | 7.718e-02 | 9.966e-01 | False |
| 0.7 | 1.983e+00 | 24 | 18 | 42 | 4.860312e-04 | 9.38510 | 7.196e-11 | 1.368e-01 | 9.967e-01 | False |
| 1 | 2.879e+00 | 24 | 18 | 42 | 4.253070e-04 | 9.40716 | 5.510e-11 | 1.975e-01 | 9.970e-01 | False |

s=0 gates:

* B/A Frobenius ratio = 0.0
  (gate < 1e-12: **True**).
* K_SLAM vs K_IS rel Fro =
  0.0
  (gate < 1e-12: **True**).

B/A ratio rule: ||B_R||_F/||A_R||_F on the realified blocks (same as the
complex ratio under whiten_realify).  Born discrepancy rule:
||A_full-A_born||_F/||A_full||_F with A_born from `hh.born_forward` (available
= True).  Contrast rank/near-degenerate transitions:
1 rank event(s); details in JSON.

Near-degenerate transitions (rho_min < 1e-6, dimensionless retention rho):
1 recorded onset(s):
| s | rho_min | theta_min deg |
| --- | --- | --- |
| 0.01 | 8.3392e-11 | 5.232191e-04 |

At s=0 the principal-angle convention for a zero Range(B) is declared as
theta_min = 90 deg (the limiting/maximal convention); no physical angle
exists between a nonempty subspace and the trivial subspace.  The first
contrast "rank event" (s=0 -> 0.01, rank(B): 0 -> 18, rank([A,B]): 24 -> 42)
is the expected exit from the exactly-zero B baseline and is recorded as an
observed machine-rank change, not a certified degeneracy transition.

## Figures

* [figures/family7_resolution_convergence.png](figures/family7_resolution_convergence.png) - resolution
  convergence of theta_min_deg, confusable_mass, retained_mass, and count
  rho<1e-6 versus N.
* [figures/family7_frequency_rank_sweep.png](figures/family7_frequency_rank_sweep.png) - frequency sweep of
  sigma_min(B), theta_min_deg, rho_min, and the smallest K_eff gap with rank
  event / near-crossing lines.
* [figures/family7_contrast_sweep.png](figures/family7_contrast_sweep.png) - theta_min_deg and
  retained_mass versus s, log-scale B/A ratio, and Born discrepancy.

## Cannot establish

* These are finite-grid diagnostics only.  They do **not** prove continuum
  transfer, limiting transversality as N -> infinity, or any theorem about
  rank/transversality stability away from the tested resolutions.
* N refinement is diagnostic; it does not certify convergence of principal
  angles, masses, rho spectra, or rank identities to continuum limits.
* The frequency sweep flags rank events and near crossings on a discrete
  21-point f grid with one fixed N=16 discretisation.  Events between grid
  points, and any N-dependence of those events, are not established.
* The contrast sweep is one 8-point scale path at N=16, f=1.0.  Rank
  "events" are machine-rank changes on that path, not certified degeneracies;
  no Born/full-wave regime boundary is claimed.
* No forced pass is applied anywhere in Family 7; all trends and events are
  recorded as observed on the exact executed grid/config.
