# Family 6: colored-noise / per-frequency SNR whitening model

Date: 2026-09-03T11:22:33.133908+00:00 UTC (SGT; experiment run 2026-09-03).
Experiment: `experiment_pose_confounding_spectral_geometry`.

Family 6 reuses the corrected Family 1/2/4 construction (helmholtz
`build_AB`, Family-1 scene and arc poses, Family-2 smooth p=24 unit-column
RBF basis and machine-rank rule) and does not modify any existing source or
results.  It replaces the raw identity-noise metric of Family 4 with a
declared per-frequency colored-noise model and reruns the algebraic,
duplicate, and frequency-stacking questions on the colored-whitened blocks.
All numbers are finite-dimensional statements about the discrete N=16
whitened/realified dense linear algebra.

## Exact command and runtime

```bash
.venv/bin/python src/family6_colored_noise.py
```

Wall runtime: 1.15 s.  Platform:
macOS-26.6.2-arm64-arm-64bit, Python 3.12.13,
numpy 2.5.2, scipy 1.18.1,
matplotlib 3.11.1.  Deterministic except for
the fixed-seed Monte Carlo whitening sanity check (seeds 2718, 4242, 5150).

Scenario: N=16, T=6 poses, n_rx=4, m_c=24, m_real=48 per frequency, p=24
smooth basis, two-blob chi0 (amp 0.3/0.5, sigma 0.09/0.07 at
(-0.15,-0.12)/(0.18,0.14)), 90-degree arc at radius 1.6.

## Colored-noise / whitening model (exactly as implemented)

For each frequency f the complex (pre-whitening, pre-realification)
Jacobians are `A_c` (24 x 256) and `B_c` (24 x 18) from
`hh.build_AB(chi0, poses, rx_offsets, tx_offset, N, 2*pi*f)`.

* Noise covariance: `C_f = sigma_f^2 * (I_T kron R_rx)`, where
  `R_rx[i,j] = exp(-|rx_i - rx_j| / ell)` is the declared 2D Matérn-1/2
  receiver covariance (row/column order matches the receiver order used by
  `build_AB`; poses are independent).
* Per-frequency SNR: `sigma_f^2 = ||A_c||_F^2 / (m_c * snr_f)` with linear
  snr (reference snr_f = 100, i.e. 20 dB).
* Whitening: `C_f = V diag(d) V^H` and
  `W_f = V diag(1/sqrt(max(d, tiny))) V^H` with tiny = 1e-30; a 1e-14*I
  diagonal jitter is added only if `cond(C_f) > 1e12` or a non-positive
  eigenvalue appears (never triggered for the ell values in this run).
* Whitened complex blocks: `A_w = W_f A_c`, `B_w = W_f B_c`, then realified
  by `hh.whiten_realify(A_w, B_w, None)` =
  `sqrt(2)*[Re; Im]` (W is already applied, so None only realifies).
* Smooth-basis blocks: `A_s = A_pix_R @ S`; multi-frequency stacks are
  vertical stacks of the real whitened blocks.

Unlike Family 4's raw identity-noise stack, each colored frequency block is
SNR-normalized by construction (its whitened row scale is set by
sigma_f = ||A_c||_F / sqrt(m_c * snr_f)), so check C compares equal-SNR
blocks rather than raw block energies.

## Claim-status table

| check | status | executed comparison | key numbers |
|---|---|---|---|
| A: analytic whitening | **PASS** | rel Fro ||W_f C_f W_f^H - I|| / ||I|| | 7.370e-16 (f=1.0), 5.550e-16 (f=1.4), 5.769e-16 (f=1.8) |
| A: Monte Carlo whitening | **recorded** | sample covariance of whitened-realified CN(0,C_f) vs identity | seed 2718 rel Fro ~ 0.107 (statistical, O(sqrt(2/n))) |
| B: kernel identity ker K_SLAM = {A u in Range(B)} | **PASS** | two-sided residual + subspace distance | nullity 0 vs 0, subspace distance 0.000e+00 |
| B: rank identity | **PASS** | r(K_IS)-r(K_SLAM) == r(A)+r(B)-r([A,B]) | residual 0 |
| B: no-prior principal-angle retention | **PASS** | rho vs 1 - sigma_i^2(Q_B^T Q_A) | max abs diff 1.776e-15; rho<1 count 18 |
| C: PSD prefix increment at snr2=100 | **PASS** | min eig of K_SLAM/K_eff increments | K_SLAM=1.220e-07, K_eff=7.424e-08 |
| C: rho movement low-SNR -> 0, finite high-SNR | **PASS** | movements on snr2 grid + tail 1e5..1e10 | low (snr2=1e-2): u0: 5.816124e-05 | u1: 9.563285e-05 | u2: 5.532846e-04; tail plateau (snr2=1e10): u0: 1.020462e-02 | u1: 9.420245e-03 | u2: 2.653840e-02 |
| D: duplicate no-prior invariance (same C_f) | **PASS** | rel Fro K_IS/K_SLAM + rho | 1.176e-15, 8.837e-14, 6.216e-13 |
| D: fixed finite prior breaks duplicate rho (expected) | **observed** | max abs rho change, rel Fro K_eff | 3.398424e-01; 2.546359e+00 |
| E: receiver-correlation effect | **recorded** | ell sweep, no forced pass | see table below |

### A raw whitening sanity rows (snr=100, ell=0.10)

| f | sigma_f^2 | cond(C_f) | rel Fro W C W^H - I | max |entry| | MC rel Fro (seed 2718) | MC max |dev| |
|---:|---:|---:|---:|---:|---:|---:|
| 1.00 | 1.011944e-06 | 4.846681e+00 | 7.370057e-16 | 6.870851e-16 | 1.070997e-01 | 6.233686e-02 |
| 1.40 | 2.021532e-06 | 4.846681e+00 | 5.549596e-16 | 6.661338e-16 | 1.070997e-01 | 5.099512e-02 |
| 1.80 | 3.368383e-06 | 4.846681e+00 | 5.769351e-16 | 9.992007e-16 | 1.070997e-01 | 5.828986e-02 |

Monte Carlo draws use 4000 samples per seed with proper CN(0,I) complex
standard normals (`standard_normal @ [1,1j] / sqrt(2)`); the empirical
whitened-realified covariance therefore estimates the identity with sampling
error O(sqrt(2/n)) per entry.  Raw per-seed and per-frequency MC numbers are
in the JSON.

### B algebraic spine (f=1.0, snr=100, ell=0.10, smooth basis)

* Kernel: nullity(K_SLAM)=0, nullity(C)=0,
  subspace distance=0.000e+00,
  min sigma(N1^T N2)=None,
  max ||C u|| over N1=0.000e+00,
  max ||K_SLAM u|| over N2=0.000e+00.
* Rank identity: rA=24, rB=18, rAB=42,
  rKIS=24, rKSL=24, residual=0.
* Retention identity: max|rho - (1 - cos^2)|=1.776e-15;
  count rho < 1 (tol 1e-8)=18;
  rho_min=1.972931e-10.

Both null spaces are trivial (nullity 0 / 0) in this smooth-basis model, as
in the Family 2 smooth case, so the kernel identity holds vacuously with
zero residuals; the near-confounded directions are instead quantified by the
retention spectrum (18 entries below 1 - 1e-8, rho_min ~ 1e-10).

### C SNR sweep: u0 rows (all directions in JSON)

| snr2 | rho_single | rho_stack | movement | shared-z residual |
|---:|---:|---:|---:|---:|
| 0.01 | 1.972932e-10 | 5.816144e-05 | 5.816124e-05 | 7.626365e-03 |
| 0.1 | 1.972932e-10 | 5.799132e-04 | 5.799130e-04 | 2.408139e-02 |
| 1 | 1.972932e-10 | 5.636099e-03 | 5.636099e-03 | 7.507396e-02 |
| 10 | 1.972932e-10 | 4.448916e-02 | 4.448916e-02 | 2.109245e-01 |
| 100 | 1.972932e-10 | 1.277042e-01 | 1.277042e-01 | 3.573573e-01 |
| 1e+03 | 1.972932e-10 | 6.358775e-02 | 6.358775e-02 | 2.521661e-01 |
| 1e+04 | 1.972932e-10 | 2.267252e-02 | 2.267252e-02 | 1.505740e-01 |

The primary grid is snr2 in {1e-2, 0.1, 1, 10, 100, 1e3, 1e4}.  A
supplementary tail grid {1e5, 1e6, 1e8, 1e10} is measured to confirm that
the movements (which overshoot at snr2 ~ 100 and then decline) approach a
finite positive plateau rather than a false endpoint at 1e4.  Tail values
and shared-pose residuals are in the JSON for every direction.

PSD prefix increment at snr2=100:

| matrix | min eig increment | tol | gate |
|---|---:|---:|---|
| K_SLAM | 1.220103e-07 | 9.622e-09 | PASS |
| K_eff | 7.424385e-08 | 9.683e-09 | PASS |

### D duplicate-block invariance (c=2, alpha=1.0, same C_f per block)

no-prior relKIS=1.176e-15, relKSL=8.837e-14, max|rho diff|=6.216e-13; fixed prior max|rho diff|=3.398424e-01, rel K_eff=2.546359e+00

### E receiver-correlation sweep (f=1.0, snr=100)

| ell | sigma_f^2 | cond(C_f) | rel F dist K_SLAM vs identity | rel F dist K_eff vs identity | theta_min deg | retained mass (eig sum) | confusable mass | rho_min | count rho<1e-6 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.001 | 1.011944e-06 | 1.000000e+00 | 9.881964e+05 | 1.402012e+05 | 4.253044e-04 | 3.881282 | 14.592838 | 5.509984e-11 | 5 |
| 0.02 | 1.011944e-06 | 1.059028e+00 | 9.904287e+05 | 1.379861e+05 | 4.357208e-04 | 3.890418 | 14.591157 | 5.783229e-11 | 5 |
| 0.05 | 1.011944e-06 | 2.011894e+00 | 1.081669e+06 | 1.153609e+05 | 5.703777e-04 | 4.252952 | 14.572540 | 9.910127e-11 | 5 |
| 0.1 | 1.011944e-06 | 4.846681e+00 | 1.391756e+06 | 9.615930e+04 | 8.047848e-04 | 5.480059 | 14.541932 | 1.972931e-10 | 5 |
| 0.2 | 1.011944e-06 | 1.189014e+01 | 2.114664e+06 | 9.142263e+04 | 1.083848e-03 | 8.338934 | 14.506195 | 3.578414e-10 | 5 |
| 0.5 | 1.011944e-06 | 3.516602e+01 | 4.330431e+06 | 1.092855e+05 | 1.362746e-03 | 17.097709 | 14.477132 | 5.656968e-10 | 5 |

Distances are relative to the raw identity-noise (Family 4 convention) K
matrices.  The absolute distance mixes the per-SNR sigma_f scale with
receiver decorrelation; rho/principal-angle/mass columns are scale-free.
`retained mass` is the energy-weighted sum of the K_SLAM eigenvalues (Family
6 requirement), while `confusable mass` is the unweighted sum of
principal-angle squared cosines between Range(A) and Range(B) (Family 4
convention); the two use different weightings and do not sum to r_A.

## Figures

* figures/family6_noise_whitening_sanity.png
* figures/family6_snr_sweep.png
* figures/family6_correlation_effect.png
* figures/family6_retention_spectra_colored.png

Figure paths in results JSON: figures/family6_noise_whitening_sanity.png, figures/family6_snr_sweep.png, figures/family6_correlation_effect.png, figures/family6_retention_spectra_colored.png.

## Cannot-establish section

* All statements are finite-dimensional and model-specific: they describe
  the discrete N=16 whitened/realified Jacobians (and their smooth p=24
  projection), not continuum-limit, exact-global-SE(2), recovery, or
  estimator guarantees.
* The covariance `C_f = sigma_f^2 (I_T kron R_rx)` and its Matérn-1/2
  receiver factor are declared synthetic noise models.  No claim is made
  that they describe a physical noise process, a measured noise covariance,
  or an asymptotic estimator.
* The Monte Carlo whitening check establishes finite-sample decorrelation
  of the model at O(sqrt(2/n)) statistical resolution; it cannot certify
  machine-precision whiteness.
* Check C movements and the high-SNR plateau are observed on the three
  confounded directions of the f1=1.0 single block only, for the declared
  block ordering and the finite SNR grid {1e-2..1e4} plus the declared tail
  grid up to 1e10; no universal frequency-diversity theorem is claimed.
* Check D invariance holds only when each duplicated block carries the
  identical C_f/W_f (the tested same-block duplicate).  It does not extend
  to arbitrary receiver resampling, different ell, or model mismatch.
* PSD prefix monotonicity is reported only for the prefix f1 -> {f1,f2}
  at snr2=100; normalized retention spectra are not interpreted as
  monotonicity because K_IS changes with the stack.
* Receiver-correlation rows mix SNR scale and decorrelation in the
  identity-distance column, so only the recorded interpretation (scale-free
  rho/angles/masses) is used for geometric statements.

## Artifacts and digests

```text
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242  src/family4_frequency_trajectory.py
5d3996a3aae094442717b51ab12062487a834313aa808b19048db251861eaa30  src/family6_colored_noise.py
34c417710a318d8157b960f86c5302b952958dc1aa137bea81732117c214e748  results/family6_colored_noise.json
944988ab0138f987b0342534bb6210ebc914b8f7f6f89569c027cbf98b33ab76  figures/family6_noise_whitening_sanity.png
771d3d634b265b5fa4a24d953a60bc63fc32e51ba471c2248de0e923889f544c  figures/family6_snr_sweep.png
8496aef5fcb2566566b29420343d04577ac61eccbb8fff7781dee29d50bb2fad  figures/family6_correlation_effect.png
37d109dd1009ac5aeea3062ed032df063e10c2d63bbb524b836c448b151baf2d  figures/family6_retention_spectra_colored.png
```
