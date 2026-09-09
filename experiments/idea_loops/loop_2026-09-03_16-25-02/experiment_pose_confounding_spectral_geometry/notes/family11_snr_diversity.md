# Family 11: SNR diversity across scenes and frequency sets

Date: 2026-09-03T14:52:12.113744+00:00 UTC.  Experiment: `experiment_pose_confounding_spectral_geometry`.

Command: `.venv/bin/python src/family11_snr_diversity.py`.  Wall runtime: 4.12 s.  Platform: macOS-26.6.2-arm64-arm-64bit, Python 3.12.13, numpy 2.5.2, matplotlib 3.11.1.

## Purpose and scope
Replicates the family 4/6 non-monotonic SNR-dependent frequency-diversity curve (movement of the three most pose-confounded single-frequency directions as a second frequency's SNR varies), including the finite high-SNR plateau, over three scenes and f2 in {1.4, 1.8}.  It then characterises equal per-frequency-SNR stacks F1/F2/F3/F5 and a duplicate-block control, plus a fixed-count bandwidth sweep.  Everything is deterministic dense linear algebra on the discrete N=16 2D scalar Helmholtz toy (identity complex noise, SNR-whitened exactly as family10/family6 ell->0); no forced pass.

## Part A summary: per (scene, f2) series

| scene | f2 | direction | rho_single | peak movement | peak snr2 | tail plateau (median) | non-monotone |
|---|---|---|---|---|---|---|---|
| two_blob | 1.4 | u0 | 5.510e-11 | 1.257126e-01 | 100 | 5.011343e-03 | yes |
| two_blob | 1.4 | u1 | 1.414e-10 | 2.334654e-01 | 100 | 4.315232e-03 | yes |
| two_blob | 1.4 | u2 | 1.096e-09 | 5.756348e-01 | 100 | 1.261219e-02 | yes |
| two_blob | 1.8 | u0 | 5.510e-11 | 4.342802e-01 | 1e+03 | 1.357545e-02 | yes |
| two_blob | 1.8 | u1 | 1.414e-10 | 1.612308e-01 | 1e+03 | 2.105576e-02 | yes |
| two_blob | 1.8 | u2 | 1.096e-09 | 4.596378e-01 | 316 | 6.895528e-02 | yes |
| ring | 1.4 | u0 | 9.962e-11 | 5.988487e-01 | 316 | 6.068653e-02 | yes |
| ring | 1.4 | u1 | 4.547e-10 | 7.177072e-01 | 1e+03 | 7.126725e-01 | yes |
| ring | 1.4 | u2 | 7.265e-10 | 4.053000e-01 | 316 | 1.989163e-01 | yes |
| ring | 1.8 | u0 | 9.962e-11 | 7.166836e-01 | 1e+03 | 6.310011e-02 | yes |
| ring | 1.8 | u1 | 4.547e-10 | 4.937739e-01 | 316 | 4.697055e-01 | yes |
| ring | 1.8 | u2 | 7.265e-10 | 4.307602e-01 | 1e+03 | 1.446517e-02 | yes |
| low_contrast | 1.4 | u0 | 8.159e-11 | 1.066362e-01 | 100 | 6.172201e-03 | yes |
| low_contrast | 1.4 | u1 | 7.357e-10 | 4.708391e-01 | 316 | 3.407395e-02 | yes |
| low_contrast | 1.4 | u2 | 7.029e-09 | 6.525455e-01 | 31.6 | 1.276499e-02 | yes |
| low_contrast | 1.8 | u0 | 8.159e-11 | 2.499464e-01 | 3.16e+03 | 1.979080e-01 | yes |
| low_contrast | 1.8 | u1 | 7.357e-10 | 4.881411e-01 | 1e+04 | 3.231883e-01 | no |
| low_contrast | 1.8 | u2 | 7.029e-09 | 3.998627e-01 | 316 | 7.156373e-02 | yes |

Peak is over the 13-point primary grid `logspace(-2,4,13)`; non-monotone = primary peak above the tail median + 1e-6 at a non-endpoint snr2.  Full per-snr2 rows, shared-z residuals, and tail rows are in the JSON.

## Part B summary: equal-SNR frequency sets

| scene | set | retained DOF (sum rho) | log-volume | K_eff alpha |
|---|---|---|---|---|
| two_blob | F1 | 11.2929456138 | -61.6645879559 | alpha=1 |
| two_blob | F2 | 13.4812396301 | -47.6731635579 | alpha=1 |
| two_blob | F3 | 15.6364688685 | -24.5372028097 | alpha=1 |
| two_blob | F5 | 15.9117483745 | -20.7996945996 | alpha=1 |
| two_blob | DUP | 11.2929456138 | -61.6645879559 | alpha=2 |
| ring | F1 | 15.8438774096 | -54.0528714486 | alpha=1 |
| ring | F2 | 17.9895170810 | -39.5748335894 | alpha=1 |
| ring | F3 | 19.2602533997 | -13.0875695451 | alpha=1 |
| ring | F5 | 20.0568946868 | -9.2771962668 | alpha=1 |
| ring | DUP | 15.8438774096 | -54.0528714486 | alpha=2 |
| low_contrast | F1 | 16.9803786027 | -24.1395544193 | alpha=1 |
| low_contrast | F2 | 17.1141978414 | -21.5175388097 | alpha=1 |
| low_contrast | F3 | 18.6183940010 | -12.0463544925 | alpha=1 |
| low_contrast | F5 | 18.9645532086 | -12.2629144947 | alpha=1 |
| low_contrast | DUP | 16.9803786027 | -24.1395544193 | alpha=2 |

K_eff spectra use alpha=1 (family2/4 retention convention).  DUP retained DOF/log-volume are reported under the invariant joint-scaled prior alpha=(1+c^2)*alpha=2, which reproduces the F1 alpha=1 generalized K_eff spectrum; the fixed-alpha=1 diagnostic is recorded separately in the JSON.  Direction movements for each set are in the JSON.

## DUP duplicate-block invariance residuals

| scene | max|movement| | no-prior DOF rel | K_eff joint DOF rel | K_eff joint logvol rel | 1e-10 gate |
|---|---|---|---|---|---|
| two_blob | 1.256e-16 | 4.305e-14 | 1.258e-15 | 1.224e-13 | PASS |
| ring | 2.682e-16 | 7.832e-15 | 7.400e-15 | 2.369e-13 | PASS |
| low_contrast | 2.187e-16 | 3.699e-13 | 7.951e-15 | 2.105e-14 | PASS |

No-prior log-volume relative residuals are dominated by logs of near-floor rho entries (absolute log-volume differences ~1e-7, relative ~1e-9) and are recorded as diagnostics; the invariance gate uses movement, no-prior max|rho| difference and retained-DOF difference, and the joint-scaled-prior K_eff DOF/log-volume relative differences (all < 1e-10 where the criteria apply).

## Bandwidth sweep (always 3 frequencies)

| scene | set | frequencies | retained DOF | log-volume |
|---|---|---|---|---|
| two_blob | BW_narrow | 1-1.2-1.4 | 14.7089136543 | -31.6527044226 |
| two_blob | BW_mid | 1-1.4-1.8 | 15.6364688685 | -24.5372028097 |
| two_blob | BW_wide | 1-1.8-2.6 | 14.5190115212 | -30.0564229605 |
| ring | BW_narrow | 1-1.2-1.4 | 18.7995832565 | -20.4486125160 |
| ring | BW_mid | 1-1.4-1.8 | 19.2602533997 | -13.0875695451 |
| ring | BW_wide | 1-1.8-2.6 | 18.6136890422 | -15.8580000478 |
| low_contrast | BW_narrow | 1-1.2-1.4 | 17.9705531225 | -16.1341232916 |
| low_contrast | BW_mid | 1-1.4-1.8 | 18.6183940010 | -12.0463544925 |
| low_contrast | BW_wide | 1-1.8-2.6 | 18.0644040625 | -14.3304575974 |

## Part C summary claims (reported, no forced pass)

| scene | f2 | directions non-monotone | peak+plateau replicates | peak movements | peak snr2s | tail plateau medians |
|---|---|---|---|---|---|---|
| two_blob | 1.4 | 3/3 | True | 1.257e-01, 2.335e-01, 5.756e-01 | 100, 100, 100 | 5.011e-03, 4.315e-03, 1.261e-02 |
| two_blob | 1.8 | 3/3 | True | 4.343e-01, 1.612e-01, 4.596e-01 | 1e+03, 1e+03, 316 | 1.358e-02, 2.106e-02, 6.896e-02 |
| ring | 1.4 | 3/3 | True | 5.988e-01, 7.177e-01, 4.053e-01 | 316, 1e+03, 316 | 6.069e-02, 7.127e-01, 1.989e-01 |
| ring | 1.8 | 3/3 | True | 7.167e-01, 4.938e-01, 4.308e-01 | 1e+03, 316, 1e+03 | 6.310e-02, 4.697e-01, 1.447e-02 |
| low_contrast | 1.4 | 3/3 | True | 1.066e-01, 4.708e-01, 6.525e-01 | 100, 316, 31.6 | 6.172e-03, 3.407e-02, 1.276e-02 |
| low_contrast | 1.8 | 2/3 | False | 2.499e-01, 4.881e-01, 3.999e-01 | 3.16e+03, 1e+04, 316 | 1.979e-01, 3.232e-01, 7.156e-02 |

Replicates both f2 values: two_blob = True, ring = True, low_contrast = False.

## Config and noise model

N=16, T=6, n_rx=4, m_real=48 per frequency, arc radius 1.6, phi -45.0..45.0, rx offsets [[-0.06, 0.0], [0.06, 0.0], [0.0, -0.06], [0.0, 0.06]], tx [0.0, 0.0].  Smooth basis p=24, unit columns.  alpha=1.0.  Scene definitions: {'two_blob': 'family1.make_chi0(points, cfg) with standard blobs', 'ring': 'chi = 0.5*exp(-((|r|-0.25)/0.05)^2), |r| over N=16 cell centres', 'low_contrast': 'chi = 0.1 * two_blob chi0'}.

## Figures

* [family11_frequency_sets.png](figures/family11_frequency_sets.png)
* [family11_snr_sweep_scenes.png](figures/family11_snr_sweep_scenes.png)

## Artifacts and digests

```text
942d479c283fdab7f3135e79350e8dce491627d62c7f611a4ed8fce9c3705944  figures/family11_frequency_sets.png
d4b90262328212fd1b0b13c71a1cb2e24d2bd280c19029a612df55cb57a4369f  figures/family11_snr_sweep_scenes.png
d9bf44749fe5c554333eda89639221ee1cd04c66fd5acde91a012cabdc07df7c  results/family11_snr_diversity.json
44a976104b0fa329a5699ce8c048e2fa3ed125115465bd6807ad117fb303e45b  src/family10_online_slam_toy.py
7b7cb6337f4ed83b88e91baf9be96d49d28ad9f6dc04001622340c106be423ad  src/family11_snr_diversity.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242  src/family4_frequency_trajectory.py
5d3996a3aae094442717b51ab12062487a834313aa808b19048db251861eaa30  src/family6_colored_noise.py
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
```
