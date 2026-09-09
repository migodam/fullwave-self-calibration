# Scale-aware audit of the retained-range pose rank

This is an independent parent-level recomputation. It preserves the locked auto-research outputs and corrects only their scientific interpretation at nuisance-range saturation.

## Rank table

| M | r | rank(H) / 2M | codim | visible pose rank | hidden pose rank | ||P_perp||2 | max sv(P_perp B) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 2 | 7/16 | 9 | 3 | 0 | 1.000e+00 | 8.267e-04 |
| 8 | 3 | 9/16 | 7 | 3 | 0 | 1.000e+00 | 7.783e-04 |
| 8 | 4 | 11/16 | 5 | 3 | 0 | 1.000e+00 | 7.458e-04 |
| 8 | 5 | 13/16 | 3 | 3 | 0 | 1.000e+00 | 5.974e-04 |
| 8 | 6 | 15/16 | 1 | 1 | 2 | 1.000e+00 | 5.679e-04 |
| 8 | 7 | 16/16 | 0 | 0 | 3 | 1.240e-15 | 6.117e-19 |
| 8 | 8 | 16/16 | 0 | 0 | 3 | 1.438e-15 | 4.150e-19 |
| 12 | 2 | 7/24 | 17 | 3 | 0 | 1.000e+00 | 1.336e-03 |
| 12 | 3 | 9/24 | 15 | 3 | 0 | 1.000e+00 | 1.314e-03 |
| 12 | 4 | 11/24 | 13 | 3 | 0 | 1.000e+00 | 1.169e-03 |
| 12 | 5 | 13/24 | 11 | 3 | 0 | 1.000e+00 | 1.006e-03 |
| 12 | 6 | 15/24 | 9 | 3 | 0 | 1.000e+00 | 9.616e-04 |
| 12 | 7 | 17/24 | 7 | 3 | 0 | 1.000e+00 | 8.187e-04 |
| 12 | 8 | 19/24 | 5 | 3 | 0 | 1.000e+00 | 4.610e-04 |
| 12 | 9 | 21/24 | 3 | 3 | 0 | 1.000e+00 | 2.887e-04 |
| 12 | 10 | 23/24 | 1 | 1 | 2 | 1.000e+00 | 2.146e-04 |
| 12 | 11 | 24/24 | 0 | 0 | 3 | 1.713e-15 | 6.163e-19 |
| 12 | 12 | 24/24 | 0 | 0 | 3 | 1.305e-15 | 1.129e-18 |
| 16 | 2 | 7/32 | 25 | 3 | 0 | 1.000e+00 | 1.527e-03 |
| 16 | 3 | 9/32 | 23 | 3 | 0 | 1.000e+00 | 1.355e-03 |
| 16 | 4 | 11/32 | 21 | 3 | 0 | 1.000e+00 | 1.184e-03 |
| 16 | 5 | 13/32 | 19 | 3 | 0 | 1.000e+00 | 1.151e-03 |
| 16 | 6 | 15/32 | 17 | 3 | 0 | 1.000e+00 | 1.107e-03 |
| 16 | 7 | 17/32 | 15 | 3 | 0 | 1.000e+00 | 9.243e-04 |
| 16 | 8 | 19/32 | 13 | 3 | 0 | 1.000e+00 | 5.167e-04 |
| 16 | 9 | 21/32 | 11 | 3 | 0 | 1.000e+00 | 2.307e-04 |
| 16 | 10 | 23/32 | 9 | 3 | 0 | 1.000e+00 | 1.841e-04 |
| 16 | 11 | 25/32 | 7 | 3 | 0 | 1.000e+00 | 1.428e-04 |
| 16 | 12 | 27/32 | 5 | 3 | 0 | 1.000e+00 | 1.322e-04 |
| 16 | 13 | 29/32 | 3 | 3 | 0 | 1.000e+00 | 5.147e-05 |
| 16 | 14 | 31/32 | 1 | 1 | 2 | 1.000e+00 | 2.643e-05 |
| 16 | 15 | 32/32 | 0 | 0 | 3 | 1.464e-15 | 1.515e-18 |
| 16 | 16 | 32/32 | 0 | 0 | 3 | 1.840e-15 | 9.855e-19 |

The transition follows a dimension obstruction, not a special Fourier/Hankel theorem. With p=3 real pose variables and K=3 real contrast variables, the generic nuisance column count is 2r+K in a 2M-dimensional real data space. At r=M-2 the orthogonal complement has dimension at most one, so at least two pose directions are hidden. At r>=M-1 a full-rank nuisance span fills the data space, so all three pose directions are hidden; the old relative-only rule mistook roundoff for three visible directions.

## Corrected nonlinear control

| setting | method | regime | success | median pose error | median map error | median residual |
|---|---|---|---:|---:|---:|---:|
| M12_k12_full_corrected | direct | direct_control | 3/3 | 7.6126e-03 | 1.2055e-02 | 2.7133e-02 |
| M12_k12_full_corrected | reduced_r10 | codim_one | 3/3 | 7.5896e-02 | 4.6321e-02 | 7.8702e-02 |
| M12_k12_full_corrected | reduced_r11 | nuisance_saturated | 3/3 | 1.1180e-01 | 2.5978e-01 | 3.1690e-01 |
| M16_k12_full_corrected | direct | direct_control | 3/3 | 1.4754e-02 | 1.5039e-02 | 2.4853e-02 |
| M16_k12_full_corrected | reduced_r14 | codim_one | 3/3 | 1.1096e-01 | 8.9523e-02 | 1.6938e-01 |
| M16_k12_full_corrected | reduced_r15 | nuisance_saturated | 3/3 | 1.1180e-01 | 2.4030e-01 | 3.1536e-01 |

A reduced solver with zero visible pose directions freezes pose at the nominal value. Therefore equality with the direct solver at r=M-1 in the locked report was a numerical-rank artifact, not evidence of lossless reduction.

## Scientific consequence

The reliable result is a pose-hiding diagnostic and a dimension bound. It is not a hidden-direction recovery algorithm. Full pose reparameterization (visible rank three) is algebraically equivalent to direct joint optimization, while truncating to the visible row space intentionally removes hidden directions.
