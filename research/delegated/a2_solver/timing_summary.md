# A2 E4 tuning-only timing microbenchmark (N16)

Raw samples are in `calibration_tune.json` (`raw_samples`). Median per
operation, conservative overhead factor 1.15 is applied when converting
reduced operations to full-RHS-equivalent work.

| quantity | measured seconds |
|---|---:|
| full RHS solve (LU, per column) | 1.117e-5 |
| n x n LU factorisation | 4.771e-4 |
| dense n x n matvec | 4.479e-6 |
| reduced QR n x 8 | 3.562e-5 |
| reduced QR n x 16 | 1.159e-4 |
| reduced QR n x 36 | 4.379e-4 |
| reduced QR n x 64 | 1.166e-3 |
| reduced QR n x 128 | 3.710e-3 |
| reduced QR n x 256 | 4.058e-3 |
| reduced RHS solve r=8 (after QR) | 8.70e-7 |
| reduced RHS solve r=64 | 2.30e-6 |
| reduced RHS solve r=256 | 9.43e-6 |
| sensing-stack SVD m=36 | 6.34e-4 |
| sensing-stack SVD m=144 | 6.05e-3 |
| RRQR n x 64 | 1.24e-3 |

Budget accounting charges full-wave RHS columns one unit each; reduced QR
and reduced RHS solves are converted using these conservative measured
ratios. LU factorisations, SVD/RRQR basis work, dense operator products
and wall time are counted and reported separately (raw counts and wall
seconds in every JSONL run record).

Typical tuning run work (budget 800): direct/phaseless/direct_control
consumed 780 units, prasc 795.6 units, fixed_rank 674.4 units; mean wall
times were 0.68/0.44/0.39/0.72/0.69 s respectively.
