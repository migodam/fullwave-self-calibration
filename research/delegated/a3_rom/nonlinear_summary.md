# A3 development nonlinear ROM runs (grid16)

Timestamp: 2026-09-07T12:29:47+08:00; total wall 200.8 s; 56 runs.

All runs use identical final multifrequency data (k=3,6,12). ``cont=1`` is the low-to-high stage schedule; ``cont=0`` is the same data added at once. Reduced acceptance counts are exact-objective accepted ROM steps, grouped by every frequency active in that stage.

| method / continuation | n | final loss mean | pose m mean | alpha RMSE mean | accepted reduced mean | fallback mean | full RHS mean | wall s mean | statuses |
|---|---|---|---|---|---|---|---|---|---|
| block_krylov::cont=0 | 4 | 213.42 | 0.0314 | 0.0293 | 6.8 (k0:6.8, k1:6.8, k3:6.8) | 2.0 | 1634 | 6.07 | no_improving_trial_and_fallback |
| block_krylov::cont=1 | 4 | 214.50 | 0.0301 | 0.0303 | 17.5 (k0:17.5, k1:9.5, k3:4.0) | 2.5 | 1809 | 6.96 | no_improving_trial_and_fallback |
| direct_adjoint::cont=0 | 4 | 213.58 | 0.0294 | 0.0270 | 0.0 | 0.0 | 1026 | 1.07 | schedule_complete |
| direct_adjoint::cont=1 | 4 | 213.72 | 0.0314 | 0.0264 | 0.0 | 0.0 | 1089 | 1.17 | schedule_complete |
| direct_gn::cont=0 | 4 | 213.42 | 0.0314 | 0.0293 | 0.0 | 0.0 | 2070 | 0.87 | no_improving_trial |
| direct_gn::cont=1 | 4 | 217.08 | 0.0403 | 0.0344 | 0.0 | 0.0 | 2652 | 1.78 | no_improving_trial |
| generic_task::cont=0 | 4 | 213.42 | 0.0314 | 0.0293 | 6.5 (k0:6.5, k1:6.5, k3:6.5) | 0.2 | 792 | 3.90 | no_improving_trial_and_fallback |
| generic_task::cont=1 | 4 | 217.08 | 0.0402 | 0.0344 | 15.5 (k0:15.5, k1:7.5, k3:6.0) | 0.2 | 1011 | 5.22 | no_improving_trial_and_fallback |
| sensing::cont=1 | 4 | 216.26 | 0.0528 | 0.0323 | 10.2 (k0:10.2, k1:2.2) | 7.5 | 3231 | 4.70 | no_improving_trial_and_fallback |
| sensing_task::cont=0 | 4 | 213.42 | 0.0314 | 0.0293 | 6.5 (k0:6.5, k1:6.5, k3:6.5) | 0.8 | 990 | 3.68 | no_improving_trial_and_fallback |
| sensing_task::cont=1 | 4 | 217.08 | 0.0402 | 0.0344 | 15.2 (k0:15.2, k1:7.2, k3:5.0) | 0.2 | 1006 | 4.63 | no_improving_trial_and_fallback |
| twofold::cont=1 | 4 | 217.08 | 0.0403 | 0.0344 | 15.2 (k0:15.2, k1:7.2, k3:5.0) | 0.0 | 908 | 2.83 | no_improving_trial_and_fallback |
| twofold_task::cont=0 | 4 | 213.42 | 0.0314 | 0.0293 | 6.5 (k0:6.5, k1:6.5, k3:6.5) | 0.2 | 792 | 3.12 | no_improving_trial_and_fallback |
| twofold_task::cont=1 | 4 | 217.08 | 0.0402 | 0.0344 | 15.2 (k0:15.2, k1:7.2, k3:5.0) | 0.2 | 1006 | 4.19 | no_improving_trial_and_fallback |

## Interpretation boundaries

- A ``no_improving_trial_and_fallback`` stop at low exact loss is the exact acceptance controller refusing further numerical improvement; it is not a solver failure by itself.
- Accepted ROM steps in the single all-frequency schedule directly include k=12 (model frequency 3); continuation accepted-step counts are split by active stage frequency.
- Chart rebuild is charged to every accepted iterate; full/reduced work are not merged into a speedup claim.
- These are development outcomes only (seeds 4101-4104). No final sweep, no SOM-falsification claim.
