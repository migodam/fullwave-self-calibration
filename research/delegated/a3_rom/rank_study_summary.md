# A3 fixed-chart rank-error study (development, grid16)

Timestamp: 2026-09-07T12:30:40+08:00
Wall: 40.77 s. Three material samples (seeds 4101-4103), k=3,6,12 rad/m.

Thresholds: state/output rel error <= 1e-3; map/pose Jacobian rel Frobenius error <= 1e-2. First pass is the lowest ladder rank that passes all four on that sample/frequency. Actual rank is shown when it differs from the target (sensing is capped by its 36-row SVD).

| Method | k=3 first-pass ranks | k=6 | k=12 |
|---|---|---|---|
| sensing | 未通过至256 | 未通过至256 | 未通过至256 |
| twofold | 96 / 96 / 96 | 128 / 128 / 128 | 160 / 192 / 160 |
| generic_task | 32 / 48 / 32 | 64 / 64 / 64 | 96 / 128 / 96 |
| block_krylov | 96 / 192 / 96 | 160 / 256 / 160 | 256 / 256 / 256 |
| sensing_task | 96 / 96 / 96 | 96 / 96 / 96 | 96 / 128 / 96 |
| twofold_task | 96 / 96 / 96 | 96 / 96 / 96 | 128 / 128 / 128 |

## Notes and limitations

- Chart construction never uses exact states/derivatives; exact forward/Jacobian solves are offline audits only.
- The Twofold seed is the implementable sequential projection ``[Vs,(I-Ps)Vd]``, not an exact intersection.
- ``sensing`` is capped by the numerical rank of the 36-row sensing matrix; the sensing/twofold rows labelled *_task include the same task-residual/tangent-RHS enrichment budget as ``generic_task``.
- D right singular blocks are frequency-fixed and cached once; the single setup SVD per frequency is reported separately, not charged to each chart.
- Raw rows are in ``rank_study_raw.json``. This is a development rank-error audit, not a final sweep and not SOM falsification.
