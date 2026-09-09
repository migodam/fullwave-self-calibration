# A4 development log

This is a new implementation, not a reproduction of unavailable A3 scripts.
The old A3 files and frozen outcomes have not been modified.

Development v1: seed 81031 (six branch scenes) and 82031 (two calibration scenes).
The complete outputs are retained. Review found three issues before freezing:

1. The deliberate two-world displacement was 12 mm, below the declared 15 mm
   task tolerance. It could not test false acceptance beyond tolerance. The
   revised adversary has 25 mm displacement within a declared 40 mm prior box.
   This is a development hard-control correction, not a changed success threshold.
2. The approximation-error mean correction had the wrong sign in the residual.
   It is now `coarse_prediction + estimated_error - observation`. A regression
   test checks exact cancellation. The v1 comparison must not support baseline claims.
3. Error-mode methods were being charged screening calls they did not use.
   Construction and screening call counts are now separate. Pilot reference
   noise is also independently drawn and charged, rather than reused at fitting.

Development v2: fresh seeds 81032 and 82032, before final protocol and hashes.
All 21 unit/regression checks pass. In these development scenes the controller
selects the full fine model and is slower than directly using it. That negative
outcome is not grounds for changing the frozen target metric or removing the
full-fine baseline. The conditional-error baseline is retained even when worse.

The two-world control defeats the conditional beta=0 coverage mechanism as
expected: it supplies data from another admissible exact full-wave response.
Its failure is outside the zero-discrepancy premise, not a contradiction of the
conditional theorem. No software gate can infer its absence from residuals.
