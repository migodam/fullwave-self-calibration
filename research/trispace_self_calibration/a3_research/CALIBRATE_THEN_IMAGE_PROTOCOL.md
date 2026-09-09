# Registered calibrate-then-image control

2026-09-08, after the 32 shape/material stress outcomes. Development on those
same four cases, not confirmatory evidence. Plug-in two-stage estimation is not
claimed as a new general method.

Question: can the useful high-band position estimate be transferred to a
low-band material fit without retaining the harmful high-band material bias?

For each seed8401–8404 and each electronic-reference condition, freeze geometry
from each completed warm_raw, isotropic and rank1 continuation. Starting the
remaining parameters from the SAME corresponding low-only pilot, fit epsilon,
common delay, four log-gains and four phases to the unchanged low=(3,6,9) data
and the unchanged optional reference. Geometry is held exactly fixed throughout.
No high data are used in the second objective, but they informed the frozen
geometry; report that dependence rather than pretending the output is low-only.

Add a separately labelled true-geometry control with all remaining parameters
initialized identically. Truth is permitted ONLY in this oracle arm, never to
select the operative geometry, nuisance initialization, weights or output.
The oracle is diagnostic, not necessarily a lower bound under model discrepancy.

32 second-stage fits: four cases × two references × four geometry sources.
Use the existing nuisance-coordinate bounds/scales, 35 evaluations and tolerances.
No additional starts, weighting, tuning or rescue fits. All input endpoints and
NPZs are immutable. Test the reduced 10-column Jacobian and exact fixed geometry
before collection. Preserve endpoints before evaluation; keep failures with stages,
atomic records, exclusive lock and provenance guards. Run serially after the
stress benchmark is complete. Charge prior pilot/calibration costs separately
from the new image stage; historical sums are not fresh wall-time measurements.

Primary contrasts: rank1-geometry second stage against direct rank1 and original
low-only. Raw-geometry and isotropic-geometry second stages test whether the
benefit depends on the calibration method. Report every case's material,
held-out low-band structural/sensor field and phase; geometry must be unchanged.
No uncertainty claim follows from fixing an estimated parameter. If material
does not improve, or new bias replaces the old one, retain the failure and do
not promote the two-stage transfer to a core contribution.
