# Minimal controlled hardware protocol -- not executed

Application: coherent near-field microwave/diffraction tomography relative to a measured, anchored object support. This is not UAV/SLAM validation.

Use a phase-coherent VNA/receiver and a positioner with independent metrology; measure effective vector receive response and the complex incident field in the target region. A radial calibration body with independently characterized radius/center is the minimum target for the present theorem. Three synthesized regular electric-l=1 fields must be validated over that region, including unwanted electric/magnetic modal leakage. A fixed real antenna orientation and a scalar gain shared across each3-by-3 tensor must be tested rather than assumed.

Reserve independent training and validation positions/records before fitting. Apply controlled receiver displacements, controlled cable-delay changes and controlled gain changes separately and in crossed combinations. Positioner/metrology truth is not passed to the inverse solver. Electronics references should bypass the scattering path and have an independently estimated noise/transfer model. Validation data must not refit the geometry.

The decisive controls are: geometry shift recovers position without fictitious material/gain compensation; cable-delay and gain shifts do not produce fake pose; poor incident-field/exterior-model fidelity triggers a qualified reject/refine rather than confidence; known dithers resolve the tensor's sign ambiguity; frequency-only repeats do not; component-gain corruption defeats the scalar-gain theorem unless independently calibrated.

Report complex-data counts, vector-component switch count, synthesized-source degrees of freedom, stage settling time, reference acquisition time and independent metrology uncertainty. A software comparison with nine tensor entries is not a nine-cheap-measurement hardware claim. Predeclare a leakage/discrepancy bound from measured fields or a separately calibrated statistical population; adjacent-grid differences are insufficient.

No apparatus, feasible VNA band, antenna hardware price, or measured error level is asserted here. The present length scales and dimensionless kR law need a consistent physical frequency/size rescaling for a real facility. A known-target calibration run does not validate arbitrary unknown-target material tomography. Existing Fresnel data without antenna displacement truth cannot substitute for this protocol.
