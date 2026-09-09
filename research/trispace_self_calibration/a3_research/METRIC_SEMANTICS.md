# Prediction and recovery metric ledger

This ledger applies to the A3 three-dimensional experiments. It does not make
the two-dimensional total-field model and three-dimensional scattered-field
model interchangeable.

| Label reserved in the manuscript | Evaluator | Fitted electronics included? | Meaning / source |
|---|---|---|---|
| Oracle-model transfer phase RMSE | Treams at estimated material/geometry, versus Treams truth | No | Old `heldout_phase_rmse_rad` in calibration records; parameter-transfer diagnostic, not deployable prediction |
| Actual sensor prediction | The inverse FFT/VIE model at frozen fitted parameters, versus independent reference sensor field | Yes | `results/prediction_audit3d.json`: `deployable_sensor_prediction`; 17 new receiver angles at radius 1.6 m, no refitting |
| Reconstructed structural field | The inverse FFT/VIE field at frozen fitted material/geometry, versus reference structural field | No | `results/prediction_audit3d.json`: `reconstructed_structural_field`; asks what field/material the inversion recovered, not how well gains fit observations |
| Non-spherical sensor/structural prediction | Inverse ellipsoid model at frozen estimates, versus ADDA-derived reference | Respectively yes/no | `results/nonspherical_calibration.json`; same distinction, known-support one-material experiment |
| Multifidelity endpoint prediction | The method's own deployed inverse resolution, versus ADDA-derived reference | Respectively yes/no | Raw `results/multifidelity_calibration.json`: `deployable_prediction`, `structural_field`; coarse-only uses coarse prediction, but its fine objective is separately audited |
| Measured held-out phase | Fitted cylindrical model, with gains frozen from training source views, versus real data | Yes | `results/measured_extension.json`; no independent clean-field or antenna-position truth exists |

Phase errors use wrapped phase differences with the evaluator's stated field
weights; field errors use normalized complex-field discrepancies. The stored
per-fit quantities are authoritative. Report scene means of stored errors,
not the square root of a newly pooled quantity under the same label.

For FFT comparisons, low-band (3,6,9) and high-band (3,9,18) errors are pooled
over different frequency sets. A larger high-band value does **not** show that
adding high-frequency data makes predictions at the shared low frequencies
worse. That causal claim would require fixed-frequency evaluation.

The previous ambiguous prose has been corrected without overwriting raw
records. A separately qualified oracle diagnostic is useful, but must never
be reported as a physically deployable estimator's held-out prediction.
