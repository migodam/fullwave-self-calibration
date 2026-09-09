# Matched acquisition-frequency development experiment

Registered 2026-09-08 after the frozen-estimate common-frequency audit. That
audit cannot isolate acquisition changes from differing shared-channel noise.
This new experiment addresses that specific defect, not a final population test.

Known-support ellipsoid, ADDA N64 reference versus N32 inverse CM+RR, one real
material coefficient, three receiver translations, common delay and four
frequency-shared complex illumination gains. Use the same development physical
scenes 8101/8102, fixed loss0.03 and prescribed material2.5/translation90mm.
Parameter conventions and noise scales otherwise follow nonspherical calibration.

Compare three data choices: low=(3,6,9); low+high=(3,6,9,18); and
low+repeat=(3,6,9,6), where the last k6 block is an independent repeated
measurement, not duplicated data. All low-band measurements, including their
noise, are bit-identical across choices. Extra high/repeat blocks use the same
standardized independent noise realization (common random numbers), with the
same absolute proper-complex noise variance fixed from the low-band field RMS
at30dB. This controls measurement count between the two four-block choices.

Test both absent and present noisy electronics reference. When present, the
same reference observations at k=(3,6,9), standard deviation0.02, are given to
ALL data choices; no free additional high-band electronics reference is added.
Same nominal start, physical bounds, least-squares tolerances and35 evaluations.
No priors or true parameters beyond the declared support/model are supplied.

Twelve fits total: two scenes × three data choices × two reference conditions.
Before fits, verify all13 Jacobian columns in a mixed-frequency adapter and
assert exact equality of low-band observation blocks and reference records.
Record every fit/status, setup/solve/evaluation costs and per-frequency work.
Repeated-frequency physical calculations may be reused, but each independent
measurement remains a separate residual block.

Evaluate only frozen estimates on the same k3/k6/k9 held-out receiver channels
against independent noiseless reference means; report sensor and structural
field/phase errors separately, plus material, position and electronic error.
Never use higher training frequency error in the low-band prediction statistic.
Compare low+high against BOTH low and low+repeat. Both helping and harmful
outcomes are retained. Two development scenes support no significance claim;
this is an acquisition/model-fidelity mechanism test, not a new algorithm or
a general theorem that high-frequency phase helps or hurts calibration.
