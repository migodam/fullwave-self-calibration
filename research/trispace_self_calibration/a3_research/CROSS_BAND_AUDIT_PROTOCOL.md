# Frozen-estimate common-frequency prediction audit

Registered 2026-09-08 before evaluating per-frequency predictions. This is
post-fit DEVELOPMENT analysis of the eight existing h=0.01 m FFT estimates;
no new fitting, tuning or claim of an unopened test population.

Question: do the (3,9,18) fits improve actual predictions at the same k=3 or 9
channels compared with the (3,6,9) fits? Previously reported group averages
mix different frequency sets and cannot answer that question.

For seeds 6101/6102 and unified/unified-reference variants, freeze every stored
estimate. Evaluate the actual inverse model with fitted electronics at the
same 17 held-out receiver locations, against the independent noiseless Treams
mean with true electronics. Also evaluate the recovered structural field
without electronics. Report field error and truth-amplitude-weighted wrapped
phase RMSE separately at every frequency. Compare only common k=3 and 9;
do not pool k=6 with k=18 under a common-channel label.

Important causal limit: the higher group REPLACES k=6 with k=18, rather than
adding high-frequency data to identical low-band measurements. Some shared
frequency noise draws also differ because frequency ordering changes. These
frozen runs cannot isolate a pure high-frequency causal effect. They can reveal
whether the existing manuscript's proposed common-frequency benefit is even
consistent with these estimates. A causal claim needs a new matched-noise,
matched-acquisition experiment, with all low-band observations retained and
an equal-count alternative acquisition control.
