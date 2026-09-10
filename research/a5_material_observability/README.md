# A5 material-observability continuation

Latest A4 closure handoff: [Chinese complete report](../../communication/A4_COMPLETE_RESEARCH_REPORT_ZH.md),
[Q5 GPT Pro total prompt](../../Theory/Questions/Q5.md), and
[English manuscript addendum](A4_CLOSURE_ADDENDUM_EN.md).
New six-scale causal control: [protocol](ROUND2_PROTOCOL.md) and
[raw results](results/interaction_ablation.json). Isolated-sphere sums retain
slightly greater profiled scale sensitivity than interacting clusters in these
six cases; inter-sphere enhancement is not supported here. Research remains active.

Read [integration and corrections](THEORY_INTEGRATION.md), then the
[registered plan](PLAN_AND_PROTOCOL.md). Source inputs: A4 and A4_2.

Executed: A4 regression21/21;9 exact modal counterexamples;18 two-sphere vector-
Maxwell cases with lmax3/4 comparisons and three gain-sharing patterns.
Results: [scale_gauge.json](results/scale_gauge.json).

Run `scale_gauge_experiment.py` with the existing A3 `.venv3d` environment.
It refuses to overwrite an existing result. A3's `maxwell3d.py` and the byte-
preserved A4 source in `public_release/` are imported, not modified. This is
not a new installed Agentic-AI-Scientist run or a full material inversion study.

Open: material-recovery design, multi-region recovery, interaction ablations,
gain constraints, source/model uncertainty, nearest-prior proof-level review,
and isolated GLM support in the installed pipeline. No TAP-readiness claim.
