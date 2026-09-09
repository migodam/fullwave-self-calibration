# Provenance, frozen integrity, and correction ledger

## Input boundary

The supplied GitHub main snapshot contains authored A3 documents and figures, not its executable source/results. The observed recursive tree SHA was `0fd9ca488f1e0b19edc27d5492a37c3130f00a16` (a tree, not a commit). The README explicitly states the reproduction gap. Connector reads retrieved the README, main A3 theory/evidence discussions and several supporting audits. Direct git clone failed due to container DNS; no A3 source was silently reconstructed and labelled original.

Six uploaded Markdown files are separately hashed in `provenance/input_manifest.json`. A3 findings quoted in this package are inherited report-level evidence, not newly reproduced raw A3 outputs. Missing original `THEORY_AUDIT.md`, A2 source/theorem/validation package and original A3/A3_2 question attachments were not read. The full Chen monograph was not reread; limited Library retrieval does not close original SOM-paper reading requirements.

## A4 freeze

`provenance/frozen_manifest.json` predates the fresh collection and hashes all frozen source files, tests and `FROZEN_PROTOCOL.md`. This is local registration, not an externally timestamped preregistration. `results/A4_AUDIT.json` checks equality. Analysis/manuscript/supplement files are not retroactively called frozen.

The final numerical tables are read directly from current raw JSON by `analyze_a4.py`; cached prose or interim summaries are not authoritative. Main collection runtime is the value in `frozen_complete.json`. All scene/policy rows, including negative controls and failures, are retained. An execution end marker alone is not scientific acceptance.

## Preserved corrections

1. Development-only mean-correction sign, independent pilot/reference use and screening accounting were corrected before the final freeze. Original development artifacts remain and are not pooled with fresh scenes.
2. Frozen reference **diagnostic** row removal was wrong: realified reference rows are noncontiguous. Its no-reference tangent metric is withdrawn; original fit/recovery values and original frozen source remain. The correction is a new source/protocol/fresh sample in `supplements/reference_audit_v2/`; it does not modify the frozen experiment.
3. A default post-freeze pytest command imported an unrelated installed package named `coverage`. Its failure trace remains. Disabling plugin autoload and using `PYTHONPATH=src` runs the intended modules and passes21 cases. This environment workaround does not change frozen code.
4. Exact exterior-model error identities were checked after freeze in a separate supplement and labelled as such. They are neither a new frozen test nor a DDA continuum certificate.
5. No raw frozen endpoint or method was altered to improve a table. No failed scene is dropped from the review.

## Source and execution scope

All A4 scientific code is a self-contained implementation using NumPy/SciPy; it does not copy a third-party solver. In-house DDA is independent of the analytic Mie implementation in numerical formulation, not a claim of external-code independence. No treams or ADDA execution occurred in A4. No hardware data were generated. User-supplied A3 documents are not redistributed inside the code ZIP; their filenames/hashes establish input provenance.

The archive is a local overlay for `research/a4_reliability_v1/`. No branch, commit or pull request was pushed to the remote repository during this delivery; existing main/frozen A3 records are untouched.

## Precise reference diagnostic correction

The withdrawn expression is `J[:-6]`. In this 30-complex-row residual, realification stacks 30 real rows followed by 30 imaginary rows, so the 27 data rows occupy realified indices `np.r_[np.arange(27), np.arange(30,57)]`. Removing the last six rows alone keeps real reference rows and drops some imaginary data rows. The correction changes only the separately versioned fresh diagnostic, not the original frozen source.

The first packaging audit failed its literal documentation check because this exact expression had not yet been written in this provenance note; numerical/hash checks already passed. That packaging report is retained in `results/packaging_first_audit.json`.
