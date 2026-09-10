# Q5 V2: audited finite-window research, not a submission-ready TAP paper

Chinese report: REPORT_ZH.md. Read MANUSCRIPT_SECTIONS.md (including Appendix A), EVIDENCE_TABLE.md and PRIOR_ART.md for the exact assumptions and unresolved gates. PROTOCOL.md is the pre-execution remote registration; IMPLEMENTATION_FREEZE.md fixes the implementation. ATTRIBUTION_DIAGNOSTIC_PROTOCOL.md explicitly labels the subsequent oracle work as post-hoc development diagnostics.

A5 and Q5 V1 are unchanged. The failed selector is not reused. Class M resolved modes, class D projected dipole tensor and class C simultaneous lossy spheres are different experiments. The code never promotes a missing class-C continuum bound to certified.

## Reproduce

From the repository root, with the tested Python 3.13 / pinned NumPy and SciPy environment:

```sh
python -m pip install -r research/q5_window_v2/requirements.txt
OPENBLAS_NUM_THREADS=1 python -m pytest -q research/q5_window_v2/tests
```

The proof programs use exact rational arithmetic; SciPy is used in the independent floating diagnostics and numerical forward solver. Existing output files are intentionally not overwritten. To repeat a proof, work in a fresh copy or use its supported alternate output option:

```sh
# Work in a separate checkout/copy; preserve the delivered evidence directory.
cd research/q5_window_v2
mv results results.delivered
mkdir results
python src/certify_modal.py --output results/modal_certificate.json
python src/certify_ratio.py
python src/certify_size.py
python src/finite_geometry.py
python src/continuum_bound.py
python src/certify_budgets.py
python src/certify_secant.py
OPENBLAS_NUM_THREADS=1 python src/validate.py
# Registered development replication, not a new final test:
OPENBLAS_NUM_THREADS=1 python src/recovery.py
```

For a fresh numerical stream with the **same registered seed and settings**, first run `validate.py`, then `recovery.py`; interrupted checkpoints can be continued with `resume_recovery.py`. Running again with that seed is replication, not a new final test. `attribution_diagnostic.py` is a separate post-hoc oracle diagnostic. All programs expect the repository's unchanged `research/trispace_self_calibration/a3_research/maxwell3d.py`.

`modal_readout.py` provides certified numerical brackets for the restricted amplitude and E/M-ratio estimators. Numerical inversion error and physical measurement error are separate. `window_cover.py` implements complete-world interval-pair covering; the caller must supply proved continuum and model-error enclosures. Its class-C status is deliberately unresolved.

## Evidence files

The delivered Q5_V2_RESEARCH_AND_EVIDENCE.zip includes all original raw complex arrays (`observations_v2.npz`), every optimizer start, interrupted/final checkpoints, and logs. The PR stores runnable source, mathematical sections, proof certificates and a lightweight audited snapshot plus cryptographic hashes. A raw archive file is not replaced by a claimed successful rerun. Reproduction may differ in floating last bits across numerical libraries; array hashes identify the exact development record.

No final confirmatory campaign, hardware-calibration certificate, continuous class-C finite material bound, or new acquisition-policy advantage is claimed.

The preserved `review_inputs/REPORT_ZH_pre_audit.md` is an earlier unverified branch narrative. Its median errors and last-decimal budgets do not override the audited snapshot. Other historical registration drafts remain untouched; PROTOCOL.md plus IMPLEMENTATION_FREEZE.md describe the executed main stream. The oracle protocol explicitly postdates inspection of that main development stream.
