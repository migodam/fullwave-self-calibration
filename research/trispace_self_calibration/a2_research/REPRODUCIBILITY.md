# Reproduction and evidence index

Project root: `/Volumes/migodam's-external-brain/Research/Inv_SLAM`.
All relative paths below are from this root unless stated otherwise. CPU-only,
one BLAS thread; no GPU or hardware claims. Python 3.13.13, NumPy 2.5.2,
SciPy 1.18.1, Matplotlib 3.11.1 in the existing experiment environment.

## Quick verification (does not rerun final optimizations)

```sh
cd "/Volumes/migodam's-external-brain/Research/Inv_SLAM"
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
A2_PY="experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python"
"$A2_PY" research/trispace_self_calibration/a2_research/verify_delivery.py
"$A2_PY" research/trispace_self_calibration/a2_research/analyze_final.py
```

`results/delivery_verification.json` contains exact executed commands, exit codes,
logs, final-code hash checks and paired larger-map count checks. Logs are local
bounded test outputs, not model reasoning or credential traces. Green checks
mean implementation integrity, not method superiority or journal readiness.

## Study evidence

| Component | Canonical records | Interpretation |
|---|---|---|
| Original Pro inputs | `Theory/Questions/A2.md`, `A2_PRASC_SOM_THEOREM_PACKAGE_EN.md`, `A2_PRASC_SOM_VALIDATION_PROTOCOL.md` | Immutable theory inputs, not acceptance authorities |
| Physical core | `research/delegated/a2_physics/physics.py`, `checks.json` | 8 bounded checks, scalar 2D volume-integral model |
| E1/E2/E3/E5 replication | `research/delegated/a2_theory_validation/replication/summary.md`, `results/` | 50 tests and fresh physical/synthetic runs; E3 MC 0/20 evaluable |
| Parent passive bound | `a2_research/results/passivity_checks.json` | 12 positive-loss enclosures, six lossless plus extreme-contrast refusal |
| Residual enrichment | `a2_research/results/enrichment_probe.json` | 16 reference probes; no nonlinear pose-recovery claim |
| Frozen E4 | `research/delegated/a2_solver_audited/final_200.jsonl`, `final_800.jsonl` | 1,200 records each; zero accepted reduced moves |
| E4 tuning/freeze | `research/delegated/a2_solver_audited/tuning_corrected.jsonl`, `frozen_parent.json` | 200 tuning records, seeds 1--10; test seeds 1001--1020 |
| Larger map | `research/delegated/a2_highdim/records.jsonl`, `frozen.json`, `summary.json` | 48 exploratory test runs, six mismatch, separate tuning/smoke |
| Literature | `research/delegated/a2_literature/EVIDENCE.md`, `references.bib` | Bounded, tiered evidence search; not exhaustive novelty proof |
| Public measurement | `research/delegated/a2_literature/data/2001_iop_17_6_301/dielTM_dec8f.exp` | Canonical 2001 data and primary descriptor; no offset truth |
| Corrected measurement pilot | `a2_research/results/measured_corrected_pilot.json` | Object/material model adequacy, NOT array calibration |
| Final physical audit | `a2_research/results/final_parent_checks.json` | Effective Fisher and independently checked propagation/passivity |

In the table, `a2_research/` abbreviates
`research/trispace_self_calibration/a2_research/`.

Frozen E4 config SHA-256:
`f72fd9735e3ff8d07cf44190db46bf60ed7ab322f40eade7bcef1b02bd821779`.
The verifier compares all frozen module hashes to the saved configuration.
The freeze checks solver modules, not every transitive external dependency;
the delivery manifest also hashes the shared physical core and parent modules.

The original E4 runner commands were `study.py tune`, `study.py final --budget
200`, and `study.py final --budget 800` under the audited solver directory, with
the interpreter above. Existing files are protected against overwrite/duplicate
append. For an independent exact rerun, use a separate copy with the same
dependency paths and **new output paths**, retaining archived records and freeze;
do not delete current results merely to satisfy a rerun guard. Do not re-tune on
the final scenes. The legacy `run_e4.py` scene-preparation entry is not supported;
see the known unexercised issues in `FINAL_REVIEW.md`.

## Parent additions / figure regeneration

```sh
"$A2_PY" research/trispace_self_calibration/a2_research/test_passivity.py
"$A2_PY" research/trispace_self_calibration/a2_research/enrichment_probe.py
"$A2_PY" research/trispace_self_calibration/a2_research/final_checks.py
"$A2_PY" research/trispace_self_calibration/a2_research/measured_corrected_pilot.py
"$A2_PY" research/trispace_self_calibration/a2_research/paper_figures.py
"$A2_PY" research/trispace_self_calibration/a2_research/build_paper.py
```

The measured pilot is post-hoc exploratory, with four fixed starts; the primary
descriptor specifies exp(+i omega t), independently cross-checked using training
incident data, and training scattering data set
per-frequency gains. Test gains are never re-estimated. Original worker fitting
results and their physics self-check verdict are withdrawn, not silently fixed.

The editable English Markdown is the manuscript source. `build_paper.py` uses
Tectonic and the macOS Times New Roman fonts to create the review-format TeX/PDF
under `output/pdf/`. Porting the build requires selecting equivalent installed
fonts and rerendering; absolute-path warnings are environmental, not proof of
bitwise portable builds. `build_manifest.json` maps the PDF to its Markdown hash.
The PDF was visually checked after rendering; this does not substitute for
scientific review.

## Automated research provenance

See `PIPELINE_AUDIT.md` and `pipeline_invocation.json`. Installed AI Scientist
development/evaluation was actually invoked; later autonomous revisions did
not provide a validated final paper. Parent-controlled replication, solver
auditing, literature assessment and final interpretation were required.
