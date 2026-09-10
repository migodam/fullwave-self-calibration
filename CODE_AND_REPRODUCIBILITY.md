# Code and reproducibility guide

This is a versioned research archive, not one application with a global install
command. Original source is byte-preserved. Do not run every historical script
indiscriminately: some generate output, call separately configured providers,
or require omitted third-party solvers / measured datasets.

## Q5: current independent development code

[Reproduction scope](research/q5_remaining_v1/Q5_REPRODUCIBILITY.md) and
[complete code/parameter map](Q5_CHATGPT_START_HERE.md) are the current entry.
Tested numerical package versions: NumPy 2.5.3, SciPy 1.16.3, Treams 0.4.7,
Python 3.13.13. These record the local environment, not a clean-install test.
The code imports the included A3 `maxwell3d.py`; no paid model is needed for
the numerical scripts. `run_native_checks.py` is a separate optional provider
entrypoint requiring a private profile and the installed external pipeline.
Private profiles, credentials, virtual environments and mpmath cache are not
published. Never run provider entrypoints merely to inspect this archive.

## A4: self-contained numerical source package

- [Requirements](research/a4_reliability_v1/requirements.txt) and
  [README](research/a4_reliability_v1/README.md).
- [Modal fields and spectra](research/a4_reliability_v1/src/modal.py).
- [Vector DDA](research/a4_reliability_v1/src/dda.py).
- [Task risk](research/a4_reliability_v1/src/risk.py),
  [coverage](research/a4_reliability_v1/src/coverage.py),
  [experiments](research/a4_reliability_v1/src/experiments.py).
- [Regression tests](research/a4_reliability_v1/tests/test_core.py),
  [test launcher](research/a4_reliability_v1/run_tests.sh).
- [Read-only stored-evidence audit](research/a4_reliability_v1/audit_package.py).
- [Fresh-run wrapper](research/a4_reliability_v1/rerun_fresh.py): use a new output
  directory; never overwrite frozen records.

Example isolated installation after cloning (requires package downloads):

```sh
cd research/a4_reliability_v1
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
sh run_tests.sh
python audit_package.py
```

The original package reports 21 passing tests. During this publication pass,
`audit_package.py` passed all stored-data/hash invariants using the existing A3
numerical environment. A pytest attempt could not start because that environment
lacks pytest; **the 21 tests were not rerun in this publication pass**. No fresh
full experiment run or theorem review is implied. See the supplied corrections:
one known historical reference-control bug is retained with a separate repaired
supplement rather than silently changing frozen source.

## A3: full-wave solvers, calibration and adjudication

Start with the original [reproducibility record](research/trispace_self_calibration/a3_research/REPRODUCIBILITY.md).
Its local paths and original environment are historical provenance, not a claim
that every dependency is bundled here.

| Purpose | Code |
| --- | --- |
| Maxwell reference and tests | [maxwell3d.py](research/trispace_self_calibration/a3_research/maxwell3d.py), [test_maxwell3d.py](research/trispace_self_calibration/a3_research/test_maxwell3d.py) |
| Non-spherical models | [nonspherical3d.py](research/trispace_self_calibration/a3_research/nonspherical3d.py), [calibration](research/trispace_self_calibration/a3_research/nonspherical_calibration.py) |
| Frozen SOM comparison | [rom_final.py](research/trispace_self_calibration/a3_research/rom_final.py), [protocol](research/trispace_self_calibration/a3_research/ROM_FINAL_PROTOCOL.md), [raw records](research/trispace_self_calibration/a3_research/results/rom_final.json) |
| Model-discrepancy stress | [shape_material_stress.py](research/trispace_self_calibration/a3_research/shape_material_stress.py), [two-stage fits](research/trispace_self_calibration/a3_research/calibrate_then_image.py) |
| Bias / risk checks | [frequency information](research/trispace_self_calibration/a3_research/check_frequency_information.py), [A3_3 risk](research/trispace_self_calibration/a3_research/check_a3_3_task_risk.py) |
| Missing phase branches | [branch_repair.py](research/trispace_self_calibration/a3_research/branch_repair.py) |
| Supporting solvers | [research/delegated](research/delegated) — A2/A3 imports retain their original relative structure |

ADDA external source/build, Treams installation and original measured-data files
are not bundled. ADDA-based reruns need a reviewed external setup; tiny linear
checks do not establish successful full Maxwell reproduction. Some provenance
audits intentionally reference original absolute paths and will need a separate,
explicit relocation layer rather than editing frozen hashes.

## A2 / A1 history

- [A2 reproducibility](research/trispace_self_calibration/a2_research/REPRODUCIBILITY.md),
  [audited solver](research/delegated/a2_solver_audited),
  [physical forward model](research/delegated/a2_physics).
- [A2 autonomous source](experiments/idea_loops/loop_2026-09-05_17-55-24_a2_development/experiment_a2_prasc_theory_validation).
- [A1 physical source and audits](experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry).
- [TriSpace equivalent-current validation](experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som).

Read adjacent parent audits and invalidation notices before choosing a historical
implementation. Old and corrected code are both retained to explain the process.
