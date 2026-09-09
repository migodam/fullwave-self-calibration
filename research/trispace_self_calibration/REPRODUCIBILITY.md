# Reproducibility manifest

Date: 2026-09-04

This manifest records the executable evidence used by the paper draft. The
locked autonomous-research directory is not modified by the parent correction.

## Environment observed in the experiment virtual environment

| component | version |
|---|---:|
| Python | 3.13.13 |
| NumPy | 2.5.2 |
| SciPy | 1.18.1 |
| Matplotlib | 3.11.1 |
| Pillow | 12.3.0 |

Platform policy: Apple Silicon CPU. No CUDA or GPU infrastructure is required.

## Paths

```text
EXP=experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som
FINAL=experiments/idea_loops/loop_2026-09-04_02-58-31/final_2026-09-04_06-31-18_geometry_lifted_retained_range_som
AUDIT=research/trispace_self_calibration/validation
PY=$EXP/.venv/bin/python
```

These lines are notation only; commands below use explicit paths and do not
depend on shell variables.

## Minimal verified reruns

From the repository root:

```bash
./experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python \
  experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/geom_som_core.py

./experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python \
  experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/make_consolidated_summary.py

./experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python \
  research/trispace_self_calibration/validation/validate_rank_saturation.py
```

The first two were independently reproduced in the mechanical audit. The
third was rerun twice by the parent after manuscript synthesis.

## Full stored experiment entry points

| scientific block | entry point | primary output |
|---|---|---|
| E1/E3 lift and derivatives | `geom_som_core.py` | `results_e1_e3.json` |
| E2/E4 coordinate drift and hiding | `run_e2_e4.py` | `results_e2_e4.json` |
| gauge control | `run_e2b_gauge.py` | `results_e2b_gauge.json` |
| E5 reconstruction | `run_e5.py`, `run_e5_final.py`, `run_e5_seeds.py` | `results_e5*.json` |
| multi-transmitter | `run_e6_multitx.py` | `results_e6.json` |
| robust/directional source | `run_e6b_multitx_robust.py` | `results_e6b_multitx_robust.json` |
| state penalty | `run_e7_state_consistent.py`, `run_e7b_soft_state.py` | `results_e7*.json` |
| grid/scene sensitivity | `run_e9_grid_scenes.py` | `results_e9_grid_scenes.json` |
| settings sweep | `run_e10_settings_sweep.py` | `results_e10_settings_sweep.json` |
| high-rank sweep | `run_e10b_m12_m16_highrank.py` | `results_e10b_m12_m16_highrank.json` |
| hard/soft VP null | `run_e11_vp_state_null.py` | `results_e11_vp_state_null.json` |
| corrected saturation audit | `validation/validate_rank_saturation.py` | `validation/rank_saturation_results.json` |

Each script is deterministic under its stored seed list and writes its own
JSON/PNG artifacts. The original run logs do not consistently preserve the
exact shell invocation or exit code, and E6/E6b JSON does not retain the SciPy
status code. Those are reproducibility limitations, not silently repaired
history.

## Stable hashes

| artifact | SHA-256 |
|---|---|
| `geom_som_core.py` | `2283d677788b140a82eb67c8e6b47d2b47c9946f7415e714cc086c473e935926` |
| `results_e1_e3.json` | `f66ad73c6c35ae12ae8a1b14bb4a543e29829e66357c108927148fd515b29142` |
| `consolidated_round3_summary.json` | `bcf00325da4be3c2d852ee1674c034a0cb5f941a27bc607b80879ca5b692e657` |
| `consolidated_round3_tables.md` | `a4da40c65d8ba5534dd89ac29e10ee5ee1db5ea5c0f6a0f8beca9d977a95f430` |
| corrected `rank_saturation_results.json` | `e80302e390b5b46b5c6853bc21c0e58a4664dcf947b3f655abbd8cd1a34130ff` |
| corrected `rank_saturation_summary.md` | `5547b9c85a8246351c0de94c6294632b6edc74cb07788fb9ca915e5239f04351` |
| corrected `rank_saturation_corrected.png` | `e5ef59532d159c1f62d2736906cf6c2cd46d2f74888f9c56c2a9637a026facde` |

The three corrected-artifact hashes were identical across two consecutive
parent reruns.

## Acceptance interpretation

- Byte-identical or low-residual identities validate deterministic mechanics.
- A SciPy positive status is not a scientific success criterion.
- The old E10/E10b relative-only high-rank interpretation is not accepted.
- The parent-corrected rank-saturation JSON and figure are the manuscript
  evidence for `r>=M-1`.
- Old E4 “joint data+state” labeling is not accepted as a verification of the
  state half of the hiding condition.
