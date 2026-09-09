# Reproduction log

Workdir for all original experiment files:
`/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som`

Scratch dir used for regenerated outputs:
`/tmp/trispace_repro.RTqw7v` (created with `mktemp -d`).

Python used: experiment-local `.venv/bin/python`
Python 3.13.13, numpy 2.5.2, scipy 1.18.1, matplotlib 3.11.1, Pillow 12.3.0.
No global/system Python packages were relied on.

## 1. Consolidated-summary generator rerun

```shell
REPRO=/tmp/trispace_repro.RTqw7v
SRC="<experiment dir above>"
mkdir -p "$REPRO"
cp "$SRC"/make_consolidated_summary.py \
   "$SRC"/results_e10_settings_sweep.json \
   "$SRC"/results_e10b_m12_m16_highrank.json \
   "$SRC"/results_e6b_multitx_robust.json \
   "$SRC"/results_e11_vp_state_null.json \
   "$SRC"/results_e5_final.json \
   "$SRC"/results_e9_grid_scenes.json \
   "$SRC"/results_e6.json \
   "$SRC"/results_e7b.json "$REPRO"/
cd "$REPRO"
"$SRC/.venv/bin/python" make_consolidated_summary.py >stdout.txt 2>stderr.txt
echo exit=$?           # 0
cmp consolidated_round3_summary.json "$SRC/consolidated_round3_summary.json"
echo cmp_summary=$?    # 0
cmp consolidated_round3_tables.md "$SRC/consolidated_round3_tables.md"
echo cmp_md=$?         # 0
```

Regenerated hashes (identical to original):

| File | SHA-256 |
|---|---|
| `consolidated_round3_summary.json` | `bcf00325da4be3c2d852ee1674c034a0cb5f941a27bc607b80879ca5b692e657` |
| `consolidated_round3_tables.md` | `a4da40c65d8ba5534dd89ac29e10ee5ee1db5ea5c0f6a0f8beca9d977a95f430` |

Exit code 0; `cmp` exit code 0 for both files.

## 2. E1/E3 core rerun (receiver/transmitter/co-moving derivative check)

```shell
cp "$SRC/geom_som_core.py" "$REPRO/"
cd "$REPRO"
MPLCONFIGDIR="$REPRO/.mplconfig" "$SRC/.venv/bin/python" geom_som_core.py \
  >stdout_core.txt 2>stderr_core.txt
echo exit=$?          # 0
cmp results_e1_e3.json "$SRC/results_e1_e3.json"
echo cmp_json=$?      # 0
```

Representative regenerated FD errors (same as stored):

- rx: receiver 2.241e-7, transmitter(data) 2.385e-7, co-moving 2.230e-7 at eps=1e-4
- theta: receiver 1.524e-8, transmitter(data) 1.536e-8, co-moving 1.358e-8 at eps=1e-4
- Each decade in eps reduces error by ~1e2, i.e. centered-difference O(eps^2).

Regenerated hash (identical to original):

| File | SHA-256 |
|---|---|
| `results_e1_e3.json` | `f66ad73c6c35ae12ae8a1b14bb4a543e29829e66357c108927148fd515b29142` |

Exit code 0; `cmp` exit code 0.

## 3. Independent lift/visible-subspace recomputation

Imported `geom_som_core.py`, `run_e2_e4.py`, `run_e5.py`,
`run_e5_final.py`, `run_e10_settings_sweep.py` from the scratch copy and
recomputed, for every E10 setting at `(alpha_init,p_init)`:

- unrestricted and retained r4/r6 lift residual ratios;
- QH R normalized orthogonality;
- `visible_pose_subspace` singular values, n_vis, hidden_rank.

Output: `max recomputed-vs-stored residual ratio diff=0.000e+00`,
`QH_R diff=0.000e+00`, `subspace diff=0.000e+00`.

Dimension-count check for full circles (independently recomputed):

| Setting | r | Hn columns | orthocomplement dimension | trace(P_perp) | B_red sv1 | stored n_vis/hidden |
|---|---|---|---|---|---|---|
| M12 | 10 | 23 | 1 | 1 | 2.146e-4 | 1/2 |
| M12 | 11 | 25 (rank 24) | 0 | ~0 | 6.16e-19 | 3/0 (`vanished_Bred`) |
| M16 | 14 | 31 | 1 | 1 | 2.643e-5 | 1/2 |
| M16 | 15 | 33 (rank 32) | 0 | ~0 | 1.52e-18 | 3/0 (`vanished_Bred`) |

The same linearized sweep was additionally recomputed for the limited
apertures: M12 half r10 -> n_vis=1/hidden_rank=2 (sv1 8.945e-6); M16
three-quarter r14 -> n_vis=1/hidden_rank=2 (sv1 1.028e-5). The report only
claims the uniform full-circle nonlinear runs, which is a fair scope
restriction, but the linear count is not unique to circular symmetry.

## 4. File and figure integrity

- All PNGs referenced by the report and consolidated file list load with
  Pillow and have non-zero size.
- The four figures copied into
  `final_2026-09-04_06-31-18_geometry_lifted_retained_range_som/figures/`
  are SHA-256-identical to the experiment-dir plots.

## 5. Full-runner commands (not rerun)

The full sweep scripts take no arguments and write into their own directory,
so the expected regeneration commands follow the pattern:

```shell
"$SRC/.venv/bin/python" "$SRC/run_e10_settings_sweep.py" >run_e10_output.txt
"$SRC/.venv/bin/python" "$SRC/run_e10b_m12_m16_highrank.py" >run_e10b_output.txt
"$SRC/.venv/bin/python" "$SRC/run_e6b_multitx_robust.py" >run_e6b_output.txt
"$SRC/.venv/bin/python" "$SRC/run_e11_vp_state_null.py" >run_e11_output.txt
"$SRC/.venv/bin/python" "$SRC/run_e5_final.py" >run_e5_final_output.txt
```

Per audit scope, the 240-case E6b sweep and other full sweeps were **not**
rerun. Command logs are not embedded in the recorded `run_*_output.txt`
files, so exact shell invocation is inferred rather than evidenced.

No requirements file or package freeze exists beside the in-experiment
`.venv`; package versions were read from that interpreter for this log.

