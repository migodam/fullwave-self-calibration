# Prioritized issues

## P1 - Report overclaims exact equivalence and trailing-rank coverage

`experiment_report.json` key results:

- "reduced_r4 matches direct to full precision in every setting"
- "reduced_r6 also matches direct everywhere pose is visible at r6"
- "Trailing r=M-1 and r=M ... reduced matches direct exactly"

The first is contradicted by the half-circle setting even at rounded precision
(direct median 0.010429502 vs reduced_r4 median 0.010209086) and is not
machine-precision elsewhere (e.g. M12_k16 full-circle max pose-error
difference ~1.6e-7). The second is contradicted by the half-circle r6 row:
stored subspace is n_vis=3/hidden_rank=0, yet direct median 0.010429502 vs
reduced_r6 median 0.027149679, with one seed at pose 1.286716. The third goes
beyond executed runs: only r=M-1 (not r=M) was run, and M12 r11 differs from
direct by max p difference 1.14e-4 and map-error difference 3.10e-5, which is
hidden by the 4-significant-figure markdown.

Suggested repair: reword to "agrees to the printed precision in settings that
reach the same basin" or "matches direct under full-rank V_vis when both
solvers converge to the same local basin", and do not claim r=M nonlinear
equivalence.

## P1 - Convergence gate hides the recorded local minima; residual columns missing

All half-circle `known_alpha` seeds and one half-circle `reduced_r6` seed pass
`status>0 and optimality<1e-7` while in a different basin:

| setting/method/seed | status | optimality | pose error | final residual ratio |
|---|---|---|---|---|
| half known_alpha seed 0 | 1 | 8.57e-11 | 0.8603 | 0.2663 |
| half known_alpha seed 1 | 1 | 9.66e-11 | 0.8612 | 0.2775 |
| half known_alpha seed 2 | 1 | 8.18e-11 | 0.8569 | 0.2732 |
| half reduced_r6 seed 0 | 1 | 1.58e-11 | 1.2867 | 0.03956 |

Direct/known_pose residuals are ~0.024-0.031 at noise ratio 0.0316. The report
correctly records the caveat, but neither `consolidated_round3_tables.md` nor
the E10 tables show final residual ratios, so a reader cannot detect the
basin failure from the tables.

## P2 - "Hidden/visible" and r=M-2 "transition" need explicit linearization/dimension-count framing

The subspace is computed once at `(alpha_init,p_init)` from SVD of G_s at that
pose; the reduced solver holds V_vis fixed. `hidden_rank=2` at r=M-2 follows
from the realified column count: Hn is 2M x (2r+3), so at r=M-2 its
orthocomplement has dimension 1 and `rank(B_red)<=1`. This is a robust
algebraic bound in the tested cases (sv2/sv3 ~1e-19, sv1 ~2e-4), not a
threshold artifact at the transition. It is also not unique to circular
symmetry: recomputing the linearized sweep for M12 half and M16
three-quarter gives the same r=M-2 count. The report's wording ("clean rank
transition ... hidden_rank=2") is easy to read as a physical/global recovery
claim; the evidence supports a local, fixed-basis, dimension-count statement
plus nonlinear reduced-runs that freeze the complement.

The threshold/numerical artifact is the trailing r>=M-1 state, where
P_perp~0 and the 1e-8-relative rule returns n_vis=3 on ~1e-18 singular
values. That is handled in code by `regime=vanished_Bred`/`VANISH_SV_ABS`,
but the printed tables show only n_vis/hidden_rank, so the transition curve
"n_vis drops to 1 then jumps back to 3" can be misread.

## P2 - Reduced-vs-direct equality is not a nontrivial result

When n_vis=3, `p=p_init+V_vis q` with a 3x3 orthogonal V_vis is an invertible
reparameterization of the same full-physics residual. Equality of solutions
proves that the reduced code evaluates the same forward model/Jacobian and
that both runs reach the same local basin in that setting. It does not prove a
"lossless" property beyond the parameterization, and it does not guarantee
basin selection (half-circle reduced_r6). Any paper use should be framed as a
consistency check, not independent evidence for a new recovery mechanism.

## P2 - E11 table residual normalization is not self-explanatory

`consolidated_round3_tables.md` Table R3-4a shows `data res med` for
`direct`/`reduced_r6` (full-physics forward residual) and for
`vp_hard_r6`/`vp_soft_r6` (VP-model data residual) without a caveat. The
stored E11 JSON/interpretation text distinguishes these and also records
`full_forward_data_residual_ratio` (~0.83 soft), but the table is apples to
oranges without the note. Hard-row medians also include both non-converged
seed-0/1 runs (success 0/2); the report text says this, but the table caption
does not.

## P3 - E6/E6b "converged" success definition is effectively status>0

In `run_e6_multitx.py`/`run_e6b_multitx_robust.py`:

```python
success = bool(result.status > 0 and (result.optimality < 1e-8
                                      or result.status in (1, 2, 3, 4)))
```

Positive scipy statuses are always in {1,2,3,4}, so the optimality clause is
redundant and status=4 xtol stops count as success. E6/E6b JSON does not store
`status`, so this cannot be re-audited from the result files alone. All 240
E6b direct/wrongpose records have low optimality in the raw JSON, so no bad
case was observed, but "all 240/240 runs converged" should be described with
the weaker gate.

## P3 - Reproducibility metadata incomplete

- No requirements.txt/package freeze; only the in-experiment `.venv` pins
  numpy/scipy/matplotlib/Pillow versions.
- `run_*_output.txt` files do not record the exact shell command or exit
  code.
- E5-final Jacobian FD checks print to stdout and are asserted during the
  run but are not stored in `results_e5_final.json`; the log
  `run_e5_final_output.txt` is the only retained evidence.

## P3 - Smaller reporting ambiguities

- R3-1 lift rows do not state that the lift perturbation is only
  `b_pert=B_r[:,0]` (receiver x-translation at alpha_init/p_init); the raw
  JSON note does.
- Report "known_pose map error 0.008-0.012" omits the three-quarter setting
  median 0.006409 if intended to cover all settings.
- Report "third realified SV 2.0e-4" (E6b dipole L1) refers to the
  unnormalized third `realified_sv` 2.0299e-4, while Table R3-3d prints the
  column-normalized third singular value 0.0583; two labels that look like
  the same quantity are not.
- E7/E7b tables count `success` as scipy success + status>0 and do not show
  `n_strict` (optimality<1e-8), so convergence strength is understated.

