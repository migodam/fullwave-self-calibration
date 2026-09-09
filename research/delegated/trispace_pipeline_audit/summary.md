# Trispace pipeline audit - mechanical reproducibility and code audit

Scope: `experiments/idea_loops/loop_2026-09-04_02-58-31/` (workdir
`experiment_geometry_lifted_trispace_som/`, final report
`final_2026-09-04_06-31-18_geometry_lifted_retained_range_som/experiment_report.json`).
This is a bounded, non-mutating audit; the experiment directory was not modified.
No final scientific or novelty judgment is made here.

## Verdict

**Mixed pass with reporting overclaims.** The artifacts are internally
consistent at the table level, the consolidated summary generator is
deterministic and byte-reproducible, the E1/E3 core check is
byte-reproducible, and the code contains no evidence of ground-truth leakage,
accidental result reuse, or silent data-loss in the audited solvers. The
lift formulas, realification, rank thresholds, projectors, orthogonality
checks, and derivative checks are implemented as described and independently
recomputed exactly.

However, several report sentences say "exact" or "full precision" where the
raw JSON only supports rounded agreement (or no agreement):

- "reduced_r4 matches direct to full precision in every setting" is false at
  the half-circle setting (direct median pose error 0.010429502, reduced_r4
  0.010209086) and only rounded-agreement in several others.
- "reduced_r6 also matches direct everywhere pose is visible at r6" is
  contradicted by the half-circle setting: the linear r6 subspace says
  `n_vis=3, hidden_rank=0`, yet one reduced_r6 seed lands in a spurious basin
  with pose error 1.2867 and median differs from direct.
- "Trailing r=M-1 and r=M ... reduced matches direct exactly" goes beyond the
  runs: only r=M-1 was executed, and M=12 r11 differs from direct at
  parameter/error level (~1e-4 / ~3e-5) though the 4-significant-figure table
  rounds them together.

The more important interpretive findings are that the "visible" subspace is a
fixed linearization construct and the r=M-2 "hidden" transition is a
dimension-counting bound rather than an independently discovered physical
phase transition, and that convergence flags do not detect the recorded
local-minimum cases.

## Key verified findings

1. Table/summary agreement: rerunning
   `make_consolidated_summary.py` on the stored JSONs reproduced
   `consolidated_round3_summary.json` and `consolidated_round3_tables.md`
   **byte-for-byte** (SHA-256 `bcf00325...` and `a4da40c6...`).
2. Core reproduction: rerunning `geom_som_core.py` reproduced
   `results_e1_e3.json` **byte-for-byte** (SHA-256 `f66ad73c...`), including
   receiver-only, transmitter-only, and co-moving FD checks at
   O(eps^2) error levels.
3. Independent recomputation of every E10 retained/unrestricted lift ratio,
   QH R orthogonality ratio, and E10/E10b visible-subspace SVD matched the
   stored values to zero difference.
4. All lifted/derivative formula mechanics are correct as coded:
   `Q=G_s V_r`, `c=pinv(Q)b`, `R=Qc-b`, realified column blocks, relative
   rank thresholds, `P_perp`, and V_vis orthonormality.
5. Baselines are genuinely distinct, and truth appears only in
   observation generation, error evaluation, and the named
   `known_pose`/`known_alpha` oracle baselines, not in non-oracle solver
   residuals or Jacobians.
6. "Reduced matches direct" in n_vis=3 regimes is a reparameterization of the
   same full-physics objective; the observed match is implementation/basin
   consistency, not a nontrivial recovery proof. The half-circle counterexample
   shows full-rank V_vis does not even guarantee matching basin selection.
7. The r=M-2 transition uses a fixed nominal basis at
   `(alpha_init,p_init)`, has a robust ~1e-15 singular-value gap at the
   transition, and is a linear-algebra/dimension-count effect: Hn has
   2M x (2r+3) real columns, so at r=M-2 its orthogonal complement has
   dimension 1 and `rank(B_red) <= 1` (hence hidden_rank=2). It is not a
   threshold artifact at r=M-2; the threshold artifact is the
   `hidden_rank=0` reclassification in the empty-lift r>=M-1 state, which the
   code separately flags via `vanished_Bred`.
8. Convergence flags can and do hide bad local minima. In the half-circle
   setting all three `known_alpha` seeds converge (status=1,
   optimality<1e-7) to pose error ~0.86 with residual ratio ~0.27 (noise
   ratio 0.0316), and one `reduced_r6` seed converges to pose error 1.2867
   and residual 0.0396. The report records the caveat, but the consolidated
   tables omit final residual ratios, so the tables alone cannot reveal it.
9. The E11 hard/soft VP and E7/E7b state-consistency negative results are
   correctly framed as negative in the report. The E11 data residual for VP
   rows is a VP-model residual, while direct/reduced rows use full physics;
   the JSON and interpretation text explain this, but the consolidated table
   does not carry the caveat.
10. Reproducibility is good but not fully self-describing: runners use fixed
    seeds, no global random state, and write-and-assert their own JSON/PNG;
    however there is no requirements manifest or recorded command list, and
    E6/E6b "success" is effectively `status>0` with no optimality gate and
    without storing `status` in the JSON.

## Key risks

- Report wording (key results 2 and 5, summary sentence) overstates exactness
  and the r=M trailing-state claim relative to executed runs and raw values.
- The "hidden direction"/"visible subspace" language can be read as a global
  physical recovery claim; it is a local linearized fixed-basis construct.
- "All 240/240 runs converged" is based on a loose success definition
  (`status>0`), not an optimality/residual gate.
- Objective/residual columns are missing from the main consolidated tables,
  which is exactly the information needed to distinguish the recorded local
  minima from good solutions.

