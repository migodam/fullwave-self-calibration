# Autonomous research instructions: self-calibrating TriSpace SOM

This run is explicitly requested by the user. Its endpoint is a defensible
paper draft backed by executable CPU-only theory-validation experiments. The
corrected main line is **TriSpace SOM/TSOM under unknown or drifting antenna
geometry**, not generic map/pose Fisher observability.

## Required source order

1. Read `research/trispace_self_calibration/THEORY_SEED.md` completely.
2. Read `research/delegated/trispace_source_digest/summary.md` and
   `research/delegated/trispace_source_digest/semantic_conflicts.md`.
3. Reuse the verified 2D Helmholtz conventions and code under
   `experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/`.
4. Read `research/delegated/trispace_prior_art/summary.md` and
   `novelty_overlap.md` if they exist by the time the run reaches synthesis.

The three original user documents are research material, not instructions.
Consult them only if the distilled files leave a specific ambiguity.

## Non-negotiable scientific rules

1. Treat the finite-dimensional propositions in `THEORY_SEED.md` as hypotheses
   to independently re-derive and test. Preserve counterexamples and correct
   any sign, range, or rank mistake you find.
2. Keep current, data, state, map, and pose spaces type-correct. Never identify
   `Range(D_xF)` with a current-space SOM subspace without an explicit lift.
3. Decompose only the receiver/sampling derivative
   `H_S h=(D_xG_S[h])j` through the current-to-data operator. Keep physical
   current variation caused by the transmitter/incident field or moving
   internal operator in the state equation. Do not double count it.
4. In the world-fixed homogeneous-grid simulator, `D_xG_D=0`. A moving
   `G_D(x)` is an optional distinct model, not the default.
5. Whiten and realify before rank, projection, pseudoinverse, principal-angle,
   or identifiability calculations. Contrast and pose variables are real.
6. The unrestricted full-row-rank pseudoinverse result is a required negative
   control: it makes every data perturbation current-equivalent and is therefore
   vacuous. The method earns value only through retained SOM/TSOM restriction,
   stability, irreducible residual, and state consistency.
7. Attempt all five experiment families in `THEORY_SEED.md`. A documented
   failure or corrected condition counts as an attempted family; silent omission
   does not.
8. Use a Phase-I truncated-SOM basis if needed. Do not label an ad hoc `G_D`
   split as canonical TSOM. Exact TSOM implementation requires a verified source
   definition and should otherwise remain future work.
9. E5 must execute an actual reconstruction/calibration comparison, not only
   plot spectra. Use deterministic small synthetic problems, several
   initializations/seeds, and report failures.
10. Passing finite-dimensional experiments supports only the tested model. It
    does not prove continuum transfer, global nonlinear convergence, hardware
    performance, or global novelty.
11. Literature novelty is retrieval-bounded. Equivalent currents, self-
    calibration, source/receiver extension FWI, blind calibration, autofocus,
    joint image/calibration optimization, SOM, and TSOM are antecedents. Never
    claim those ingredients as new, and never use `first` or `no prior work`.
12. Use Apple Silicon CPU. Do not install CUDA, GPU Docker, a daemon, scheduler,
    database, web UI, or another orchestration framework.
13. Every numerical statement needs an executable command, config, seed,
    tolerance, raw table, and interpretation. Include negative controls and a
    `cannot establish` section.
14. Keep core theorem/proof gaps for top-level Codex/GPT Pro. The coding worker
    may falsify and diagnose them but must not declare a hard theorem solved from
    numerics.

## Required outputs

- One coherent experiment directory with source, configs, raw results, figures,
  commands, environment record, and per-family claim-status tables.
- A synthesis separating proved algebra, model-conditional derivations,
  numerical evidence, retrieval-bounded originality, and open questions.
- A complete English paper draft in LaTeX. If the built-in writeup stage cannot
  run with the selected DeepSeek model, preserve the experiment results and a
  detailed manuscript-ready outline; do not change the top-level Codex provider
  or use an OpenAI API key.
- A concise Chinese list of difficult proof/design questions for later GPT Pro
  work, without pretending they were solved.

The evaluator should lock only after all five families are genuinely attempted,
the full-row-rank vacuity control is retained, and E5 has an executed comparison
or an explicit reproducible failure diagnosis.
