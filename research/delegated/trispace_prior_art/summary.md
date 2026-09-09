# TriSpace SOM self-calibration — bounded prior-art cleanup summary

Date: 2026-09-04. Inputs: the three user source documents via the local
digest (`research/delegated/trispace_source_digest/`), the targeted prior-art
addendum, the prior ScholarQA verification bundle, plus targeted exact-title
abstract lookups for anchors lacking local descriptions. No broad re-search
was run. All novelty statements are retrieval-bounded.

## Findings

1. The 16 anchors cleanly split into two groups: current-space SOM/TSOM
   machinery with fixed known geometry (#1, #2, #13, #14, #15) and
   unknown-geometry/calibration/extension work in other representations
   (#3-#12). No anchor bridges the two by making SOM/TSOM spectral subspaces
   themselves pose-dependent.
2. Generic self-calibration, joint image/contrast + geometry/phase
   estimation, source/receiver extension, virtual experiments, and
   data+state residual combination in contrast-source inversion are all
   confirmed prior art; none of these can be a novelty claim.
3. The claimed narrow object — a canonical pose-equivalent induced-current
   lift inside a retained SOM/TSOM current space, with an irreducible data
   residual plus a state-equation inconsistency used jointly for
   self-calibration — was not found in the verified set (retrieval-bounded).
4. A bare "geometry error can be written as a current" statement is vacuous
   when the unrestricted sensing map is full row rank (single-snapshot
   representability); any defensible formulation must state the
   retained-subspace or stacked/shared-pose restriction.
5. Receiver-extension FWI (#5-#7) is the closest genus to geometry DOFs
   absorbing data misfit, but it operates on seismic wave-equation
   source/receiver coordinates, with no contrast-source current space and no
   residual-pair construction.
6. The evidence for almost all anchors is abstract/landing level; no anchor
   was full-text inspected in this cleanup, and the earlier ScholarQA/Semantic
   Scholar broad search failed on DNS/429, so several retrieval gaps remain
   (e.g., full text of Karthik-Ghosh, Huang et al., Metivier-Brossier).
7. Current-space vs map/pose data-space objects must stay labeled; the
   digest confirms `V_S/V_D` are current-space and `Ran A/Ran B/K_eff` are
   data-space, so overlap claims must not equate SOM modes with map
   observability or pose-retention modes.

## Artifacts

- `research/delegated/trispace_prior_art/screened_papers.md`
- `research/delegated/trispace_prior_art/claim_ledger.md`
- `research/delegated/trispace_prior_art/verified_references.json`
- `research/delegated/trispace_prior_art/novelty_overlap.md`
- this file, `research/delegated/trispace_prior_art/summary.md`

Pre-existing files in the folder (`TASK.md`, `SEARCH_PLAN.md`,
`transport_error_run1.json`, `transport_run1.err`) were left untouched; the
transport record shows the earlier ScholarQA broad search failed without
candidates.

## Limitations

- Evidence tier for every row is stated in `screened_papers.md`; most rows
  are `abstract` or `local-source`, none is `full-text` in this run.
- Metadata fields were not guessed; `verified_references.json` omits any
  volume/issue/page field that was not verified.
- Targeted lookups on 2026-09-04 verified identity/abstracts only; they do
  not establish how each method handles multiple scattering, convergence, or
  identifiability in full text.
- Semantic Scholar/ScholarQA transport was unavailable (DNS error earlier;
  anonymous 429 in prior runs), so "not found" reflects a bounded retrieval
  set, not a literature-wide absence.

## Paths suggested for the main thread

1. Fix the restriction language before any draft: retained subspace of which
   operator (nominal `G_S(X0)` vs pose-dependent `G_S(X)`), how the
   irreducible residual is computed, and which state residual is measured.
2. Decide and state whether per-channel phase/gain calibration and antenna
   phase centers are inside the pose variable; #3/#4 show a separate
   measurement-calibration track.
3. If the paper keeps a self-calibration claim, phrase it as "we formulate
   and test a SOM/TSOM-restricted pose-lift self-calibration ... " and cite
   #5-#7 plus #8 for the joint estimation framing.
4. For final wording, distinguish the two degenerate cases: full-row-rank
   single snapshot (equivalence trivial) and non-closed range in continuum
   sensing (asymptotic-only compensation).

## Where Codex judgment is needed

- Whether the restricted lift + residual-pair construction is a supported
  recombination or a genuinely new narrow formulation (not decided here).
- Whether "irreducible data residual" and "state-equation inconsistency" are
  the right paired statistics for pose coordinates vs per-channel
  calibration, and whether a retained-subspace leakage residual alone is
  enough to avoid the full-row-rank vacuity.
- Whether the current-space `V_P` definition problem (open in the source
  digest) is solvable before any overlap claim can be finalized.
- Final claim language must remain retrieval-bounded and must not say
  "first", "no prior work", or imply exhaustive coverage.
