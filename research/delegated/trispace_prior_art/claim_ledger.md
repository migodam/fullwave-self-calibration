# Claim ledger: TriSpace SOM self-calibration anchors

Date: 2026-09-04. Purpose: classify candidate claims that could appear in a
paper against the 16 verified anchors. This is evidence cleanup; it does not
decide final novelty.

Classes:

- **Antecedent/generic** — established prior art; claiming it as new would be
  wrong.
- **Supported recombination** — the evidence supports a conditional,
  model-specific synthesis of prior ingredients, phrased with
  retrieval-bounded language; not an independent new ingredient.
- **Genuinely unresolved (retrieval-bounded)** — no anchor or retrieved text
  contains the exact restricted statement; this only means "not found in this
  bounded set", never "first" or "absent from all literature".
- **Prohibited overclaim** — vacuous, conflating, or globally absolute
  phrasing that the evidence forbids.

| Candidate claim | Class | Basis | Consequence / required phrasing |
|---|---|---|---|
| "Self-calibration of an imaging/measurement system is new." | Antecedent/generic | Hanabusa #3, Cathers #4, autofocus #10-#12, blind calibration #9; source C itself flags array/radar/bilinear self-calibration as established. | Do not award novelty to self-calibration. Say "we formulate/test a self-calibrating ... scheme under these restrictions". |
| "Joint contrast/image reconstruction with unknown transmitter or receiver geometry is new." | Antecedent/generic | Karthik–Ghosh #8 (contrast + transmitter localization), receiver extension #5-#7, radar autofocus #10-#12. | No "first joint imaging and localization" claim. The joint-unknowns problem is prior art. |
| "Source or receiver extension (adding artificial source/receiver DOFs) is new." | Antecedent/generic | Huang–Nammour–Symes #5; Metivier–Brossier #6; da Silva et al. #7. | Credit extension/relocalization; describe only the difference in forward model and in the residual semantics. |
| "Geometry error can be represented by an equivalent current, hence a pose subspace exists." | Prohibited overclaim (as stated) | With an unrestricted full-row-rank sensing map `G_S`, every data perturbation is representable at a single snapshot; equivalence is then vacuous (noted in the TriSpace math-sanity review and source-conflict rules). | Must state the restriction: inside a retained SOM/TSOM current subspace, across stacked/shared pose DOFs, with closed-range assumptions, and with the residual kept irreducible. |
| "The retained/irreducible data residual plus state-equation inconsistency of a canonical SOM-restricted pose-lift jointly drive self-calibration." | Genuinely unresolved (retrieval-bounded) | None of #1-#16 defines this exact restricted object; #14 shows joint use of data/state residuals in CSI but for fixed geometry and as a cost function; #5-#7 use extension residuals without a contrast-source current lift. | Say "within the retrieved set, we did not find ..."; provide the restriction and residual semantics; let Codex decide novelty. |
| "SOM/TSOM current-space subspace decomposition is new." | Antecedent/generic | Chen #1, Zhong–Chen #2; related SOM/TSOM variants in the prior ScholarQA bundle. | SOM/TSOM current-space reduction is prior art; reuse it as machinery. |
| "Current-space SOM modes equal map-observability modes or pose-retention modes." | Prohibited overclaim | Source A and the semantic-conflict digest: `V_S`, `V_D` are current-space objects; map/pose tangent geometry lives in data space (`Ran A`, `Ran B`, `K_eff`); no intrinsic Grassmann comparison without a pullback. | Keep every object labeled by its ambient space (current/map/data/pose); never identify spectra across spaces. |
| "Jointly using the data equation and the state equation is new in contrast-source inversion." | Antecedent/generic | CC-CSI #14 (cross-correlated data/state error) and the broader CSI family. | Only the calibration-specific role of the two residuals (not their joint use in inversion) can be a candidate contribution, and only if no retrieved prior does it. |
| "Adding virtual antenna samples improves SOM imaging; virtual measurements are new." | Antecedent/generic | Zhang et al. #13 synthesizes virtual antennas and uses SOM. | Aperture/data augmentation with virtual antennas is prior art and is not unknown-geometry self-calibration. |
| "Virtual experiments / equivalent contrast-source constructions are new." | Antecedent/generic | Bevacqua et al. #15 and the equivalent-source literature more broadly. | Distinguish a pose-equivalent current lift for unknown `G_S` from solver-level virtual-experiment recombination. |
| "Limited aperture / asymmetric operator factorization explains our geometry-dependence result." | Supported recombination (only if used carefully) | Audibert–Haddar #16 analyzes non-symmetric factorization of the far-field operator under limited aperture. | It is a qualitative, fixed-geometry result; cite only as factorization/aperture-structure context, not as calibration evidence. |
| "Bilinear identifiability up to gauge applies to our joint unknowns, so generic recovery is guaranteed." | Prohibited overclaim (without model-specific check) | Li–Lee–Bresler #9 gives general bilinear identifiability conditions; they are not a Helmholtz/full-wave specialization. | Claim only that the general theory governs the bilinear structure; verify rank/identifiability conditions for the actual operators and state which gauge remains. |
| "Microwave/FWI data calibration with 2-port error models or DNN conversion removes all model error." | Prohibited overclaim (if generalized) | Cathers #4 and Hanabusa #3 calibrate specific systematic/measurement errors under known geometry. | Keep hardware calibration and pose/geometry uncertainty separate; neither demonstrates removal of full-wave pose-induced operator error. |
| "We scan the literature exhaustively; no prior work exists." | Prohibited overclaim | All screening was bounded; Semantic Scholar/ScholarQA transport had DNS/429 failures; only targeted lookups were available. | Never use "first", "no prior work", or "absent from all literature"; state the retrieval boundary in the paper. |
