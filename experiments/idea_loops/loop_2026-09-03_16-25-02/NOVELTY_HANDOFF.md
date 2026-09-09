# Novelty-stage handoff

Date: 2026-09-03

## Automated-stage outcome

The automatic deep-novelty verifier did not return a novelty verdict. Its
literature-tool calls repeatedly received HTTP 429 responses and the agent
raised `MaxTurnsExceeded: Max turns (12) exceeded`. This is an infrastructure
failure, not evidence for either `is_novel=true` or `is_novel=false`.
`read_paper_in_depth` was also degraded because the helper rejected the model
name `deepseek-v4-pro`.

## Codex scientific gate

Development is conditionally allowed only for the narrowed contribution in
`research/literature/TARGETED_PRIOR_ART_ADDENDUM.md`: the combination of
whitened/realified map and physical-pose tangents for nonlinear full-wave
volumetric contrast, the finite-rank Schur information defect, the no-prior
principal-angle retention law, the Born-empty-background second-order split,
and executable frequency/trajectory/rank-event falsifiers.

The following are established prior-art territories and must not be claimed as
novel: joint inverse-scattering reconstruction and transmitter localization,
joint image/motion-error estimation, Fisher information for SLAM, generic
bilinear identifiability, or SOM/current-subspace decomposition. In particular,
Karthik and Ghosh (PIERS 2023, DOI
10.1109/PIERS59004.2023.10221374) explicitly reconstruct contrast while
localizing transmitters.

No use of “first”, “no prior work”, or any global/exhaustive novelty wording is
permitted. The manuscript must report this conditional gate and retrieval
limitations.
