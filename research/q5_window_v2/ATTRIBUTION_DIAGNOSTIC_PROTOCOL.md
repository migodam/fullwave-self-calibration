# Post-hoc attribution diagnostic, frozen before this diagnostic's execution

This diagnostic is motivated by the already inspected V2 development outcomes.
It is not confirmatory evidence, a new method, or a new final test stream.
Use the unchanged observations_v2.npz, main L=3 solver, fixed three material/pose
starts and original 80-evaluation/tolerance settings. No acquisition is selected.

Fit the materials under: true complex gain fixed; true receiver shift fixed;
both true gain and shift fixed. No reference sample enters these oracle objectives.
All other nuisance parameters retain the original annulus/domain. Record every
start and all failures. Oracle results diagnose possible nuisance compensation;
they are not deployable baselines or a proof of continuum causal attribution.
