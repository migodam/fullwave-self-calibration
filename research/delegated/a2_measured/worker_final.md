Bounded pilot complete. All outputs under `research/delegated/a2_measured/`; no solver/physics/highdim code or raw data was modified.

**Findings**

1. Column semantics resolved from the fulltext intro (Belkebir–Saillard 2001, IP 17:1565–1571): view index, receiver index, frequency in GHz, then Re/Im **total** E and Re/Im **incident** E; scattered = total − incident; `exp(+iωt)`; TM = E ∥ cylinder axis.
2. Published geometry: radius 15 mm, εr = 3±0.3, centre “about 30 mm” off-axis (direction unpublished), emitter 720±3 mm, receiver 760±3 mm, 36 views × 49 receivers × 1–8 GHz.
3. Ingestion QA **MET**: 14,112 rows, all finite, zero duplicates, full frequency coverage, hash matches the literature worker’s `SHA256SUMS.txt`; 4/4 executable checks pass.
4. Measured incident field is rotationally consistent across views (<1%) and phase-confirms `exp(+iωt)`; the data itself is clean.
5. Line-source and full-ring equivalent-source incident models were tried and rejected (≈99% residual / near-field ill-conditioning) — the horn beam pattern dominates; failures preserved.
6. Adopted model: analytic Mie series under plane-wave illumination from the emitter direction, with one complex gain per frequency profiled identically for all variants. Core self-verified: optical theorem 1−1e-15, boundary continuity 7e-16, PEC limit 1.2e-4, series tail 2e-54.
7. Joint-fit identifiability validated by the injected-perturbation pilot: centre (24,−8) mm and εr 2.7 recovered to 8 µm / 3e-4 under 5% noise. Injected errors are synthetic, not natural truth.
8. Real data: recovered centre is stable at ~26.5–27 mm from the axis across band choices — consistent with the published “about 30 mm”; direction is a gauge choice, not claimed physical.
9. εr is **not** identifiable under this model: fitted value drifts 5.2 → 2.5 → 1.2 as the band grows, and the 4-param “loss” (tan δ = 0.33) is an effective parameter absorbing model error, outside the published 3±0.3.
10. Residuals are 2–8% at 1–2 GHz but 0.4–0.9 at 4–8 GHz: plane-wave incident modeling breaks down where the horn beam narrows — model misspecification, not a data defect.

**Gates:** download MET · ingestion QA MET · model self-verification MET · model fit BOUNDED/PARTIAL · actual unknown-geometry recovery **NOT MET** (no ground truth exists; only consistency with “about 30 mm”).

**Artifacts**

- [SUMMARY.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_measured/SUMMARY.md), [METADATA.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_measured/METADATA.md)
- [ingestion.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_measured/ingestion.py), [test_ingestion.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_measured/test_ingestion.py), [model.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_measured/model.py)
- [results.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_measured/results.json) (plus `results_ingestion.json`, `results_model.json`)
- [plots/](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_measured/plots) — 5 figures (QA, incident consistency, fit comparison, stability, robustness basin)
- [sources/Ip01Introduction_Belkebir.pdf](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_measured/sources/Ip01Introduction_Belkebir.pdf)

**Limitations:** unknown field units and global scale (absorbed per frequency); unpublished offset direction and loss tangent; horn pattern/receiver antenna factors unmodeled (dominant at ≥4 GHz); no uncertainty quantification beyond multi-start agreement; single-target single-file scope — no pixel-domain reconstruction and no superiority claim.