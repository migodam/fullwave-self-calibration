# Bounded adversarial review of the current A3 manuscript

Date: 2026-09-08. Reviewer: isolated Codex worker; **not DeepSeek Pro and not an independent external scientific review**. Parent Codex retains mathematical validity, novelty adjudication and final scientific interpretation.

Scope: read AGENTS.md and the five requested manuscript/protocol/ledger/audit files. No external literature lookup, source-code audit, new experiments, workers, pipelines or credential inspection. Findings below concern internal evidence/claim alignment; proposed tests were **not executed**. References to paper lines refer to the version read during this review.

## Overall assessment

The manuscript is unusually explicit about development scope, failures, additional-data interventions and absence of SOM superiority. Those disclosures should be preserved. The numerical multifidelity conclusion agrees with the supplied audit; I found no numerical contradiction in that comparison. Five consequential objections remain, including one concrete terminology contradiction and one important bridge between experiments that is missing from the manuscript. Other objections are honestly acknowledged validation gaps, not concealed failures.

## 1. Phase metric semantics remain inconsistent outside the correction paragraph

**Type: genuine presentation/interpretation defect; repair before sharing the draft.**

The abstract says increasing frequency “worsens held-out phase predictions” (paper line 7). Section VI.D correctly states that the original phase statistic uses the independent reference solver, estimated material/pose and **no fitted electronics**, and explicitly calls this an oracle-model parameter-transfer diagnostic (line 136). Nevertheless, lines 149–151 again call the original 0.267/0.256, 0.0164, 0.0540 and 0.0375 rad statistics “held-out phase” without the oracle qualifier. Actual inverse-model sensor-phase errors appear separately at line 155 (low-band 0.00481/0.00371, high-band approximately 0.0173/0.0170 rad). A reader can reasonably treat these as conflicting reports of one deployable prediction metric. The abstract's qualitative claim may be supported by the actual prediction audit, but its attribution is currently ambiguous.

**Minimal repair and check:** name every old statistic “oracle-model transfer phase RMSE”; reserve “sensor prediction” for fitted-electronics inverse-model evaluation. Tie the abstract frequency claim directly to the actual matched low/high sensor-prediction comparison. Make a compact metric ledger listing evaluator, fitted electronics inclusion, held-out channels, aggregation and source artifact. Cross-check each numerical occurrence against that ledger. This is principally a frozen-artifact audit, not a reason to refit.

## 2. The nonlinear acquisition results do not test the nuisance model motivating the graph rule

**Type: genuine missing experimental linkage; already disclosed in the critic ledger, insufficiently explicit in the paper.**

Section III-B assumes unrestricted separable complex gains (lines 56–60), and Section V combines graph screening, nuisance-profiled information and branch separation into an intended procedure (line 100). The four-scene nonlinear result then presents fixed/random/local-information/branch outcomes (line 177). The critic ledger explicitly says this experiment has **known electronics** and does not validate the unrestricted-gain sparse-graph policy (CRITIC_AND_REPAIR.md line 17). The manuscript does not repeat that decisive nuisance distinction alongside the nonlinear result. Its warnings against policy superiority do not fully explain which mechanism was actually tested.

**Minimal repair and test:** state the known-electronics condition next to the four-scene result now. To establish the missing bridge, use a small predeclared nonlinear case with the stated separable gain freedom and material profiling; compare equal-edge acquisitions that do versus do not close the necessary cycles, charging all initial fits and refits. Record null directions, bank coverage, accepted parameter error and rejection separately. This would be a mechanism check only; a population advantage still needs final scene-level comparisons.

## 3. “Reliability-controlled calibration” currently describes a direction more than a validated continuum-fidelity decision rule

**Type: honestly disclosed substantive gap; do not treat discrete numerical checks as its repair.**

The discrete controller uses exact discrete-model acceptance and residual/output tolerances (line 88). The paper correctly separates this from continuum discrepancy (lines 86, 124–130), and the branch guarantee requires coverage and a justified discrepancy bound (lines 92–102). Yet line 149 says high frequency “is admitted only after adequate forward-model validation,” without a quantitative admission rule or a demonstrated operational accept/refine/reject experiment. The critic ledger agrees that a tested fidelity/admission controller is still needed (line 15). The abstract wisely says the observations *motivate* reliability control; retain that conditional posture.

**Minimal repair and test:** either label the frequency-admission sentence explicitly as a proposed requirement, or freeze one implementable discrepancy/admission criterion before a small validation exercise. Use an independent reference only for final assessment, not an inaccessible online oracle. Include one deliberately underresolved and one adequately resolved case, and test whether the criterion selects refine/reject/admit consistently with its declared target. If no defensible error bound is available, call the rule heuristic and report false admissions. A converged GMRES solve or a small fine/coarse difference alone is insufficient.

## 4. No demonstrated incremental decision benefit yet joins the established ingredients into the paper's core contribution

**Type: openly acknowledged contribution/validation gap, not a new allegation of false novelty.**

The paper attributes closure/projection/reduction to prior work (lines 15–21, 60, 86), does not assert first-discovery priority for the Maxwell specialization (line 68), finds no ROM time advantage (lines 114–118), and rejects multifidelity correction beyond coarse warm starts (lines 165–167). The acquisition population result remains open. Consequently, the evidence supports a careful diagnostic/mechanism paper, but does not yet identify a validated new algorithmic benefit. The conditional far-field mechanism may be scientifically useful; its originality and venue sufficiency are for the parent and external scientific reviewer, not this worker.

**Minimal repair and test:** select one explicit additional claim and attach a decisive test to it. For example, test whether the predeclared nuisance-aware acquisition/fidelity decision improves a stated recovery or false-acceptance outcome against a strong equal-resource baseline on unopened scenes. If the intended article instead centers on a physical mechanism, explicitly frame its scope and complete the closest-prior comparison; do not demand an algorithmic novelty claim it does not make. No amount of repeated algebra tests substitutes for either route. ACCEPTANCE.md correctly leaves these gates open.

## 5. The phaseless comparison controls the observations but not optimization reliability

**Type: honestly disclosed limitation with a locally overbroad inference.**

The Rice likelihood uses magnitudes of the same noisy complex data and removes exactly insensitive phase/delay parameters (line 153), which is an appropriate observation-model control. But the same paragraph infers support for coherent data from two direct fits using one common nominal initialization. A shared starting point does not equalize success probability in two different nonconvex landscapes. The reported gap could combine information loss and optimization failure. The paragraph discloses this initialization limitation and explicitly avoids claiming dominance over classical phaseless SOM, so this is not a hidden baseline omission.

**Minimal repair and test:** phrase the existing result as a comparison of these two fitted procedures from the declared initialization. For an information-based interpretation, use a predeclared matched multistart/wall-budget experiment on the same data with selection by each method's training likelihood, and compare held-out predictions plus pose/material error. Preserve failed starts and report all optimization cost. A local profiled-information comparison can supplement, but cannot by itself settle nonlinear recovery.

## Multifidelity audit: checked and internally consistent

- Paper mean times 20.93, 3.51, 17.13, 34.97 and 20.52 seconds agree with the JSON summaries for fine-direct, coarse-only, coarse-warm, value-corrected and tangent-corrected methods.
- Tangent/warm paired ratios 1.2117107373 and 1.1845730806 support the stated 18.5–21.2% slowdown. Mean fine RHS 120 versus 144 is correctly reported and does not imply a wall-time benefit.
- Coarse-warm and tangent endpoints have relative fine-objective gaps of order 1e-14, supporting discrete endpoint equivalence. Coarse-only gaps exceed 3.1 and its mean material error is 6.926%, versus 2.287% for fine-direct.
- The top-level JSON `passed: true` is not evidence of superiority; its explicit decision rejects the correction advantage. The paper and acceptance ledger interpret it correctly.
- The protocol promises sensor-prediction and structural-field metrics, while this audit JSON exposes timing, work, pose/material and objective summaries only. Those metric records and optimizer termination flags cannot be verified from the supplied audit alone. In particular, the paper's value-only outer-limit statement is consistent with the ledgers but not independently substantiated by fields in this JSON. This is an audit-coverage limitation, not evidence that the statement is false. Link the raw per-fit records before final reproducibility review.

## Explicit non-objections

Do not reopen already corrected issues as if the draft still asserted them: numerical rank is separated from free-current nuisance; extra electronics references are additional data; ADDA/VIE share a method family; noise replicates are not independent scenes; bound-saturated Fresnel radii are not measured antenna recovery; value-only and tangent corrections are negative ablations. Final statistical and hardware gates remain open rather than falsely marked complete.
