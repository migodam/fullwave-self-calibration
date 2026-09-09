# Bounded audit of the frozen ROM collection and analysis

Date: 2026-09-08. Reviewed only ROM_FINAL_PROTOCOL.md, timed_rom.py, rom_final.py and analyze_rom_final.py. No results/interim outcomes were opened; no code was executed, edited or interrupted. This is source-level review, not execution verification. Parent owns statistical and scientific adjudication.

## Disposition

**No unconditional collection-invalidating error found in the inspected files.** The paired success-bound, median interval and sign-test implementations match the frozen restricted estimands. Deadline checks reject ordinary completed trials after the deadline and charge completed forward work before raising the deadline exception. Two conditional failure-accounting issues and several verification limits should be checked after the fixed collection ends, without changing frozen methods or tuning against outcomes.

## Statistical implementation

1. **Paired success difference: consistent and conservative.** Twofold-only wins and Krylov-only wins are counted over the same 40 scenes (analysis lines 34–44). The binomial marginal CP formulas handle zero/all counts correctly. Taking a 97.5% one-sided upper limit for the win probability minus a 97.5% one-sided lower limit for the loss probability gives at least 95% one-sided coverage by the union bound; independence of those two counts is not required. The complementary lower bound uses the same valid construction. These separately valid one-sided 95% limits should not be relabeled together as an exact central 95% interval.
2. **Speed eligibility and estimand: consistent.** Both methods must succeed and satisfy the relative objective-gap threshold (analysis lines 38–42). The denominator is Krylov's objective, a concrete one-sided normalization of the protocol's relative-gap language. The reported statistic is a conditional median of consumed-time ratios, not mean speedup or time to first successful recovery. It includes nonpreemptive overrun by design. Conditional selection is explicitly disclosed, and zero eligible pairs yield no speed rejection.
3. **Exact median interval: correct order statistics.** The largest k satisfying twice the Binomial(n, 0.5) lower tail <=0.05 gives [X_(k), X_(n-k+1)] (analysis lines 18–23). For positive runtime ratios, [0, infinity] is a valid uninformative fallback at insufficient n. The JSON uses `None` for that infinite upper endpoint; reports should label it as unbounded, not a missing finite estimate. Discrete/tied distributions can make coverage conservative.
4. **Sign/binomial test: correct direction.** The count of ratios <=0.8 has lower-tail p-value under p=0.5 (line 45). Under a median-at-most-0.8 null, P(ratio<=0.8)>=0.5, so the lower-tail rejection is valid and conservative at ties. Dropping ties would change this particular implementation and is not necessary. The protocol's requirement to reject both component nulls before excluding the union of the two prespecified advantage routes is an intersection-union test; no Bonferroni division between these two component tests is required for that exclusion claim.
5. **Inference scope is honestly narrow.** Seeded scene independence, runtime stationarity and conditional eligibility are assumptions, not proven by these scripts. Cyclic method order plus serial execution mitigates order effects but cannot by itself establish independent identically distributed machine timing. The script explicitly states the stationarity assumption. No population-wide, q49 or classical SOM claim follows.

## Deadline and work accounting

- `forward` checks before and after the physical solve, charges returned work before the post-solve check, and `record` checks immediately before selecting an accepted state (timed_rom.py lines 26–41). Adjoint calls likewise charge returned work before checking time (lines 47–56). Chart construction/evaluation and GN steps have surrounding checks. I cannot verify the imported helpers' internal work counters because they are outside this task's allowed files.
- Setup time is included through the start offset (timed_rom.py line 20; rom_final.py lines 50–53), and the timed solver preserves only the last accepted state. Truth-dependent errors and final physical-loss evaluation occur after `online` is captured (lines 98–112), so they do not directly affect admission.
- The implementation is an ordinary nonpreemptive software deadline, not hard real time. `record` timestamps immediately before copying/publishing the state, and ledger construction precedes initialization of `start`; small unclocked bookkeeping gaps also occur between setup measurement and entry. These are precision limitations rather than evidence of accepting expensive numerical outputs after deadline. Do not describe eight seconds as a hard wall-clock execution cap.
- The reported `atomic_overrun_seconds` is total positive elapsed overrun; Python bookkeeping and time-check latency can be included, not only one isolated linear-algebra primitive. The protocol already disclaims exactly equal consumed time.
- Initial nominal parameters are the fallback endpoint when there are no recorded events. They are chosen before the online solves and need no truth information; the offline evaluator must not be mistaken for an accepted late optimization update. If model setup alone consumes the budget, report that no evaluated iterate was available online.

## Conditional defects requiring post-collection accounting checks

### A. Offline evaluation exceptions are conflated with online optimization failures

`timed.solve` evaluates the frozen accepted state after timing has ended (lines 102–113). If this offline forward evaluation or error computation raises an exception, the enclosing catch in rom_final.py lines 55–58 writes `joint_success=False` and only the generic exception fields. It loses the accepted state/events and labels the run identically to a failure during optimization. This would bias the success estimand if an otherwise valid decision endpoint were classified as a failure solely because the separate evaluator failed. The exception-row elapsed time would also include offline activity.

**Disposition:** conditional coding/accounting defect; no evidence that it occurred because outcomes were not inspected. After collection, identify exception stages from retained evidence. If any audit-only failures occurred, do not silently treat them as algorithm failures or rerun optimization selectively. Parent should decide and document whether frozen-endpoint recovery is possible without changing decisions, or whether the affected comparison is incomplete. If there are no relevant exceptions, this path does not invalidate the completed normal rows.

### B. Generic exception rows cannot meet the complete-work/endpoints reporting promise

The broad catch preserves exception type/message and elapsed time but omits ledger counters, last accepted parameters/events, setup time and overrun. An exception inside an atomic model/helper call can also occur before its work dictionary is returned. Therefore “all work” and “all endpoints” cannot be verified for such rows from this record structure alone. The analysis treats absent overrun as zero in its maximum calculation (line 52), which understates or obscures unknown overrun for exception rows.

**Disposition:** reporting defect conditional on exceptions. Report unavailable work/overrun as unavailable and separate exceptions from normal deadline exits. Retaining exceptions as prespecified unsuccessful runs is otherwise a defensible solver-benchmark policy when the failure truly occurs online; the concern is lost accounting and stage ambiguity, not a demand to delete failures.

## Verification boundaries and provenance

- The runner enforces serial method calls and cyclic order; one numerical thread and absence of concurrent heavy jobs are not enforced in these four files. They may be set by the launch environment. Verify from launch records, not by assuming the Python script supplies them.
- Model frequencies, seed/rank details and imported chart guards/work counters are inherited from modules outside this review. `case` calls the truth forward without explicit frequency IDs while background/noise scaling supplies explicit IDs; actual compatibility depends on the imported model defaults. No mismatch is established here, and checking those defaults would exceed the requested scope.
- Resume checks bind the protocol and listed solver sources and reject mixed hashes. The analyzer checks complete seed/method coverage and a single pair of hashes. Its own source is not included in the recorded solver manifest and it does not compare stored hashes to current files. Preserve the frozen analyzer separately; these checks establish internal consistency of the collection records, not complete independent provenance authentication.
- This review does not determine whether 40 scenes will be informative enough. The frozen protocol correctly requires “unresolved” when bounds cannot exclude the prespecified benefit and does not permit enlargement or outcome-driven substitution of secondary comparisons.

## Parent handoff

Continue the frozen collection without changing the reviewed algorithms on the basis of this audit. Once the full run ends, check exception stages, launch-thread conditions and provenance; report any resulting limitations explicitly. No interim performance finding or numerical acceptance judgment has been made by this reviewer.
