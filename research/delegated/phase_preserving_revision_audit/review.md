# Phase-preserving revision audit — detailed skeptical review

Date: 2026-09-05
Auditor: isolated deepseek-pro review worker
Scope: bounded review of the three listed files. No source file was edited. No
experiment artifact was re-run in this audit.

Files reviewed:

- `research/trispace_self_calibration/PAPER_DRAFT.md` (1054 lines)
- `research/trispace_self_calibration/GPT_PRO_HARD_PROBLEMS_ZH.md` (894 lines)
- `research/trispace_self_calibration/FINAL_SCIENTIFIC_REVIEW.md` (342 lines)

## Overall verdict

The draft is unusually disciplined: its novelty firewall, claim ledger, and
explicit separation of "proved" from "proposed" from "falsified" are in order.
Within this bounded audit I found **no fatal type or dimension error, no false
statement about phaseless SOM, and no unjustified novelty or TGRS/TAP readiness
claim**. The central method claim (PRASC-SOM beats direct joint inversion) is
correctly and repeatedly labeled open, not proven.

The defects found are refinements, not reversals: one under-specified provenance
in the common data tangent, one rhetoric that risks implying an unproven duality,
a GPT Pro prompt whose decisive question is present but ranked below a large
standard-machinery program, and an experiment plan that names matched baselines
but not a pre-registered superiority rule. Direction survives; the method claim
remains unvalidated, exactly as the draft says.

## 1. Type and dimension audit

Spot checks performed and verified:

| Object | Check | Result |
|---|---|---|
| eq (1)–(3) state differentiation | `M δj − D_Etot δχ − H_D h = 0` follows from differentiating `j = D_χ{e + Dj}` at fixed implicit dependence | correct |
| eq (5)–(7) realification | `R(A)=[[ReA,−ImA],[ImA,ReA]]` for complex coefficients; `E(A)=[ReA;ImA]` for real coefficients | dimensionally correct (2M×2r, 2M×q, 2M×p) |
| §5.1 lift (18)–(19) | `H_S = SU_r C_r + R_r`, `(SU_r)* R_r = 0`, `rank C_r ≤ p` | correct; full-row-rank `R=0` vacuity correctly labeled negative control |
| Proposition 1 congruence | `w = R_A v` turns the pencil into `[I − Q_A^T P_B̃ Q_A] w = ρ w`; singular values of `Q_A^T Q_B` are cosines of principal angles | correct; "supplemented by ones" handles the rank imbalance |
| Theorem 1 projector difference | `P_{N_r^⊥} − P_{N_{r+1}^⊥} = P_{E_r}` for nested `N_r ⊆ N_{r+1}`; equality iff `P_{E_r} B̂ = 0` | correct |
| Corollary 1 | rank-nullity on `N_r^⊥` of real dimension `m − dim N_r` | correct; uses gauge-quotiented `p_g` |
| eq (21) intensity Jacobian | `D_θ|y|²[h] = 2Re{diag(ȳ)D_θ y[h]}`; pure elementwise phase rotation annihilated | correct |
| eq (31) stacked quotient Gram | `B̃^T(I−P_Range(Ã))B̃` with `Ã,B̃` already current-projected; `Range(Ã) ⊆ C^⊥` makes the projection compose correctly | correct |
| §5.5 shared-map warning | per-frame map elimination can over-eliminate a shared map | correct |

Findings:

- **T1 (minor, under-specification).** Eq (4) at lines 239–257 mixes two
  formulations. In the reduced state-consistent form, `A_χ` and `B` are total
  full-wave derivatives. In the independent-current form, `χ` affects the data
  only through the state constraint (3), so writing `A_χ δχ` as an additional
  independent data nuisance is conservative only because the nuisance space is
  artificially enlarged. The draft says this at line 257, but does not spell out
  what `A_χ` is in the independent-current case. This is the exact double-counting
  failure mode the draft itself warns about at line 237, so the provenance of
  `A_χ` per formulation should be pinned down.
- **T2 (trivial edge case).** `g_S(ℓ,r) = (σ_{ℓ,r} − σ_{ℓ,r+1})/σ_{ℓ,1}` in eq (8)
  is undefined at `r = rank(S_ℓ)`; state `r ≤ rank−1` or a convention.
- **T3 (harmless but loose asymptotic).** Eq (23) writes
  `δ arg G_k = k δr + O(|δr|/r)`. The actual correction from the
  `1 + O((kr)^{-1})` factor is `O(|δr|/(k r²))`, so the bound is valid but weaker
  than necessary. Stating the sharper bound preempts a reviewer who notices the
  `O(|δr|/r)` term can look comparable to `k δr` at small `r`.

## 2. Theorem statements and assumptions

- **Theorem 1** (line 576): the preamble "fixed pose, map, frequency, whitening,
  and gauge quotient" is sufficient — it pins `Â, B̂, S`, which is exactly what the
  nested-nuisance argument needs. Residual risk: in a *reduced* formulation whose
  total derivatives `A_χ, B` depend on the retained parameterization `U_r`, the
  fixed-`Â` premise silently fails. Recommendation: add one sentence restricting
  Theorem 1 to the independent-current/data-block formulation with rank-independent
  `Â, B̂`, and explicitly decline to assert monotonicity for rank-dependent reduced
  formulations.
- **Proposition 1**: assumptions are adequate (PD `K_0` on the selected support, no
  pose prior). The restriction-to-support is only in the proof ("`R_A` nonsingular
  on the selected support"); promote it to the statement.
- **Proposition 2**: "differentiation and integrability conditions that allow score
  conditioning" is an acceptable draft-level gesture, but a submission must name
  the conditions (score integrability, exchange of ∇_θ and expectation). The
  equality claim (`s_Y` Z-measurable a.s.) is stated correctly.
- **Corollary 1**: correct. Uses the gauge-quotiented `p_g`, which is the right
  version; see cross-document note in §8.

No theorem statement was found whose conclusion overreaches its stated hypotheses.

## 3. Statements about phaseless SOM

No false statement found. The draft explicitly rejects the strawman
("phaseless SOM has no current-space constraint") in §2.1, in the abstract, in the
claim ledger (status "rejected"), and again in the GPT Pro prompt's red lines.
The characterization — Pan–Zhong–Chen–Yeo phaseless SOM still uses spectrum
analysis and a deterministic/ambiguous contrast-source partition, with the
deterministic recovery modified for missing phase — is the correct, non-strawman
account. This is a positive finding; prior drafts of this project had this point
as a live risk and it has now been correctly closed.

Uncertainty: attribution was checked against project knowledge and memory notes,
not by re-reading the 2011 TGRS primary source in this audit.

## 4. Conflation of ρ/K_eff with B_vis/J_x

Explicit separation is present and correct:

- line 442–448: "They are not the same matrix and need not have the same
  dimension," with three distinct gates in box (17);
- the non-monotonicity of `ρ(r)` in rank is stated and enforced in the algorithm
  (feasible-set search rather than "everything below a threshold is feasible");
- the GPT prompt P2.3 forbids equating the spectra and demands a precise duality
  theorem, and the red lines repeat the prohibition.

One residual risk: "two sides of the same local confounding geometry" (line 442)
is rhetorical and could be misread as a proven duality. The only precise relation
established in the draft is that `J_x(r)` is the pose block Schur complement of the
joint Gram after current+map elimination, while `K_eff(r)` is the map block Schur
complement after current+pose elimination — different blocks, different nuisance
orderings. Recommendation: either state that relation in one sentence or mark the
duality as an open theorem (it already sits in P2.3). No actual error found.

## 5. Novelty and TGRS/TAP readiness

No unjustified claim found. The contribution table honestly grades the Loewner
monotonicity as "derivable linear algebra, physically specialized here"; joint
pose+image, Schur/Fisher/principal-angle machinery, and frequency continuation are
all quarantined as prior art or textbook; the "narrow candidate contribution" is
three SOM-specific hypotheses. TGRS and TAP are discussed only as routes with
explicit remaining gates (Section 9), not as readiness claims, and the limitations
section restates bounded prior-art retrieval. Claim-ledger statuses match the
evidence descriptions. No action needed.

## 6. Does the GPT Pro prompt prioritize the decisive problems?

Directionally yes, with three weaknesses.

The genuinely decisive question — whether PRASC-SOM has any non-trivial advantage
over a well-tuned direct joint inversion ("why must it be SOM") — is present as
P4.3, is explicitly labeled non-avoidable ("不能回避"), and the quality red lines
order a downgrade to theory/diagnostic result if no advantage can be proven or
designed. That is the right fail-safe, and it is in place.

Weaknesses:

- **Ranking.** P4.3 sits at the third sub-position of the "second highest
  priority" block, behind the very large P0–P3 program. Much of P1.1 (Fisher
  information data processing) re-derives Proposition 2, which is already proved
  in the draft. The genuinely open parts of P1 are P1.2 (nuisance-eliminated
  efficient information, correctly flagged as not following automatically from
  blockwise Loewner order) and P1.3 (the phase-basin gate). Risk: budget is spent
  on textbook re-derivation before the gating P4.3 verdict.
- **Anchoring.** Lines 6 and 29 frame the deliverable as "theory and design
  sufficient to support a TGRS/TAP method paper." Even with the red-line
  guardrails, this invites a confirmatory package. Recommendation: re-anchor to
  "produce the strongest theory package **and** an explicit pass/fail verdict on
  the SOM-specific advantage, including the downgraded positioning if it fails."
- **Volume.** P0–P7 plus S1–S3 plus a nine-part deliverable format in one prompt
  risks shallow per-item treatment. A hard stop after the P4.3 verdict, before
  expanding P6/P7 prose, would protect depth where it matters.

Net: the prompt asks for the right things; the ordering and framing could let the
right thing be crowded out.

## 7. Failure tests against direct joint inversion

What is already present:

- E4 lists direct joint full-wave inversion with the same optimizer and priors,
  plus phaseless SOM, wrong-geometry SOM, fixed-rank SOM, and known-pose/known-map
  oracles; metrics include failure rate, runtime, forward/adjoint counts, and
  "matched cost or accuracy";
- the inherited loop already contains the decisive negative control: the old
  restricted solver is either a reparameterization of the direct objective or
  freezes hidden pose and performs *worse* (Section 7, item 4; review §5.1 table);
- the draft states the method claim fails absent a statistically meaningful
  advantage.

Gaps:

- **No pre-registered superiority rule.** "Statistically meaningful advantage at
  matched cost or accuracy" (line 917) is not operationalized: which metric is
  primary, what effect size is decisive, hypothesis test, seed budget,
  multiplicity correction. Without this, post-hoc interpretation remains possible.
- **No cycle-skipping failure test of the frequency gate.** The known half-circle
  pseudo-minimum (pose error ≈ 0.86, review line 221) is the natural adversarial
  control for gate (36); E4 should include a run where the gate either prevents or
  fails to prevent a wrong branch.
- **No expected-loss test.** A bound on the claim requires at least one case where
  PRASC-SOM is predicted to lose — e.g., phase corruption that is not
  pose-dominated (clock drift, mutual coupling, unmodeled phase noise). E4 lists
  baselines but no such adversarial condition.
- **Matched-cost is metric-only.** The protocol lists forward/adjoint counts but
  does not fix an equal forward-solve budget or a tie-break rule between accuracy
  and cost.

## 8. Cross-document consistency

- Corollary 1 uses the gauge-quotiented `p_g = p − dim G_x`; FINAL_SCIENTIFIC_REVIEW
  §4.2 uses the unquotiented `p = 3` with the global SE(2) gauge discussed
  separately (lines 139–140). Both are internally consistent, but the two
  conventions should be harmonized so a reader does not mix `p` and `p_g`.
- The `r ≥ M−1` nuisance-saturation result is consistent between paper §7.4 and
  the review's corrected table (`M=8,12,16`, codim `3,1,0`).
- The review's state-hiding iff (`R_U h = 0` and `T_U h = 0`) is presented as a
  correct solvability statement that never became binding; the paper demotes the
  state witness to a negative constraint. Consistent.

## 9. Verified-correct ledger (for the record)

Correct as stated: independent-current/state differentiation; realification
conventions; retained lift and full-rank vacuity; Fisher data-processing
inequality and its equality condition; intensity-Jacobian phase cancellation;
large-`kr` phase/amplitude scaling and `O(k^{-1})` ambiguity spacing;
principal-angle form of `ρ_i`; Loewner monotonicity with explicit equality and
strict-loss conditions; dimension obstruction; stacked shared-map quotient Gram;
Schur form with pose prior; the non-monotonicity of `ρ(r)` in rank; the claim
that per-frame map elimination can over-eliminate a shared map.

## Uncertainties

- Bounded to the three listed files. The numeric evidence (residual `5.5e-16`,
  leakage `0.784`, projector jump `1.424`, median error tables, `r ≥ M−1`
  saturation) was read as reported and not re-derived from raw loop artifacts,
  `validation/rank_saturation_results.json`, or figures.
- Phaseless-SOM attribution and the reference list were checked against knowledge
  and project memory notes, not against freshly opened primary sources; DOIs were
  not live-verified.
- No nonlinear or numerical re-run was performed, so the claim that the old
  negative controls transfer to PRASC-SOM is taken on the draft's own word.

## Bottom line

The phase-preserving direction survives this audit as a correctly-hedged theory
foundation. The five items worth parent attention (eq (4) provenance; Theorem 1
scope sentence; §4.4 duality rhetoric; prompt re-anchoring/ranking; E4 decision
rule and adversarial controls) are refinements, not blockers. The central
algorithmic claim remains open and must not be described as validated.
