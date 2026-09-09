# Independent severe reviewer audit — A2 paper draft

Role: isolated independent reviewer (no child agents). Evidence sources read:
`PAPER_DRAFT_A2.md`, `THEORY_AUDIT.md`, `A2_PRASC_SOM_THEOREM_PACKAGE_EN.md`,
`a2_solver_audited/{AUDIT.md,frozen_parent.json,cert.py,calibrate.py,reduced.py,
objective.py,final_200.jsonl,final_800.jsonl}`, `a2_highdim/summary.md`,
`a2_theory_validation/replication/summary.md`, `a2_literature/EVIDENCE.md`,
`a2_physics/summary.md`, `a2_measured/SUMMARY.md`,
`A2_PRASC_SOM_VALIDATION_PROTOCOL.md`.

Independent verifications executed this pass (not inherited):

1. `partial_r log G_2 = ik - 1/(2r) + O((kr^2)^{-1})`: confirmed numerically
   against scipy Hankel functions (real part converges to -1/(2r)). Paper sign
   and error order correct.
2. Passive-bound identity `Im D = Q - delta I` with
   `delta = k^2h^2/4 - pi k a J_1(ka)/2`: residual ~1e-15 for
   N in {8,16}, k in {pi,4pi}; `Q` Gram matrix numerically PSD (min eigenvalue
   >= -2.8e-16). Full Proposition 1 chain re-derived: correct.
3. Nonmonotone `rho` example A=(1,1,0), B=(1,0,1), C0/C1/C2: recomputed
   directly, retention sequence exactly 3/4, 1, 0. Correct.
4. Gradient-error bound `eps_J ||r|| + (||J||_F+eps_J) eps_r` and conditional
   descent condition `kappa_g < h_-/h_+`: re-derived, correct.
5. `cert.py`/`reduced.py`/`certificates.py` implementation matches paper
   Section 5; the 12 positive-loss dense checks enclose independent errors
   (asserts pass).

No severe mathematical error found in the checked theory core. The serious
problems are experimental validity, novelty, statistics, and application
evidence — not the algebra.

## Verdict

The paper is exceptionally honest about scope, but its empirical core does not
yet support publication: the frozen E4 comparison exercised the reduced method
zero times, the two "passing" statistical gates are degenerate artifacts, the
only positive nonlinear evidence uses the exact-direct method, and the
measured-data pilot does not validate geometry recovery. The theory section is
defensible as a limits-and-design contribution; the experimental claims are not
publication-ready at TAP or TGRS.

## Findings

### F1 — SEVERE. The frozen E4 comparison is void with respect to the reduced method.

Evidence: `final_200.jsonl`/`final_800.jsonl` ledger fields; `E4_RESULTS.md`;
paper Sec. 7.3. Every method has 0/240 joint successes at both budgets.
`prasc` records 0 certificate probes and 1 fallback event at BOTH budgets —
the guarded method is a pure direct fallback in 480/480 runs. `fixed_rank`
performs 0 objective evaluations at budget 200 and exactly 1 rejected
certificate probe plus 0 objective evaluations at budget 800
(`finished_budget_policy=False` everywhere). `reduced_moves=0` in all 2,400
records. The paper discloses most of this in prose ("the experiment did not
adequately exercise that pathway"), but the abstract's "frozen matched-budget
comparison" phrasing and Sec. 7.3's framing still present the comparison as an
evidence category.

Fatal to: any claim that guarded SOM-informed reduction was tested or compared.
Nonfatal to: the algebraic identities and the paper's own withdrawal of
superiority. Revision: state in the abstract and in 7.3, first sentence, "no
reduced iteration was accepted in any of the 2,400 runs; the frozen comparison
tests only the cost-accounting protocol, not the method."

### F2 — SEVERE. The causal failure mechanism is cost accounting plus stage scheduling, not the method.

Evidence: one record shows `ledger.units=158.667` with
`stage_evaluations=[0,0,0,1]` — a single full-band f/g evaluation at the final
stage consumes the 200 cap. Fixed stage fractions force zero low-frequency
evaluations; the skinny-to-dense matvec conversion inflates basis/probe costs.
At budget 800 direct/prasc/phaseless all stop at exactly 714.00 units
(~4.5 f/g) and direct pose error worsens 0.27083 -> 0.32192 m. Paper Sec. 7.3
identifies this, but it belongs in the abstract as protocol invalidation.

Minimum discriminating test (cannot be fixed retrospectively): re-derive cost
units dimension-aware (rank-dependent), pre-register stage budgets that
guarantee >=1 accepted low-frequency step and a completed final-data stage,
then require the original compound gate on fresh scenes (no seeds 1-10 or
1001-1020). The comparison is only meaningful if at least one reduced step is
accepted.

### F3 — SEVERE (method viability). The certificate is quantitatively too loose where it matters.

Evidence: `results/passivity_checks.json` bound/true ratios: state 13x-210x,
derivative 41x-9,051x, increasing with wavenumber. `results/enrichment_probe.json`:
at the top frequency even rank-128 enriched bases (half the 256-cell space)
fail admission (data-error bounds 1.16-18.2 vs threshold 0.6; relative
gradient bounds 0.94-49 vs 0.1). Non-vacuity holds only at low frequency at
rank 72 (28% of full rank), where no computational advantage is demonstrated.
Paper Sec. 7.2 says "non-vacuous in some regimes" without quantifying
tightness. Fatal to the implicit hope that certificates enable affordable
reduction; nonfatal to the bound's correctness. Revision: publish the
bound/true ratio table and state the high-frequency admission failure as an
explicit method limitation.

### F4 — MAJOR. Novelty is thin against verified prior art.

Evidence: `a2_literature/EVIDENCE.md` (fulltext) — Idriss/Raj 2025 already
jointly optimizes complex per-transmitter calibration with multifrequency SOM
on Fresnel measured data; Bellomo 2014 handles phase centers in diffraction
tomography; Huang/Nammour/Symes 2017 source-receiver extension is the standard
joint data/geometry FWI framework. The paper's own tools (data-processing
inequality, canonical correlations, Schur complements, bias identity) are
textbook, and the paper concedes this. The defensible core is the typed
rank/bias/acquisition budget plus the model-conditional passive resolvent
bound; the passive bound is an elementary Gram-kernel + diagonal-correction
argument whose novelty the paper itself leaves unestablished. Nonfatal to the
typed-budget framing; fatal to any broader "SOM self-calibration" novelty.
Revision: scope all claims to the typed budget framework and add an explicit
positioning of the passive bound against VIE discretization-norm literature.

### F5 — MAJOR. Two of the three primary gates are vacuous, not passed.

Evidence: `E4_RESULTS.md` prints gates `[False, True, True]`. With identical
outputs across methods, the normalized pose/map upper bounds are exactly 0 and
the success-difference lower bound is exactly 0 — degenerate artifacts, not
noninferiority evidence. Paper Sec. 7.3 explains this in prose but the result
table and figure caption do not label the two True gates vacuous. Revision:
render them as "vacuous (identical outputs)" and never count them toward any
acceptance summary.

### F6 — MAJOR. The efficient Fisher contraction is physically untested.

Evidence: replication summary — E1 physical tangents (N=8, n=10, sigma^2=1)
verify only the FULL-joint Loewner order J_coh >= J_ph; the efficient
(post-nuisance-elimination) inequality in Theorem 1 is checked only on
synthetic models. Physical E3 bias-risk MC is evaluable in 0/20 records
(correctly declared empty). The abstract's "even after nuisance elimination"
is theorem-level on physical data. Revision: state explicitly which level is
physically verified. A physical efficient-Fisher check is cheap and should be
run.

### F7 — MAJOR. Application/measured evidence does not validate geometry recovery.

Evidence: the only positive nonlinear result (larger-map study,
`a2_highdim/summary.md`) uses the exact-direct solver, not the guarded method;
6 scenes x 2 radii, descriptive cluster bootstrap, restricted 49-Gaussian
family. The measured pilot (`a2_measured/SUMMARY.md`, completed during this
review) has no position-error ground truth; eps_r drifts 5.2 -> 2.5 -> 1.2 as
the fitted band grows (outside the published 3 +/- 0.3) and 4-8 GHz residuals
reach 0.4-0.9 under plane-wave misspecification. The paper's
`MEASURED_PILOT_PENDING` must be replaced with this non-validation result; it
cannot be cited as measured-data readiness.

### F8 — MINOR. Identity tolerances are loose and quoted without scale.

Evidence: Sec. 7.1 quotes "maximum backward-scaled residual 0.274" against a
test threshold of <100 (`test_e2.py`), a roundoff-level gate. Not a math error
(0.274 is relative to the backward-error scale), but quoting the number
without the threshold invites misreading as a 27% identity failure. Revision:
state the tolerance whenever the residual is quoted.

### F9 — MINOR. "Certificate" terminology overstates unless qualified.

Evidence: title and keywords ("Physical-State Certificates",
"numerical-physics certificate"). The bound is correct, but (a) it is
evaluated in floating point with an ad hoc 1e-12 margin, not interval
arithmetic — the paper concedes this in Sec. 5; (b) availability is
model-conditional (positive common loss, u>0 on all cells, homogeneous
background, fixed quadrature); (c) it bounds numerical error at a candidate,
not branch or scene confidence — also conceded. The prose is disciplined; the
title/keywords are not. Revision: "Model-Conditional State/Residual Error
Bounds" in the title, or qualify every occurrence of "certificate".

### F10 — MINOR. Policy completion must not read as success.

Evidence: at budget 200 direct/prasc/phaseless have
`finished_budget_policy=True` with 0 successes and one final-stage evaluation
and no accepted iterate; at 800 fixed_rank has `False` with 585.39 charged
units. Paper states "status alone can never establish success", but the
mapping between `finished_budget_policy` and the recovery criteria is not
printed in the results. Revision: report completed-policy counts alongside
success counts in the E4 table.

## Resolved vs unmet gates

Resolved manuscript issues: scaling-bug chain rule and square-sum stacking
(audited copy), pre-test initialization amendment, exact Rice intensity law on
the same parent, matched direct fallback, honest separation of r_free/r_num/
L_det, honest labeling of the empty E3 MC and the all-zero E4 result.

Genuinely unmet publication gates: (i) a discriminating matched-budget
comparison that exercises the reduced method (F1/F2); (ii) a certificate with
a demonstrated useful operating point (F3); (iii) any measured/3D validation
with known position error (F7); (iv) novelty positioning against
Idriss-Raj/Bellomo/extension literature (F4); (v) the protocol's own TAP/TGRS
evidence gates (validation protocol line 177) — unmet.
