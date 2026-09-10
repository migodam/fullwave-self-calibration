# FINDINGS — increment 2 (Part A q' arrangement sweep, Part B finite-risk checks)

Interpreter `/Volumes/.../a3_research/.venv3d/bin/python`, mpmath 1.3.0 via `PYTHONPATH=/Users/migodam/.cache/uv/archive-v0/VjSSQm31390egyfy`, numpy 2.5.3. Part A 20.5 s, Part B 5.6 s.

Artifacts
- Part A: `results/qprime_variants_20260910T094847Z.json` + `logs/qprime_variants.log` (authoritative)
- Part B: `results/finite_risk_20260910T095306Z.json` + `logs/finite_risk.log` (authoritative)
- Superseded, step-2 cross-check label only, same numbers: `results/finite_risk_20260910T095123Z.json`, `...T095240Z.json`
- Quarantined: `results/interval_20260910T093248Z.json` -> `results/PARTIAL_interval_20260910T093248Z.json.bad` (truncated JSON at char 3690, content unchanged, sha256 `3158515f...c516`); noted in the Part A log.

## Part A — min-over-boxes lower bound of dq/deps, iv.dps=35, same boxes (2500 / 3000)

| variant | I1 q'_lower | rel vs supplied | I2 q'_lower | rel vs supplied |
|---|---|---|---|---|
| supplied | 4.3262426032187571e-04 | — | 1.3086650394735435e-04 | — |
| V1 A5 chain (used by `mie_reactance`) | 4.293019592886e-04 | -0.768 % | 1.294368394441e-04 | -1.092 % |
| V2 closed-form dD_j/dz | 4.282819918549e-04 | -1.004 % | 1.290762161586e-04 | -1.368 % |
| V3 log-derivative | 4.266188517697e-04 | -1.388 % | 1.287354660101e-04 | -1.628 % |
| **V4 fully expanded** (tightest) | **4.431050126934e-04** | **+2.423 %** | **1.357809867428e-04** | **+3.755 %** |
| V5 dps=25 / dps=50 | 4.293019592886e-04 | -0.768 % | 1.294368394441e-04 | -1.092 % |
| V7 centered (eps = mid + delta) | 4.293019592886e-04 | -0.768 % | 1.294368394441e-04 | -1.092 % |
| V6 midpoint - 0.5h\|q''\| / -1.0h\|q''\| | 4.551452649037e-04 / 4.550614172678e-04 | +5.2 % | 1.399780554883e-04 / 1.399551402336e-04 | +7.0 % |
| V6 midpoint only (non-rigorous) | 4.552291125396e-04 | +5.2 % | 1.400009707429e-04 | +7.0 % |

Envelope: I1 `[4.266188517697e-04, 4.552291125396e-04]` (range 2.861e-05); I2 `[1.287354660101e-04, 1.400009707429e-04]` (range 1.127e-05).

- **No variant reproduces the supplied value within 1e-9 relative** (V1..V7 all fail). The supplied value lies strictly *inside* [V4, V6], so it is a valid but differently-grouped enclosure; V1/V2/V3 are looser, V4 is the tightest rigorous bound found.
- V5 (dps 25/35/50 identical) and V7 (identical to V1) show the answer is arrangement-driven, not precision-driven.

## Part A — location of the infimum (rigorous)

- Dense `mp.dps=50` scan, 20001 points: **argmin = eps 4.0** (I1, min `4.551544685356e-04`) and **argmin = eps 5.0** (I2, min `1.399811270778e-04`) — both at the **right endpoint**.
- Interval q'' via the section-8 chain rule applied twice, dps=35, all 5500 boxes: **strictly negative** (max q''_upper `-1.315836e-04` I1, `-3.380447e-05` I2) => q' strictly decreasing => inf q' at the right endpoint. ✔
- Point check: max rel diff analytic q'' vs Richardson FD of q' = `4.478e-19` (step-dependent 4.5e-11 / 4.5e-15 / 4.5e-19 for h=1e-2/1e-3/1e-4); all 5500 analytic points lie inside their dps-35 interval enclosure.
- Dense scan and V6 agree to ~2e-5 relative, so the true infimum is 4.5515e-04 / 1.3998e-04; both the supplied value and my V1 value are conservative lower bounds on it.

Why the two implementations differ: raw interval enclosures of q/Den are loose at the low-eps end. Box `[1.5,1.501]`, x=0.2 gives `q in [6.27254306194357e-04, 9.00501369381922e-04]` while the true q range there is `[7.59249e-04, 7.60574e-04]` (~72x wider than the true variation). `q_lower`, `q_upper`, `denominator_lower` reproduce to <=1.6e-16 only because both codes share the same loose enclosure, not because it is tight.

## Part B — independent reproduction (exact steps at mp.dps=50)

1. **Scale c**: `c_1 = 1.2499999999999992939`, `c_2 = 1.2499999999999992753`; `c-1.25 = -7.06e-16 / -7.25e-16`; `gain2 - gain1/c = -4.52e-16`. **Reproduced** (~1e-15; 1.25 is exact, the residual is the 16-digit rounding of `world2_eps`).
2. **Section 5 identity**, two ways, per channel: `|delta_mu|` direct vs closed form agree to rel `1.70e-12` / `2.57e-12`; `delta_mu_1 = -4.4435604897130745e-07 - 1.3329430520715469e-09 i`, `delta_mu_2 = -2.0358066259898024e-07 - 4.133504510466580e-10 i`; `|delta|/|g t| = 3.33301e-04 <= |c-1| q = 3.33301e-04` (tight to 1.4e-6 rel). **Reproduced.** Cross-check `a1(xi=-1) = -t` exactly (`8.7e-19`) and `q = i a1/(1-a1)` reproduces `q` to 2.4e-16 rel.
3. **sigma = 0.2 % of |g1 t1|**: implied factor `0.001999999999999996` / `0.001999999999999995` (rel 2.1e-15 / 2.4e-15). **Reproduced.**
4. **D and p_***: `D = 0.20123662403751861` vs supplied `0.20123662403751783` (abs 7.8e-16, rel **3.86e-15**); `p_* = 0.44342318925625731` vs supplied `0.4434231892562575` (abs -1.9e-16, rel **4.25e-16**). **Reproduced.**
5. **Monte-Carlo, Theorem 4** (seed 20260910, N=4e6 per world, known-covariance quadratic discriminant): mean error `0.44360825`, 95 % CI `[0.44326398, 0.44395252]`; world 0 `0.44389175`, world 1 `0.44332475`; counts 1775567 / 1773299. `Phi(-D/sqrt2) = 0.44342319` **inside the CI**. Empirical whitened distance from sample means `0.20071240` (rel 2.6e-3, matching the 4e6-sample error). **Reproduced.**
6. **Section 7 bounds** with the supplied `m_i`: `0.0493981482 / 0.0926437516` vs supplied `0.0493980458 / 0.0926437138` (rel 2.1e-6 / 4.1e-7; the gap is the 4e-21 rounding of `b_a` in the supplied file). **Reproduced.** With my more conservative `m_i = (4.2929719577e-04, 1.2943650641e-04)`: `0.0497804328 / 0.0936670266`, growth **+0.774 % / +1.105 %**; **both stay below 0.10**. ✔ Union-bound failure probability `0.003` => coverage lower bound `0.9970000000000001` == supplied 0.997. **Reproduced.**
7. **Monte-Carlo, section 7 estimator** (seed 20260911, N=2e6 per channel-world, chunked): good-event fraction `0.9980015 / 0.997967 / 0.997997 / 0.9980675` for (ch1 w0/w1), (ch2 w0/w1). Each simulated unit has ONE channel eta plus the shared xi, so the correct per-unit theoretical lower bound is `(1-1e-3)^2 = 0.998001` (ch1/w0 matches to 5e-7); the supplied `>=0.997` is the joint 3-event union bound and is **met with margin**. Max `|eps_hat - eps|` on the good event `1.023e-02, 1.693e-02, 2.671e-02, 5.174e-02` against bounds `4.978e-02, 4.978e-02, 9.367e-02, 9.367e-02` => **0 violations**. Unconditional empirical max error `1.088e-02, 1.875e-02, 2.879e-02, 5.474e-02` — the excluded bad events produced no visible outlier at N=2e6.

## Not reproduced / cautions

- `qprime_lower` (-0.768 % / -1.092 %) and `amplitude_derivative_lower` (-0.774 % / -1.105 %) are the only supplied constants that do not reproduce; everything else in sections 5-7 reproduces to <= 4.3e-15. Both my value and the supplied value are conservative lower bounds on the true infimum (4.5515e-04 / 1.3998e-04); the difference is interval-arrangement tightness, not a disagreement about the mathematics.
- `code/mie_reactance.q_mie()` (and `q_mie_both_signs`) uses the prompt-literal conversion `-i a1/(1+a1)`, which is NOT equal to q (relative 2.0 for xi=-1, 2.7e-3 for xi=+1). Use `q_closed`/`mp_q`/`iv_quantities`, or `i*a1/(1-a1)` with `xi_sign=-1`.
