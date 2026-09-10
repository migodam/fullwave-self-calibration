# FINDINGS — increment 3 (Part C fixed-loss Born §3, Part D reference GLS §9, Part E a1 sign)

Interpreter `/Volumes/.../a3_research/.venv3d/bin/python`, numpy 2.5.3, mpmath 1.3.0 via
`PYTHONPATH=/Users/migodam/.cache/uv/archive-v0/VjSSQm31390egyfy`. Fixed seed **20260910** in every JSON.
Wall: C 20.7 s, D 3.5 s, E 0.2 s. No installs, single thread, no figures, `inputs/` untouched.

Artifacts (authoritative)
- Part C: `results/fixed_loss_born_20260910T100202Z.json` + `logs/fixed_loss_born.log`
- Part D: `results/ref_gls_20260910T100641Z.json` + `logs/ref_gls.log`
- Part E: `results/mie_sign_convention_20260910T100913Z.json` + `logs/mie_sign_convention.log`
- Superseded, kept for provenance (do not cite): Part C debug `…T095919Z / T095944Z / T100011Z / T100048Z / T100122Z.json`;
  Part D `ref_gls_20260910T100317Z.json` (grid-objective bug, below); Part E `mie_sign_convention_20260910T100743Z.json`
  (mis-coded mpmath/Hankel cross-checks).

## Part C — A5 §3 fixed-loss Born (B: C^{m×2}, ell = (0.03,0.05), 200 draws)
1. **det M(w) = -|g|^2 det[u,ell]**: VERIFIED. Max rel err 1.56e-13 (m6 well), 1.82e-13 (m12 well),
   8.30e-10 / 1.89e-10 (cond(B)=1e4). Direct check: det M(w)=+0.0408921 with det[u,ell]=-0.0408921, and
   det M + |g|^2 det[u,ell] = 1.4e-17 — the identity is exact, the *sign* just follows det[u,ell].
2. **Exact reconstruction** (w = B^dag y, M(w)[a,b] = ell, h = a+ib, u = Re(hw), g = 1/h):
   max rel |du|/|u| and |dg|/|g| = 1.56e-13 (m6 well), 1.82e-13 (m12 well), 8.30e-10 (m6 cond 1e4),
   1.89e-10 (m12 cond 1e4). Explicit inverse == 20-random-start nonlinear LS to 2.4e-14 (well-cond).
3. **Uniqueness / ambiguity**: 20 000 random (u,g) perturbations of size 1e-3 give min ||delta y|| = **3.42e-04**
   (median 5.25e-03) — strictly positive, so generic perturbations change the data. ell=0: (c*u, g/c) for
   c = 0.3, 2, 7 reproduces y to rel 1.6e-16 ⇒ real scale ambiguity confirmed.
4. **Parallel counterexample** (u=30*ell, u'=45*ell, c=(45+i)/(30+i), g=1, g'=1/c):
   g' = 0.6668311944718659 + 0.007403751233958541i vs quoted 0.666831+0.007404 (abs 3.16e-07, rel 4.7e-07,
   matches to 1e-6); same y to rel **2.00e-16 (< 1e-15)** over 9 independent random B. eps worlds (1.9,2.5) and
   (2.35,3.25). |g| = 1.0 inside gain_allowed_modulus [0.75,1.25]; |g'| = 0.666872 is OUTSIDE — legitimate only
   under the shared rescale lambda = 1.2 in [1.1247,1.25], which puts |g| = 1.2 and |g'| = 0.80025 both inside
   (residual 2.29e-16). So the raw pair is a data-equivalence fact, the rescaled pair is the legal counterexample.
5. **Near-parallel u(s) = q*ell + s*perp(ell)**: log-log slope of measured ||delta u|| vs |D| = **-1.0052**
   (B=I), **-0.9978** (random m6), **-0.9986** (cond 1e4) ⇒ error grows ~1/|D| as predicted. Max measured/bound
   ratio: B=I random 0.9987, adversarial delta-h **1.00000016** (overshoot 1.6e-7 = numerical search resolution:
   the bound is ATTAINED, not violated), random m6 0.6512, cond-1e4 m6 0.2205. All <= 1 to within 1.6e-7.
   mu = sigma_min(M(w)) matches the analytic prediction to max rel 6.8e-8.
6. **Local lower bound**: lambda_min(H) = (1-sqrt(1-4D^2/S^2))/2 VERIFIED over 20 000 draws, max abs err
   **3.17e-16** (mpmath dps=50). sigma_min(A_vis) >= sqrt2|g|sigma_min(B)sqrt(lambda_min(H)):
   B=I ratio in [1-1.14e-11, 1+6.4e-12] ⇒ **TIGHT**; random m8 [1.0023,1.2523] med 1.066; m12 [1.052,1.480]
   med 1.375; cond-1e4 m8 [1.042,2.072] med 1.166 ⇒ **SLACK by up to 2.07x**. Labelled probe: adding one extra
   geometry-nuisance direction drops the ratio to 0.9907, i.e. the bound can fail once geometry is nuisance
   (exactly the caveat in the text).

## Part D — A5 §9 reference GLS (m=5, sigma=0.4, sigma_r=0.35, 5000 noise draws)
7. **Scalar identity VERIFIED** (real and complex z): L(g*) == RHS(Woodbury) == RHS(direct solve) == augmented GLS
   to rel 1.34e-16 / 0.00e+00; g_* formula == dense-grid argmin. **No log-determinant**: adding
   log det(sigma^2 I + sigma_r^2 f f*) moves theta-hat by max |dtheta| = 2.24e-02 (median 9.8e-04).
   **Annulus break**: |g*| = 2.6439 outside [0.75,1.25]; L_uncon = 2.8774e+02 vs L_con = 1.2347e+03 (rel 3.291).
8. **Profiling equivalence** (single draw): theta_profiled = 0.613447871725, joint VarPro = 0.613447872570,
   augmented GLS = 0.613447873088 ⇒ rel 1.38e-09 / 2.22e-09. Plug-in theta = 0.46010 (diff 0.153).
9. **Null calibration** (T = 2 x concentrated weighted residual): at the TRUE theta, T ~ chi^2_{2m} exactly
   (mean 10.048 vs 10, KS p = 0.53). Profiled: T ~ chi^2_{2m-1} (mean 9.050 vs 9, KS p = 0.20).
   Gate rejection at the 95% chi^2_{2m-1} quantile (16.919): correct reference GLS **5.48%** (nominal 5%,
   MC se 0.31%) vs naive plug-in **78.18%** at the same threshold (and 75.7% at chi^2_{2m}). Plug-in mean
   T = 47.92 (5.32x dof), q95 = 131.2 vs 16.9. Estimator quality: RMSE 0.0604 (GLS) vs 0.1256 (plug-in),
   ratio **2.08x**; median |err| 0.0405 vs 0.0733. Mean-bias difference gls-plugin = -4.1e-04 +- 1.6e-03
   (not resolvable at N=5000) — the damage is variance/gate inflation, not a mean shift.
10. **Full residual derivative VERIFIED**: max rel err vs mpmath dps=60 central FD (h=1e-15) = **1.53e-31**
    (limited only by the h^2 truncation of the reference), max |v^*r| = 5.3e-60. Dropping
    -v((dv)^* r)/(v^* v) gives max rel err **0.9871** ⇒ the term is NOT negligible. Q_v = I - vv^*/(v^*v)
    (the orthogonal projector onto the complement of v) is the one used; the "2I - vv^*/(v^*v)" reading is
    not a projector and is excluded by the derivation in DERIVATIONS.md §8.

## Part E — a1 -> q sign convention (eps=2, x=0.2)
11. **Branch identified**: z*h1^(1)(z) == psi1 - i*chi1 == 1.4e-16, i.e. z*h1^(1) is the code's
    `mie_a1(..., xi_sign=-1)`. The code DEFAULT (xi_sign=+1) is psi1 + i*chi1 = z*h1^(2) (ingoing for
    e^{-i omega t}). Verified with an exact closed form -e^{iz}(z+i)/z, mpmath, scipy, and hankel1 routes.
12. Textbook a1 (h1^(1), e^{-i omega t}) = **1.7774330813186167e-06 - 0.0013332028810538396i** (4 routes agree
    to 6.5e-16); Re(a1) + |a1|^2 = **3.5548661626372334e-06** and Re(a1) - |a1|^2 = -2.1e-22 (lossless
    optical theorem); t = i q/(1-i q) = -1.7774330813186243e-06 + 0.0013332028810538424i = **-a1** (rel 2.1e-15).
    q(reactance) = 0.0013332052507369593.
13. Candidate conversions on that a1: `-i*a1/(1+a1)` = -0.001333195772055032 - 3.5548e-06i (**rel 2.0 vs q — wrong**),
    `+i*a1/(1-a1)` = 0.0013332052507369564 (**rel 2.1e-15 vs q — reproduces q**). So **candidate B is correct**;
    candidate A is the alternating (wrong-sign) geometric series and errs at O(|a1|^2) = 1.8e-06 (rel 2.7e-3).
    `code/mie_reactance.q_mie()` default (xi_sign=+1, form A) returns 0.001333195772055032 - 3.5548e-06i,
    rel 2.67e-03 vs q. On that conjugated branch the correct form is `-i*a1/(1-a1)` (= q, rel 2.1e-15).

## Not reproduced / bugs found and fixed
- **Part D `J1_grid` bug (fixed)**: the Woodbury numerator was computed as f^*y instead of f^*(y - z f), so the
  profiled grid minimised a different objective and landed at theta = 0.5333 instead of 0.6134 (15% spurious
  mismatch vs VarPro). After the one-line fix all three profilers agree to ~1e-9. The old JSON T100317Z must
  not be cited.
- **Part E first run (superseded)**: my mpmath chi1' had a sign error and the hankel1 derivative used an
  underflowed central complex step, so the "routes" disagreed by rel 2.0 / 1.5e3. Rewritten with exact closed
  forms; all routes now agree to 6.5e-16.
- **Part C ill-conditioned LS**: with cond(B)=1e4 the *naive unweighted* nonlinear LS loses precision
  (max rel diff 6.2e-2, only 0-4 of 20 starts reach the global optimum). The explicit inverse itself stays
  exact to ~1e-9. This is a scaling/solver artifact, not a failure of the §3 map.
