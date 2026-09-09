# Parent audit gates for Family 4

The automatically generated Family-4 task is largely consistent with the
theory, subject to the following binding interpretation and control rules.

1. **Shared nuisance.**  Frequency blocks must be stacked vertically with one
   common 18-dimensional pose increment.  Separate pose copies would test a
   different model and cannot support a frequency-diversity claim.
2. **Absolute versus normalized information.**  Adding a block is tested by
   PSD monotonicity of `K_IS` and `K_eff/K_SLAM`.  Generalized retentions need
   not be coordinatewise monotone because `K_IS` changes.  Do not call a
   retention decrease a contradiction of absolute monotonicity.
3. **Map direction reconstruction.**  If `w` is an eigenvector of
   `Q_A^T P_perp Q_A`, the corresponding map parameter direction is
   `u=V_A diag(1/s_A) w`, normalized so `u^T K_IS u=1`.  Using `w` directly as
   a 24-vector map coefficient is only valid if this transformation has been
   applied.
4. **Duplicate control.**  Use `[A;cA]`, `[B;cB]`, with factor
   `gamma=1+c^2`.  No-prior and `J_aug=gamma J` predictions are exact; fixed
   prior is intentionally not invariant.
5. **Noise/row scaling.**  Frequencies can have different raw amplitudes.  The
   result must state the whitening/noise convention.  A raw identity-noise
   stack establishes the theorem for that metric, but an empirical
   frequency-diversity magnitude may partly reflect block energy.  Include or
   defer a declared block-normalized/SNR-matched control; do not silently
   conflate it with the raw result.
6. **Equal-path trajectory confound.**  The proposed construction uses radii
   `0.8`, `1.6`, and `3.2` for the circle, 180-degree arc, and 90-degree arc,
   while the straight path ranges still farther from the target at its ends.
   Thus coverage angle, standoff distance, and per-pose signal strength change
   together.  Record for every path: continuous design length, open polyline
   length through sampled poses, closed length if closure is counted, min/max
   target range, and row-block norms.  A circle using `endpoint=False` has a
   full circumference only if the last-to-first closing segment is explicitly
   part of the path.
7. The tested trajectory ordering is therefore a reproducible scenario result
   or counterexample, not a causal theorem that angular coverage alone creates
   the order.  A same-standoff or per-pose-energy-normalized secondary control
   is needed before making that narrower interpretation.
8. Preserve any failure of the chosen `1e-4` separation threshold.  The
   algebraic statement is the shared-compensation residual; the number of
   directions moved by a particular frequency is scene-dependent.
