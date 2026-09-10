# Parent derivation: the reference theorem needs an implementation-error budget

This is a corollary of the supplied A5 restricted-class inverse bound, not a new
electromagnetic novelty claim. No extension to unknown-geometry simultaneous
sphere imaging is asserted.

Let h_i(epsilon)=|t_i(epsilon)|, h_i'≥m_i>0, h_i≤H_i on the declared intervals.
Observed modal data are Y_i=g t_i+e_i+n_i and a real gain-modulus reference is
Z=|g|+e_a+n_a. Suppose |g|≥a_min, |e_i|≤b_model,i, |e_a|≤b_ref,bias,
and on a specified noise event |n_i|≤b_noise,i, |n_a|≤b_ref,noise.
Put B_i=b_model,i+b_noise,i and B_a=b_ref,bias+b_ref,noise<a_min.
Reverse triangle inequality, division by Z≥a_min-B_a, clipping, and the inverse
Lipschitz bound give

$$
|\widehat\epsilon_i-\epsilon_i|
\le \frac{B_i+H_iB_a}{(a_{\min}-B_a)m_i}.
$$

For target tolerance delta, this sufficient certificate requires

$$
b_{\rm model,i}+(H_i+\delta m_i)b_{\rm ref,bias}
\le \delta(a_{\min}-b_{\rm ref,noise})m_i
-b_{\rm noise,i}-H_i b_{\rm ref,noise}.
$$

The right side is the remaining implementation-error budget. A negative value
means this particular certificate cannot achieve the target even before adding
model bias; it does NOT prove impossibility for every estimator. A positive value
is useful only with justified measurement/model/reference error bounds. An
adjacent-grid discrepancy is not automatically such a bound. Relative modal
calibration errors must be propagated into e_i before applying the inequality.

Likewise, for two Gaussian worlds in a common complex-whitened metric, arbitrary
deterministic error balls can reduce mean separation D to max(D-beta0-beta1,0).
The corresponding binary testing lower bound is Phi(-D_eff/sqrt(2)). This is a
least-favourable bound for the DECLARED error balls; without showing that the
physical errors admit those directions, it is not a new Maxwell impossibility
theorem. Conversely, a finite candidate distance is only an upper bound on the
minimum separation and cannot establish a positive uniform guarantee.

These deductions explain the next gate: quantify achievable modal/reference
precision, rather than only showing a nonzero material derivative or re-running
optimizers. The supplied interval constants remain inherited until the separate
independent reconstruction is audited.
