# Parent theory audit — A3

Status: mathematical audit, not a novelty certificate or experimental acceptance.

## A3 statements retained with exact scope

1. **Fixed-chart first-order completeness is correct.** Invertible C1 M and fixed full-column-rank complex U imply uniqueness of the residual-minimizing coefficients. At zero residual, differentiating C* C c=C* b gives C* C c_v=C*(b_v-M_v j). Exact state and declared tangent columns are obtained iff they lie in Range(U). This is a complex-linear fixed-chart statement, not a lower bound for a moving manifold or output-only interpolation. For multiple illuminations concatenate columns.
2. **Nonzero residual correction is mandatory.** With C=MU and r=b-Cc, C* C c_v=C*(b_v-M_v Uc)+(M_v U)*r. If U moves replace C_v by M_v U+M U_v, and j_v=U_v c+Uc_v. Frozen-chart numerical derivatives cannot silently represent moving-chart derivatives.
3. **Gain quotient claim is correct after transforming covariance.** For real raw noise Sigma and whitening W=Sigma^(-1/2), the cycle differential in whitened coordinates is L=L_raw Sigma^(1/2). Therefore L^T(LL^T)^dagger L=I-P_(W G_raw) when its kernel is exactly that gain range. Using untransformed cycle covariance is generally wrong. No claim of finite-noise statistical sufficiency is made.
4. **Independent-branch Gaussian bound is correct and conditional.** Condition on training and acquisition selection. Candidate predictions are frozen. A representative within beta of the true prediction must be present. Projection of new standard Gaussian noise onto each candidate difference gives the stated tail/union bound. A fitted residual is not a certified beta. Re-fitting on validation noise invalidates this derivation. K=1 gives no global coverage guarantee.
5. **Dual-residual identity and inverse bound are correct.** The correction is p_tilde* r and the remainder is r_d* M^-1 r. Its norm requires a verified inverse bound. The triangular augmented inverse is diag(M^-1,M^-1) plus one block -M^-1 M_v M^-1, giving gamma+gamma^2||M_v||. An offline exact singular value is not a cheap online certificate.
6. **Exact clock/gain and phase-center gauges are retained.** Free per-frequency phase absorbs delay. Two nonidentical planar rotations resolve the particular translation/phase-center gauge only if effective positions are themselves observable. Two 3D rotations leave their relative rotation axis; three noncoaxial orientations can remove that gauge. These statements do not establish full-wave global identifiability.

## Extension for sparse acquisitions: gain-graph obstruction

This is a graph/gain-quotient specialization of established incidence and closure ideas. Originality has NOT been established. It is useful because the A3 rectangular all-channel cross ratios do not directly cover missing Tx–Rx links.

Let E be acquired nonzero entries of H on a bipartite graph with active receiver and transmitter vertices. Let B have one +1 at the receiver and one -1 at the transmitter per edge. Changing the sign of transmitter log gains is just a reparameterization. Let n be the active vertex count and c the component count. The complex log-gain tangent is diag(h) B, with rank n-c, so the complex local gain-invariant dimension is |E|-n+c (twice that in real coordinates). This counts available invariant directions, NOT the rank of the geometry Jacobian after material elimination.

**Global tree obstruction.** On a forest, arbitrary nonzero h and h' are gain equivalent: for each component choose one root gain, then assign the remaining vertex gains successively so that a_r b_t=h'_e/h_e on each newly encountered edge. There is no cycle consistency equation. Thus any geometry/material change preserving nonzero edge values can be absorbed by free separable gains. No numerical basis, SOM or otherwise, can restore geometry information without adding assumptions or measurements. On a cyclic graph the signed products around fundamental cycles impose exactly the missing compatibility conditions.

**Acquisition consequence.** Adding a new leaf edge with its new independent gain yields no additional quotient dimension. Adding a chord between existing vertices in a component yields one complex cycle constraint. A two-edge batch through a new receiver between two already connected transmitters creates one cycle although either individual edge creates none. Therefore single-edge information-greedy selection can have zero marginal utility for both members of a useful pair. This is a topology-based screening rule, not proof that each cycle identifies geometry.

**Whitened implementation without ratios.** Form complex gain Jacobian Gc=diag(h)B, realify as [[Re Gc,-Im Gc],[Im Gc,Re Gc]], whiten by the raw covariance, and project physical parameter derivatives using a rank-revealing SVD. Then jointly profile material and other electronic nuisance. This avoids unstable divisions by small measured fields. A zero edge invalidates the nonzero graph dimension theorem; the raw Jacobian remains the operative numerical diagnostic. Discarding weak measurements is not part of the likelihood.

**Physical source of positive information.** Nonseparable full-wave transfer can produce nonconstant cycle invariants, but unknown material may still reproduce geometry derivatives. Hence a cycle-rich graph is necessary in the fully free-gain model but not sufficient. Frequency-shared hardware gains, reference measurements or constrained gain models change the graph/parameter model and must be declared.

## Algorithm adoption gates

- Graph screening can reject a structural zero-information design before expensive simulation; subsequent selection uses actual whitened physical tangents and independent branch separation.
- A reported smallest positive singular value must not hide missing target dimensions. Report full target rank, including zeros, in physically scaled coordinates.
- Gaussian validation guarantees require externally justified beta or a declared model-conditional beta=0. If unavailable, output empirical diagnostics rather than a certificate.
- No derivative test, graph argument or algebra identity proves SOM superiority, high-dimensional basin coverage, hardware calibration, or publication readiness.

## Parent far-field receiver/gain specialization

See FAR_FIELD_GAIN_THEOREM.md for the full conditional proof and C1 remainder assumptions. The leading (ik-1/R)(n dot h)H receiver translation derivative is a gain tangent when receiver complex gains are free across illuminations at fixed frequency. After profiling gains, the remaining raw derivative is O(R^-2), or O(R^-1) relative to field amplitude. Fixed relative-SNR Fisher information thus has an O(R^-2) upper bound at fixed k. This is not a universal high-k asymptotic, not a lower bound, and not applicable without modification to gains tied across frequency.

Ten independent vector Treams checks at k9/k18 and R1.3–20.8 support the mechanism: unprofiled singular values remain approximately constant at fixed SNR, while the smallest gain-profiled singular value scales with fitted distance slopes -1.010/-1.023. These fitted slopes are numerical observations, not the proof. Material was fixed in this test; allowing material nuisance cannot enlarge the visible norm. Nearest-prior novelty verification remains open.
