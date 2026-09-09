**[accepted]** Conditions: whitened/realified finite-dim, no
pose prior; pad extra angles to `pi/2` if `rank A>rank B`; at most `rank B`
differ from 1; ordinary eigenvectors may rotate broadly. Finite prior gives
weighted shrinkage: `D=BJ_X^{-1/2}`,
`I-B(B^*B+J_X)^{-1}B^*=(I+DD^*)^{-1}`, `R_X=Q_A^*(I+DD^*)^{-1}Q_A`, spectrum
in `(0,1]`, generally not `sin^2`. **[conditional]**

(f) `rank L_X<=rank B` **[accepted]**; `tr L_X=||C^{+/2}B^*A||_F^2`,
`lambda_max(L_X)=||C^{+/2}B^*A||_2^2` (`C=B^*B+J_X`); affected
(`#{theta_i<pi/2}`) vs destroyed (`#{theta_i=0}`) DoF are distinct.

(g) Along fixed-rank regular path `J_X(alpha)=alpha I` with support covering
every data-coupled pose direction: `alpha->0` recovers `K_SLAM`,
`alpha->infinity` recovers `K_IS`. **[conditional]** Singular/unanchored graph
gauges never converge to `K_IS`.

(h) Duplicate `(A_2,B_2)=c(A_1,B_1)` scales `K_IS`,`K_SLAM` by `1+|c|^2` and
leaves every `rho_i` unchanged **[accepted] only no-prior or jointly scaled
prior** (prior scaled by `1+|c|^2`). Caveat (req. 13): fixed finite `J_X`
changes data-to-prior weighting, so normalized retention generally changes.
