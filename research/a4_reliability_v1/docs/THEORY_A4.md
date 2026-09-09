# A4 theory: a restricted positive electromagnetic condition

Status: independently derived and numerically checked in this package. The algebra is not a claim of scientific priority. Closest-prior comparison remains a publication gate. No theorem below establishes arbitrary-shape microwave tomography or universal phase attribution.

## 1. Physical experiment and nuisance model

Use the time convention $\exp(-i\omega t)$ and a homogeneous background with real wavenumber $k>0$. The target is an isotropic, nonchiral, concentric radially stratified dielectric supported in $|x|\leq a$, with known center. The material can vary radially and may have loss; the full Maxwell scattering problem is assumed uniquely solvable. A nonzero electric $l=1$ scattering coefficient is required at each used frequency. This is not a weak-contrast assumption. The known center is a physical spatial anchor, not a coordinate trick.

Illuminate the target with three calibrated regular electric dipole modes:

$$
E^{\rm inc}_p(x)=\left(I+k^{-2}\nabla\nabla\right)j_0(k|x|)p,
\qquad p=e_1,e_2,e_3.
$$

An equivalent realization is the transverse Herglotz integral

$$
E^{\rm inc}_p(x)=\frac1{4\pi}\int_{\mathbb S^2}(I-dd^T)p\,e^{ikd\cdot x}\,dS(d).
$$

Thus three arbitrary incident plane waves are NOT the assumed acquisition. Finite-aperture modal synthesis has additional error and cost; neither is certified in these experiments. Relative source amplitudes/phases and the three receive-vector components are calibrated. A single unknown complex electronics coefficient per complete $3\times3$ response, independently at every position and frequency, is allowed. Arbitrary componentwise or entrywise gains are not allowed by the positive theorem.

Spherical symmetry makes the electric $l=1$ scattering operator scalar on its three-dimensional irreducible mode space. Regular electric $l=1$ excitation therefore produces exactly an outgoing electric $l=1$ field, even when internal scattering is strong. This classical modal fact supplies the physical restriction; it is not a SOM numerical reduction.

For a receiver $r=Rn$, $R>a$, $|n|=1$, the response is

$$
Y(r,k)=c(k)\frac{e^{ikR}}R M(z,n),\quad z=kR,
$$
$$
M=a_0I+b_0nn^T,\quad
a_0=1+\frac{i}{z}-\frac1{z^2},\quad
b_0=-1-\frac{3i}{z}+\frac3{z^2}.
$$

The target's exact modal coefficient, electronics gain and common delay are absorbed in $c(k)$. All allowed material derivatives lie in the same real two-dimensional span as complex scalar-gain derivatives. Profiling that span therefore also profiles the admitted material nuisance. This is the critical reason the theorem is nonvacuous and also the reason it is restrictive.

Observe nine complex entries with independent proper complex noise $\mathcal{CN}(0,\sigma^2)$. Realification multiplies real and imaginary parts by $\sqrt2/\sigma$. Set $q=\|Y\|_F/\sigma$, the total Frobenius signal-to-noise amplitude ratio, not a per-entry SNR.

## 2. Theorem 1: exact material/gain-profiled geometry spectrum

Let $B$ be the real-whitened Cartesian receiver-translation Jacobian. Let $N$ include all admitted material derivatives and the two real scalar-gain derivatives. Then the singular values of $(I-P_N)B$ are

$$
s_t=s_t=\frac{\sqrt2q}{R}
\sqrt{\frac{z^4+3z^2+9}{z^4+z^2+3}},
$$
$$
\boxed{s_r=\sigma_{\min}((I-P_N)B)
=\frac{2qkz\sqrt{z^2+4}}{z^4+z^2+3}>0.}
$$

The two tangential directions are perpendicular to $n$; the radial direction is $n$. Geometry is locally identifiable despite unrestricted frequency-dependent scalar gains and unrestricted radial material nuisance, provided $q>0$ and the calibrated tensor acquisition is available.

### Proof

The nuisance span contains $Y$ and $iY$. Differentiating the common factor $ce^{ikR}/R$ contributes only to this span and disappears on profiling. Work in coordinates with $n=e_3$, so $M=\operatorname{diag}(a_0,a_0,s_0)$, where

$$
s_0=a_0+b_0=2/z^2-2i/z.
$$

A tangential Cartesian unit displacement $v$ yields $b_0(vn^T+nv^T)/R$. It is Frobenius-orthogonal to $M$ and to radial derivatives, with squared norm $2|b_0|^2/R^2$. The two tangent derivatives are mutually orthogonal.

The radial derivative is $k\operatorname{diag}(a_0',a_0',s_0')$. Its complex-orthogonal projection away from $M$ has squared norm

$$
\frac{2k^2|a_0s_0'-s_0a_0'|^2}{H},\qquad
H=2|a_0|^2+|s_0|^2.
$$

Direct computation gives

$$
H=\frac{2(z^4+z^2+3)}{z^4},\quad
|b_0|^2=\frac{z^4+3z^2+9}{z^4},\quad
W=a_0s_0'-s_0a_0'=\frac{2i(z+2i)}{z^3}.
$$

Since $|ce^{ikR}/R|/\sigma=q/\sqrt H$, real-whitening produces the claimed values. For $x=z^2$, the denominator minus numerator in $s_r^2/s_t^2$ is

$$
(x^2+x+3)(x^2+3x+9)-2x^2(x+4)
=x^4+2x^3+7x^2+18x+27>0.
$$

Thus $s_r<s_t$. Rotational covariance extends the calculation to arbitrary $n$. This proves the theorem.

### Corollary 1: there is a restricted optimal range, not monotone near-field superiority

At fixed $k$ and fixed total relative SNR $q$,

$$
s_r\sim\tfrac43qk^2R\quad(kR\downarrow0),\qquad
s_r\sim\frac{2q}{kR^2}\quad(kR\to\infty).
$$

The maximum is attained at $kR=1$, where $s_r=2qk/\sqrt5$. Indeed, after squaring and setting $x=z^2$, the derivative of $x(x+4)/(x^2+x+3)^2$ equals

$$
-\frac{2(x+6)(x-1)(x+1)}{(x^2+x+3)^3}.
$$

This optimum is feasible only when the receiver can remain outside the target. If $ka\geq1$, the unconstrained optimal radius lies inside the target. Fixed transmitter power and fixed absolute noise change $q$ with range and can change the optimal design entirely. A family allowing $R\downarrow0$ also requires the target support to remain inside $R$; the asymptotic is not a prescription to put a sensor inside a fixed sphere.

### Corollary 2: gain sharing and repeated blocks

Independent blocks add profiled information when their scalar gains are independent. For a common translation, the minimum eigenvalue is at least the sum of their squared $s_r$ bounds, evaluated at the respective receiver-target vectors. Constraining gains to be shared reduces the nuisance space, so cannot reduce local information on the same data. This monotonicity does NOT prove material identifiability or removal of every clock/translation alias.

## 3. Theorem 2: complete noiseless geometry ambiguity and a physical acquisition remedy

Under Theorem 1, a nonzero tensor at one frequency determines $R$ and the unoriented axis $nn^T$. The complete geometry ambiguity is precisely $r$ versus $-r$.

### Proof

Since $\operatorname{tr}M=2$, scalar normalization gives

$$
Z=Y/\operatorname{tr}Y=M/2.
$$

Its double transverse eigenvalue is $a_0/2$. Its simple longitudinal eigenvalue is

$$
\lambda_L=s_0/2=z^{-2}-iz^{-1}.
$$

The eigenvalue is simple because $b_0\ne0$ for every positive real $z$. Consequently

$$
z=-1/\operatorname{Im}\lambda_L,\qquad
\operatorname{Re}\lambda_L=(\operatorname{Im}\lambda_L)^2.
$$

Its eigenspace determines $nn^T$. Both signs give identical $M$, and no other positive radius or axis can give the same normalized tensor. This proves completeness for this special model, not for general coherent inversion.

Now repeat the measurement after a known nonzero receiver displacement $d$, keeping both positions outside the target. If two receiver locations fit both tensors, the first implies $r_2=\pm r_1$. For the negative alternative, the second requires

$$
-r_1+d=\pm(r_1+d).
$$

The plus sign would require $r_1=0$; the minus sign would require $d=0$. Both are excluded. The second acquisition therefore removes the exact twofold ambiguity. Its gain may be independently unknown. A finite-noise design still needs a lower signal bound and sufficient separation; nonzero $d$ alone is not a uniform stability guarantee.

A calibrated electronics scalar reference cannot break $D(r)=D(-r)$: it leaves that geometric symmetry intact. It can instead help identify material through the scalar coefficient. The type of reference must match the ambiguity.

## 4. Proposition 3: material identifiability is a different task

With one unrestricted complex $g_f$ per frequency,

$$
c_f=g_f t_1(\alpha,k_f),
$$

any two material parameters having nonzero $t_1$ can be interchanged by replacing $g_f$ with $g_f t_1(\alpha,k_f)/t_1(\widetilde\alpha,k_f)$. Material is unidentifiable, even though Theorem 1 gives a strictly positive geometry bound. Adding more frequencies with new unrestricted gains does not fix it.

Independent electronics observations $z_f=g_f+\eta_f$ constrain exactly these complex scale directions. At finite reference noise the gain parameters are not literally deleted: the joint data/reference Jacobian is used and the uncertainty is propagated. A fresh 16-scene reference experiment tests this mechanism with known geometry. It does not establish simultaneous global material/geometry recovery.

## 5. Proposition 4: task-specific model fidelity

Suppose two radial-target solvers produce $t_1^{(c)}D(r)$ and $t_1^{(f)}D(r)$ with nonzero scalar coefficients. After free complex-gain profiling their geometry objectives are identical. Refining only the interior scalar coefficient cannot improve that profiled geometry objective. It can improve quantitative material recovery when gains are constrained or independently referenced. This is a physical use of classical variable projection, not a new variable-projection algorithm.

In contrast, replacing the exact exterior dyad by a radiative or static approximation changes the geometry-carrying tensor structure. With $D_\mathrm{rad}$ retaining $I-nn^T$ and $D_\mathrm{stat}$ retaining $(-I+3nn^T)/z^2$,

$$
\frac{\|D-D_\mathrm{rad}\|_F^2}{\|D\|_F^2}
=\frac{3(z^2+1)}{z^4+z^2+3},
$$
$$
\frac{\|D-D_\mathrm{stat}\|_F^2}{\|D\|_F^2}
=\frac{z^4+3z^2}{z^4+z^2+3}.
$$

Both approximations have exactly zero scalar-profiled radial information. Thus a small relative field error can coexist with complete loss of one task direction. These exact approximation-error identities apply only to these specified kernels under the modal model. They do not enclose unknown target asymmetry, modal-synthesis error, antenna coupling or DDA discretization error. The identities were checked in a separately labeled post-frozen algebraic supplement; no benchmark was tuned.

## 6. Counterexamples and boundaries

**Electronics:** At $n=e_3$, unknown independent receive-component gains absorb the diagonal radial change, so the minimum visible singular value becomes zero. Independent gains for all nonzero tensor entries remove all geometry information. The relative channel calibration is an essential resource.

**Excitation/projection diversity:** Three calibrated modal illuminations and three components are sufficient, not proven minimal. Some generic one-illumination vector measurements have full local rank, while axial configurations lose rank. Known orthogonal rotation of a complete measured vector is an invertible unitary change of data and does not manufacture information. Additional scalar projections may help only when the gain relations between them are controlled.

**Target class:** Unknown target translation reinstates an anchor ambiguity. Ellipsoids, boxes and anisotropic targets generally mix modes and do not satisfy scalar material factorization. The frozen DDA cases are out-of-class stress tests, not proofs of extension. Even an apparently small field mismatch cannot justify continuity of the profiled spectrum: a tiny new material derivative can add a new nuisance direction after column normalization. A valid perturbation theorem needs control of nuisance subspaces, not merely $\|F_\epsilon-F_0\|$.

**Signal zeros:** A vanishing electric modal coefficient removes all available information. No uniform lower bound over all dielectric parameters is possible without a signal lower bound.

**Coverage:** Local full rank does not distinguish $r$ and $-r$. The algebraic reconstruction does not remove finite-noise uncertainty, model inadequacy, or an incorrectly declared search domain.

## 7. Conditional covered-region theorem

Let $y_b$ be independent complex-whitened nine-entry blocks and $f_b(r)$ the unit-norm vectorized dyad at $r+d_b$. Let $P_b=f_bf_b^*$ and

$$
\ell(r)=\left(\sum_b\|(I-P_b(r))y_b\|^2\right)^{1/2}.
$$

Suppose the true geometry lies in the declared domain $\Omega$, every block follows the admitted modal family up to an actual deterministic whitened error of joint norm at most $\beta$, and the noise is standard proper complex Gaussian. For $m$ complex observations use

$$
\tau=\sqrt{\chi^2_{2m,1-\alpha}/2},\qquad
S=\{r\in\Omega:\ell(r)\leq\tau+\beta\}.
$$

Projection is contractive, so the true geometry belongs to $S$ with probability at least $1-\alpha$. This event is uniform over the search procedure, avoiding a presupposed finite candidate bank containing truth.

For a rectangular cell $C$ centered at $r_c$, let $\rho$ be its half-diagonal and let $R_{b,\min}$ be its minimum distance after shift $d_b$ to the target center. The derivative calculation in Theorem 1 implies

$$
\|DP_b(r)[v]\|_2\leq\sqrt3\|v\|/|r+d_b|.
$$

Integration along each cell segment, and $\|P-Q\|_2\leq1$ for rank-one orthogonal projectors, yields

$$
\operatorname{LB}(C)=\max\left\{0,\ell(r_c)-
\left[\sum_b\|y_b\|^2\min\left(1,\frac{\sqrt3\rho}{R_{b,\min}}\right)^2\right]^{1/2}\right\}
\leq\inf_{r\in C}\ell(r).
$$

A zero minimum distance is treated conservatively by using variation bound one. Cells are pruned only when this lower bound exceeds $\tau+\beta$. Unprocessed cells must remain in the outer set when a budget is exhausted. If the retained outer set lies inside a ball of task tolerance $\varepsilon$ around a feasible estimate, accepting that estimate has false-accept probability at most $\alpha$, conditional on the physical assumptions and exact bound arithmetic. An independent validation gate can only reduce this acceptance set.

**Implementation boundary:** The supplied code uses double precision and a small padding, not directed-rounding interval arithmetic. It is an implementation of an analytical bound, not a machine-verified certificate. The nonlinear material profiling in a general scattering problem is not solved by this result. A local material optimizer gives an upper bound, not a valid cell lower bound.

## 8. An exact obstruction that the controller cannot diagnose away

If the admitted discrepancy family contains $F(x_1)-F(x_0)$, the worlds

$$
(x,e)=(x_0,F(x_1)-F(x_0))\quad\hbox{and}\quad(x_1,0)
$$

have identical data distributions. When this equivalence persists over available validation and acquisition channels, no data-only selector can distinguish them. A method that accepts the clean second world must accept the first with the same probability. For an estimator, the larger squared task risk is at least $\|T(x_1-x_0)\|^2/4$. Abstention can avoid claiming correctness but cannot recover the missing attribution. The A4 two-world control enforces this exact equality at a 25 mm displacement, exceeding the frozen 15 mm task tolerance.

This is not a new A4 impossibility theorem; it is the operational boundary inherited from A3 and deliberately tested. A universal six-action reliable controller under unrestricted discrepancy is therefore not a supportable claim.
