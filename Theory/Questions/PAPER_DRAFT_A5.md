# Material-Gain Attribution in Full-Wave Electromagnetic Measurements
## Modal Mechanisms, Finite-Noise Limits, and Cross-Discretization Failures

**Working manuscript, 10 September 2026. Not submission-ready.** This draft reports a restricted theorem and executed numerical results. It does not establish a general self-calibrating tomography algorithm, a state-of-the-art advantage, or the completion of the parent Q5 research task. The detailed theory, protocols, code, raw outcomes, and literature-access ledger accompany this manuscript.

### Abstract

Quantitative electromagnetic imaging with uncertain calibration requires distinguishing material changes from changes in electronics and geometry. We study this distinction for two independently parameterized dielectric regions. We first show that precisely known nonzero losses can themselves remove a Born-scale ambiguity, except on an explicitly characterized exceptional set; consequently, such information must not be credited to full-wave nonlinearity. Vector-Maxwell ablations then separate single-object internal response, electric/magnetic modal diversity, distinct-object response ratios, and inter-object rescattering. In the examined configuration, rescattering strengthens the weakest joint material direction but weakens a previously studied common-scale direction. A restricted exact result considers two separately measured, lossless spheres with a common unknown gain. Their Mie reactance coefficients generate finite material changes with only a radiatively small data difference. At specified noise levels, a concrete pair gives a 44.34% worst-case material-failure lower bound. One additional real gain-amplitude reference yields a uniform material-error bound on declared intervals. However, the main simultaneous-sphere numerical task remains unresolved: a reference-assisted method meets joint task tolerances in twelve same-formulation scenes but in only four of six independent-discretization scenes, with only one accepted by a residual gate. These findings identify useful mechanisms and limits without equating local observability with reliable material recovery.

**Index terms:** electromagnetic inverse scattering; material identifiability; gain calibration; Mie scattering; finite-noise ambiguity; model discrepancy.

## I. Introduction

A coherent scattered field may be fitted accurately while its material interpretation is incorrect. An unknown complex electronics factor can compensate for a changed material response; uncertain geometry introduces further alternatives. Retaining phase therefore raises a sharper question than whether phase is sensitive to position: which material-dependent changes cannot be reproduced by the admitted calibration variables, and are those changes large enough to survive the error budget?

Calibration and identifiability are established topics. Classical array self-calibration analyzes geometry and gain ambiguities [L8]. Amplitude-independent electric-dipole localization is also known [L5], and electromagnetic near-field localization has an extensive local-information framework [L6]. Antenna-aware diffraction tomography explicitly addresses incident-field and effective-height calibration [L9]. Recent quantitative radar imaging includes complex calibration factors within MFSOM-based reconstruction [L2]. General blind gain/phase calibration supplies structured bilinear identifiability results [L3]. None of these topics can be renamed as a contribution of the present work.

Likewise, vector spherical waves, T-matrices, and multiple-scattering equations are established computational tools [L4]. Reactance matrices and radiative corrections enforce familiar energy-conservation constraints [L1]. Approximation-error methods already address the discrepancy between accurate and simplified scattering models [L7], with modern connections to projection and prior-informed subspaces [L11]. Radiation-operator singular spectra concern a different unknown object from a nonlinear material derivative [L10].

The present candidate contribution has two parts. First, we connect explicit electromagnetic material constraints to finite-noise material/gain ambiguity and a narrowly defined minimal reference. Second, we provide a causal numerical audit showing which modal and multi-object effects survive calibration profiling, and where independent-discretization recovery fails. Neither part is presented as an established priority claim over all optical or microwave material-calibration literature. The accompanying ledger records the limits of the closest-prior review.

## II. Physical Models and Statistical Tasks

### A. Simultaneous-sphere task C

Use the time convention exp(-i omega t), unit relative background permittivity and permeability, and real background wavenumber k. Two known disjoint spheres have relative permittivity

$$
\epsilon_r(r)=1+\sum_{i=1}^2(\alpha_i-1+i\ell_i)1_{D_i}(r),
\qquad \alpha\in\mathbb R^2.
$$

The exact field satisfies

$$
\nabla\times\nabla\times E-k^2\epsilon_r E=0,
\qquad E=E^{inc}+E^{sc},
$$

with tangential E/H continuity and outgoing exterior radiation. The measured quantity is scattered electric field, not total field. With a common receiver-frame translation x and an electronics gain g,

$$
y=gF(\alpha,x)+e+n,\qquad n\sim\mathcal{CN}(0,\sigma^2I).
$$

Here e denotes a deterministic model discrepancy whose valid bound is not yet available for the main task. The gain satisfies a nonzero modulus constraint. An optional reference is z=g+n_r, with proper complex variance sigma_r squared.

The main experiment uses four fixed plane-wave/polarization illuminations and twelve three-component receivers. Moving these nonperturbing receivers changes observation but not the internal resolvent or the fixed illumination. Transmitter motion would require differentiating the incident field and is not simulated here. Internal polarization currents are constrained physical states, not independent nuisance currents.

One common gain is an acquisition assumption: relative illumination and receive-channel responses have already been calibrated, and remaining drift is approximately common over a switched measurement block. This condition has not been established for a real apparatus. Per-illumination and frequency-dependent gain families are considered as separate diagnostics, not silently identified with the common-gain model.

### B. Restricted modal task M

For the finite-noise theorem, two lossless homogeneous spheres are measured separately. Their geometry and modal normalization are known, and one normalized electric-dipole coefficient is recorded for each sample under one common complex gain. The absence of inter-object coupling follows from separate measurement, not an approximation to task C. The two permittivities are independently unknown on declared intervals.

This task is deliberately narrower than simultaneous imaging with unknown geometry. Its theorem cannot be transferred to task C without controlling the measurement operator, modal extraction, calibration drift, and geometric uncertainty.

### C. Task-specific profiling

For proper complex noise, real whitening is sqrt(2)/sigma times concatenated real and imaginary parts. Let A denote the two material columns and N the real-whitened geometry and electronics columns. Material information is assessed through

$$A_{vis}=(I-P_N)A.$$

When assessing one material separately, the other material column must also be included in the nuisance matrix. Geometry instead uses the material and gain columns as nuisance. A complex gain contributes two real columns. The squared singular values are local Gaussian information eigenvalues; they are not a substitute for finite separation or recovery.

The primary task tolerance is an absolute error of 0.10 in each real relative permittivity, 0.5 mm in x, and 5% in gain. Bare derivatives are not compared with field errors: a declared material step or tolerance scale is applied first.

## III. What Information Is Actually Added?

### A. Precisely known losses can already anchor a Born model

At known geometry, let B have two linearly independent complex columns and write

$$y=gB(u+i\ell),\qquad u\in\mathbb R^2,\quad g\ne0,$$

where u is real contrast and the loss vector ell is exactly known.

**Proposition 1.** If det[u,ell] is nonzero, the noiseless pair (u,g) is uniquely determined. If ell=0, admissible real contrast scaling remains ambiguous. If ell is nonzero but u=q ell, there is an exact exceptional family.

Equality of two data vectors implies u'+i ell=c(u+i ell), with c=a+ib. The imaginary part gives b u+(a-1)ell=0. Linear independence forces a=1 and b=0. On the exceptional line, c=(q'+i)/(q+i) provides the alternative material and compensating gain, subject to admissibility. Appendix A gives an explicit inverse and a finite perturbation bound.

Thus a common scale ambiguity in an unconstrained complex-contrast model cannot automatically be used for two real permittivities with fixed nonzero losses. In our numerical geometry, the Born model has a small second visible direction with the prescribed losses, whereas that direction vanishes to numerical precision at zero loss. This is prior information supplied by the material class, not information created by full-wave scattering.

### B. Full-wave response ratios, not merely additional modes

For one homogeneous sphere, all Born modal coefficients share the same contrast factor, even if the field contains several electric, magnetic, and angular components. At small size x=ka, exact nonmagnetic sphere coefficients in the present outgoing-plus-regular convention satisfy

$$
t_E=i\frac23\frac{\epsilon-1}{\epsilon+2}x^3+O(x^5),\qquad
 t_M=i\frac{\epsilon-1}{45}x^5+O(x^7).
$$

Uniformly on compact material intervals away from degeneracies,

$$
\frac{t_M}{t_E}=\frac{\epsilon+2}{30}x^2+O(x^4),\qquad
\partial_\epsilon(t_M/t_E)=\frac{x^2}{30}+O(x^4).
$$

The Born limit is a weak-contrast limit; the displayed asymptotics are a small-size limit. They must not be conflated. Electric and magnetic dipoles both have angular degree l=1. The useful full-wave effect is a nonproportional material dependence of their responses, not the claim that magnetic dipoles are absent from a Born field.

For data gVt, with known full-column-rank V and t=(t_E,t_M), a direct projection calculation gives

$$
\|Q_{Vt}Vt'\|\ge\sigma_{min}(V)
\frac{|t_E|\,|\partial_\epsilon(t_M/t_E)|}
{\sqrt{1+|t_M/t_E|^2}}.
$$

This is a nonempty known-geometry sufficient mechanism, but its value depends on modal readout conditioning and noise. It does not yet profile unknown geometry. Implementations use the original complex data, avoiding division by a small noisy coefficient.

### C. Distinct objects and rescattering are different effects

Even with one retained electric mode per object, the sum of two exact isolated responses may have a material-dependent response ratio. That structure is not evidence that inter-object rescattering is necessary. Conversely, adding rescattering changes illumination at each object and can strengthen or weaken a particular material combination.

The numerical controls therefore cross four axes: physical Born response, isolated-object response, angular/electric-magnetic truncation, and complete interaction. They also compare a single object with a coherent sum. A nonzero internal second-order field is not sufficient; its measured component must survive all admitted nuisance directions.

## IV. An Exact Finite-Noise Boundary and a Minimal Reference

### A. Exact Mie near-equivalence

For a lossless single mode, energy conservation in the adopted T convention implies Re(t)=-|t| squared. Away from t=-1, the classical reactance representation [L1] is

$$t(q)=\frac{iq}{1-iq},\qquad q\in\mathbb R.$$

On intervals where q_i(epsilon_i)>0 is strictly increasing, choose c>0 such that both transformed material values defined by q_i(epsilon_i')=c q_i(epsilon_i) and the gain g'=g/c remain admissible. Then

$$
g't(cq_i)-gt(q_i)
=-\frac{g(c-1)q_i^2}{(1-icq_i)(1-iq_i)},
$$

and consequently

$$
\frac{|g't(cq_i)-gt(q_i)|}{|gt(q_i)|}
\le |c-1||q_i|.
$$

These are exact identities for the nonlinear sphere response, not a truncated Born model. Since q scales as (ka) cubed for small spheres, the relative phase structure that excludes material scaling may be radiatively small.

This does not assert generic noiseless nonidentifiability. With measured coefficients h_i=gt_i, losslessness gives two real equations Re(h_i conjugate(g))=-|h_i| squared. For q_1 and q_2 distinct and nonzero, the coefficient matrix is nonsingular. Indeed its determinant, before multiplication by |g| squared, is

$$
\frac{q_1q_2(q_2-q_1)}{(1+q_1^2)(1+q_2^2)}.
$$

Thus noiseless uniqueness and severe finite-noise instability coexist. The case q_1=q_2 is a further exact degeneracy, not the source of the nondegenerate numerical pair below.

### B. Finite statistical risk

For two proper complex Gaussian distributions with equal covariance Sigma, define

$$D^2=(\mu_1-\mu_0)^*\Sigma^{-1}(\mu_1-\mu_0).$$

The optimal equal-prior binary error is Phi(-D/sqrt(2)). If the two required material-error sets are disjoint, every material estimator has failure probability at least this value in at least one of the two worlds. Appendix B supplies the testing reduction.

For ka=(0.2,0.15), material intervals [1.5,4] and [2,5], and gain modulus [0.75,1.25], consider

$$
\epsilon=(2,3),\quad g=1;
\qquad \epsilon'=(2.35974058,3.98931942),\quad g'=0.8.
$$

The alternative is defined exactly through the reactance inverse; displayed decimals are rounded. The material-vector distance is 1.05269, while relative data difference is 0.0303606%. Use the same absolute complex noise standard deviations in both worlds, namely 2.6664081e-6 and 1.8047986e-6. These equal 0.2% of the first-world modal amplitudes. Then D=0.2012366 and the failure lower bound is 0.4434232. The componentwise 0.10 material-error boxes are disjoint.

### C. A single real amplitude reference

Add Z=a+xi, where a=|g|>=a_min>0 and xi is real Gaussian. Define h_i=|t_i|=q_i/sqrt(1+q_i squared). Suppose, throughout each material interval, h_i<=H_i and h_i'>=m_i>0. Estimate

$$
\widehat\epsilon_i=h_i^{-1}
\left[\operatorname{clip}_{h_i(I_i)}(|Y_i|/Z)\right].
$$

On the event |n_i|<=b_i and |xi|<=b_a<a_min,

$$
|\widehat\epsilon_i-\epsilon_i|
\le \frac{b_i+H_i b_a}{(a_{min}-b_a)m_i}.
$$

For the stated intervals, interval evaluation of the exact Mie formula gives conservative derivative lower bounds greater than 4.32619e-4 and 1.30866e-4. The reference standard deviation is 0.001 in gain-amplitude units. Allocate probability 0.001 to each of the two complex-noise tails and the real-reference tail. The resulting simultaneous probability is at least 0.997, with material-error bounds below 0.04940 and 0.09265. These are uniform on the declared intervals under the fixed absolute noise model.

Within the class of interventions that add independent real scalar gain references while leaving the original data, noise, and target tolerances unchanged, zero references cannot meet a 1% worst-case failure target, whereas one suffices. Therefore one is the minimal number in that restricted class. This is not a minimum over all new frequencies, magnetic modes, repeated acquisitions, or hardware costs. It also requires a substantially accurate reference and does not establish practical modal extraction.

## V. Numerical Methods and Evaluation

### A. Independent formulations and physical controls

The new vector-spherical solver enforces each sphere's Maxwell interface conditions, projects regular/outgoing translations on angular quadrature, and solves s=(I-TU)^(-1)T a_inc. The inverse uses L=4, giving 96 complex outgoing coefficients for two spheres. A direct dyadic-Green first-Born volume integral uses spherical quadrature. Isolated and interacting responses share illumination, receiver positions, material freedom, and gain structures.

Independent data are generated with the supplied volume-integral electric-dipole discretization using Clausius-Mossotti cells and a radiative correction. It does not use sphere T coefficients. Both formulations nevertheless assume the same ideal electromagnetic geometry; no antenna-system modeling independence is claimed. Four DDA spacings are checked, and their differences are not treated as certified continuum bounds.

### B. Constrained gain profiling and references

At fixed material and geometry, the gain minimizes

$$L(g)=\|y-gf\|^2/\sigma^2+|z-g|^2/\sigma_r^2.$$

The interior solution is the usual scalar weighted least-squares expression; a modulus-annulus constraint is handled by radial clipping. Material and geometry are optimized with three fixed starts, bounded least squares, complete residual finite differences, and fixed stopping rules.

The reference-first comparator that correctly propagates reference noise has objective

$$
(y-zf)^*(\sigma^2I+\sigma_r^2ff^*)^{-1}(y-zf),
$$

which equals min_g L(g) at an interior gain. This equality is established algebraically and numerically. The simpler plug-in comparator fixes gain to z and omits this propagation; its excessive rejection cannot establish an advantage of a new optimizer. Appendix D records the exact residual derivative and equivalence.

### C. Frozen random protocol

Twelve same-formulation scenes use spectral L5 generation and L4 inversion. Six independent scenes use DDA generation at 0.009 m spacing and spectral inversion. The two real permittivities are independently sampled, as are gain and a small receiver-frame translation. Each method has three identical material starts. Proper complex noise is 1% in aggregate RMS; the optional complex reference has standard deviation 0.01.

Task success requires both absolute permittivity errors <=0.10, translation error <=0.5 mm, and relative gain error <=5%. A conventional 99% chi-square residual gate is reported only as a heuristic. Unused receiver positions assess prediction. A separate set of 36 fixed points at radius 0.16 m assesses the physical scattered E field without electronics gain. This discrete exterior structural field is neither the training residual nor an internal-domain field certificate.

## VI. Results

### A. Implementation checks

At the nominal configuration, L4-to-L5 field change is 5.17e-8, increasing translation quadrature changes the field by 5.08e-15, and refined physical Born quadrature changes it by 1.88e-14. Interface checks are modal algebra residuals, not an independent angular-grid boundary validation. Weak-contrast convergence toward Born and lossless/passive coefficient identities are verified. Twenty-one inherited A4 tests and nineteen new regression assertions pass.

The new solver also reproduces the six-scale trend from the supplied evidence: isolated-sphere sums retain slightly more of the old common-scale sensitivity than interacting clusters. This is an independent implementation check, not a new execution of Treams.

### B. Causal material spectra

At real permittivities (2,3), losses (0.03,0.05), one unknown translation and a common complex gain, field-normalized material singular values are:

| Model | Larger | Smaller |
|---|---:|---:|
| Physical Born | 0.427404 | 0.000870689 |
| Isolated sum, electric dipoles only | 0.316485 | 0.001590198 |
| Isolated sum, electric and magnetic l=1 | 0.306233 | 0.005395396 |
| Interacting, electric and magnetic l=1 | 0.306872 | 0.007308733 |
| Isolated sum, L=5 | 0.303429 | 0.005365905 |
| Interacting, L=5 | 0.304113 | 0.007281149 |

These values are per unit permittivity and are not yet noise-whitened. The weaker interacting material combination is approximately (-0.3713,-0.9285). Profiling the other material separately yields sensitivities 0.019575 and 0.007841 for the two components.

With zero loss the Born weak value is 2.39e-13. A single sphere restricted to one electric dipole remains material-invisible under its gain, whereas adding the magnetic dipole makes its material direction visible. Thus known loss, multiple objects, electric/magnetic diversity, and rescattering are distinct resources. Per-illumination gains reduce the isolated-to-interacting weak-direction increase from approximately 35.7% to 5.4%. No universal benefit of rescattering follows.

### C. Recovery and wrong attribution

| Data | Method | Task success | Wrong acceptance | Rejected |
|---|---|---:|---:|---:|
| Same formulation, 12 scenes | No reference | 11/12 | 1 | 0 |
| Same formulation | Joint reference | 12/12 | 0 | 0 |
| Same formulation | Plug-in reference first | 12/12 | 0 | 7 |
| Independent DDA, 6 scenes | No reference | 0/6 | 4 | 2 |
| Independent DDA | Joint reference | 4/6 | 0 | 5 |
| Independent DDA | Plug-in reference first | 4/6 | 0 | 5 |
| Independent DDA | Biased reference | 0/6 | 0 | 6 |

Only one of the six independent scenes is accepted by either correct-reference method. The seven same-formulation plug-in rejections use a gate that does not propagate reference noise; they must not be interpreted as a fair superiority result for joint profiling.

For independent data without reference, median absolute material errors are 0.333 and 0.479, but median translation error is only 0.075 mm. Median structural-field error is 15.1%. All 180 optimization starts report optimizer success, and the maximum material spread among the three starts of a fit is only 3.34e-6. Hence convergence and repeatable endpoints do not establish material correctness.

DDA field discrepancies against the spectral reference are 1.158%, 0.567%, 0.663%, and 0.678% at spacings 0.016, 0.013, 0.011, and 0.009 m. This trend is nonmonotonic. The finest discrepancy after gain profiling is 0.485%, larger than the approximately 0.0728% nominal linear signal from a 0.1-unit weakest material step. This is a warning about the error floor, not a certified decomposition of every failure.

### D. Finite candidates, frequency structure, and range

Eight constrained competitor searches are performed after the frozen recovery study. For example, a nominal material (2,3) and a feasible competing material (2.1,3.25521), with reoptimized gain and translation, differ by only 0.195% in sensor field. Such pairs upper-bound the nuisance-optimized separation. A large separation among the found pairs does not lower-bound all unsearched alternatives.

Three-frequency diagnostics at k=12,18,24 compare electronics classes using the same 432 complex sensor data. Independent frequency gains, shared gain, shared gain plus one delay, and affine complex frequency response give different local sensitivities. These are not same-budget comparisons with the 144-data single-frequency study, and no frozen multi-frequency recovery claim is made.

A fixed-absolute-noise range scan from 0.14 to 2 m yields decreasing local weakest material and one-dimensional geometry sensitivities at the sampled distances. No universal intermediate material optimum is observed or proved. This is not the fixed-relative-SNR radial geometry question of the earlier A4 model.

### E. Cost and statistics

The recovery script records 12.256 s wall time in a single-threaded CPU configuration, including 1.837 s recorded data generation and 9.789 s total multistart fitting. These are computational execution measurements, not total research or acquisition cost. Scene counts are too small for a general success-rate or speed-superiority claim. Reference construction, modal calibration, and real measurement duration are unmeasured costs.

## VII. Scope, Limitations, and Publication Gates

The main positive recovery claim fails to close. There is no independently justified continuum or hardware discrepancy envelope for task C, no covering lower bound on finite material separation, and no reliable end-to-end acceptance across independent data. Additional initializations do not repair an information/model problem.

The restricted task M does supply exact admissible worlds, a finite-risk theorem, explicit full-interval reference bounds, and a minimal reference count within a defined intervention class. However, separate samples, known geometry, losslessness, ideal modal access, and stable shared electronics are material assumptions. Their hardware implementation and the closest-prior originality comparison remain incomplete. This restricted result therefore does not fulfill all realism and scope requirements of the parent Q5 task.

The work does not restore SOM as a core contribution. It does not establish superiority over correctly weighted reference calibration, Bayesian approximation-error methods, or strong acquisition-design baselines. No measured data, unknown support, anisotropic material, general distributed reconstruction, or SLAM trajectory is reported.

## VIII. Conclusion

Full-wave nonlinearity can supply nonproportional modal and object-response dependencies absent from a homogeneous Born scale model. Whether these dependencies support quantitative recovery depends on the admitted material class, gain sharing, geometric nuisance, and error floor. Known loss can itself supply information; electric/magnetic diversity is distinct from higher angular orders; rescattering may have opposite effects on different material combinations. A restricted exact Mie example shows that noiseless material/gain identifiability can coexist with a substantial finite-noise failure lower bound, and quantifies how one amplitude reference changes that limit. The independent-discretization failures prevent a general positive method claim. The appropriate next step is a physically justified error and reference budget, not another local-information label.

## Appendix A. Explicit Born Inverse and Stability

Let w=B-dagger y=g(u+i ell) and h=1/g. Writing h=a+ib gives

$$M(w)(a,b)^T=\ell,\qquad M(w)=[\Im w,\Re w].$$

Its determinant is -|g| squared det[u,ell]. After solving the real 2-by-2 system, u=Re(hw). If the compressed error norm is at most eta and mu=sigma_min(M(w))>eta, then the perturbed matrix remains invertible and

$$|\delta h|\le\frac{|h|\eta}{\mu-\eta},$$

$$\|\delta u\|\le |h|\eta\left(1+\frac{\|w\|+\eta}{\mu-\eta}\right).$$

Indeed, the matrix perturbation has norm at most eta, and subtraction of the two linear systems gives the first bound; expansion of the product h w gives the second. A loss-reference error adds its norm to the numerator of the first bound. Unknown geometry requires additional treatment of B and is not included in this inverse.

For B=I and g=1, gain profiling gives the real material Gram matrix

$$H=I-\frac{uu^T+\ell\ell^T}{\|u\|^2+\|\ell\|^2}.$$

This follows by expanding the squared norm after projection from the complex line u+i ell. Its minimum eigenvalue is

$$\frac{1-\sqrt{1-4\det[u,\ell]^2/(\|u\|^2+\|\ell\|^2)^2}}2.$$

For general B, applying its smallest singular value before minimizing the gain increment yields the corresponding conservative lower bound. The exceptional line and its neighborhood are therefore explicit, rather than hidden behind a rank assumption.

## Appendix B. Near-Equivalent Worlds and the Testing Reduction

Substitute t(q)=iq/(1-iq) and g'=g/c. Bringing the two responses to a common denominator cancels the linear term and leaves -g(c-1)q squared. Division by |gt(q)| yields the bound in Section IV. Monotonicity of q converts scaled reactances into genuine independent real material values; membership in the material/gain domains is checked separately.

After complex whitening, the real projection along the mean difference has variance 1/2 and mean displacement D. The likelihood-ratio threshold lies halfway between the means, giving Phi(-D/sqrt(2)). If an estimator lands in the correct one of two disjoint material-tolerance sets, it can be decoded to a correct binary decision. Therefore the average estimation failure is no smaller than the optimal binary error, and the maximum of the two failures is no smaller than that average lower bound. The argument applies to arbitrary estimators, not merely least squares or locally unbiased ones.

For deterministic whitened error balls of radius beta, a mean separation at most 2 beta permits a common observation. No rule can then guarantee both disjoint material targets under every allowed discrepancy. This statement requires a genuine error-set bound, not an observed mesh difference.

## Appendix C. Reference Bound and Exact Interval Formula

The reverse triangle inequality gives ||Y_i|-a h_i|<=b_i. On Z>=a_min-b_a>0,

$$\left||Y_i|/Z-h_i\right|
\le\frac{b_i+H_i b_a}{a_{min}-b_a}.$$

Clipping to the known image interval cannot increase this error. The inverse mean-value bound supplies division by m_i. For proper complex noise, P(|n_i|>b_i)=exp(-b_i squared/sigma_i squared). Taking b_i=sigma_i sqrt(log 1000) and b_a=sigma_a Phi-inverse(0.9995), a union bound gives failure at most 0.003. Independence of the tail events is not needed for this bound.

For x=ka, m=sqrt(epsilon), and D_j(z)=(j_1(z)+z j_1'(z))/z, the exact electric reactance is

$$q(\epsilon,x)=
\frac{m j_1(mx)D_j(x)-j_1(x)D_j(mx)}
{m j_1(mx)D_y(x)-y_1(x)D_j(mx)}.$$

Use j_1(z)=sin(z)/z squared-cos(z)/z and

$$D_j(z)=\sin(z)/z+\cos(z)/z^2-\sin(z)/z^3.$$

Differentiation with m'=1/(2m) gives q'=(N'D-ND')/D squared. Finally h'=q'/(1+q squared) to the power 3/2. The implementation covers [1.5,4] and [2,5] with 2,500 and 3,000 interval boxes at 35-digit interval precision. It checks positive q, positive q', and a nonzero denominator. This is an interval-enclosure implementation, not an independent formal verification of the numerical backend or of a physical instrument. Numerical constants and source hashes are retained in `results/modal_boundary.json`.

## Appendix D. Correct Gain Profiling and Reference GLS

Set d=g-z and r_0=y-zf. Completing the square in

$$\|r_0-df\|^2/\sigma^2+|d|^2/\sigma_r^2$$

gives the rank-one GLS expression in Section V. No log determinant is added: that would correspond to a different marginalized model rather than this profile objective. A gain bound can activate, in which case the unconstrained equivalence is applied only when its solution is admissible.

For augmented weighted v=(f/sigma,1/sigma_r), b=(y/sigma,z/sigma_r), and r=b-vg with g=(v*b)/(v*v), differentiation along any real parameter gives

$$dr=-Q_v(dv)g-v\frac{(dv)^*r}{v^*v}.$$

The residual-dependent term cannot generally be deleted. Numerical tests against full-residual finite differences give relative errors at most 5.15e-9 for the tested parameters.

## Appendix E. Nuisance-Space Perturbations and Unclosed Gates

Suppose a real nuisance matrix N has full column rank, sigma_min(N)=nu, and its perturbation norm is eta_N<nu. Then

$$\|P_{N+E_N}-P_N\|\le\frac{\eta_N}{\nu-\eta_N},$$

because Q_N(N+E_N)=Q_NE_N and the perturbed pseudoinverse has norm at most 1/(nu-eta_N). Equal-rank projector geometry gives the asserted projector difference. Consequently, for material perturbation norm eta_A,

$$\sigma_{min}(Q_{N+E_N}(A+E_A))
\ge\sigma_{min}(Q_NA)-\eta_A
-\frac{\eta_N}{\nu-\eta_N}\|A\|.$$

The assumptions are essential. With N_delta=[e_1,delta e_2] and A=e_2, an arbitrarily small nonzero delta changes the nuisance rank and destroys the visible material direction. A small field error alone therefore supplies neither a derivative bound nor a nuisance-subspace bound.

For a fixed task separation delta, the relevant global quantity minimizes whitened data distance over all separated material pairs and all admitted nuisance variables. Found competitor pairs supply upper bounds; reliable positive acceptance requires valid lower bounds, coverage, and error sets. These ingredients remain unavailable for task C. Accordingly, the parent Q5 task remains incomplete.

## References

[L1] E. C. Le Ru, W. R. C. Somerville, and B. Auguie, "Radiative correction in approximate treatments of electromagnetic scattering by point and body scatterers," *Physical Review A*, 87, 012504, 2013. doi:10.1103/PhysRevA.87.012504. https://arxiv.org/abs/1210.0936v2

[L2] Z. Idriss and R. G. Raj, "Data-Driven Calibration Technique for Quantitative Radar Imaging," arXiv:2503.07316v2, 2025. https://arxiv.org/html/2503.07316v2

[L3] Y. Li, K. Lee, and Y. Bresler, "Optimal Sample Complexity for Blind Gain and Phase Calibration," arXiv:1512.07293v1. https://arxiv.org/html/1512.07293v1

[L4] D. Beutel, I. Fernandez-Corbaton, and C. Rockstuhl, "treams--A T-matrix scattering code for nanophotonic computations," arXiv:2309.03182v1. https://arxiv.org/html/2309.03182v1

[L5] V. Epp and J. G. Janz, "Spectral approach to the inverse problem for the field of arbitrary changing electric dipole," arXiv:1308.1662, 2013. https://arxiv.org/abs/1308.1662

[L6] A. de Jesus Torres, A. A. D'Amico, L. Sanguinetti, and M. Z. Win, "Cramer-Rao Bounds for Near-Field Localization," arXiv:2104.14825v2. https://arxiv.org/html/2104.14825v2

[L7] J. P. Kaipio et al., "A Bayesian approach to improving the Born approximation for inverse scattering with high contrast materials," *Inverse Problems*, 35, 084001, 2019. doi:10.1088/1361-6420/ab15f3. https://arxiv.org/html/1901.00909v2

[L8] P. Comon and L. Deruaz, "Array Self Calibration: Identifiability Issues," EUSIPCO, 1996. https://www.eurasip.org/Proceedings/Eusipco/1996/paper/ap_5.pdf

[L9] L. Bellomo, S. Pioch, M. Saillard, and K. Belkebir, "An Improved Antenna Calibration Methodology for Microwave Diffraction Tomography in Limited-Aspect Configurations," *IEEE Transactions on Antennas and Propagation*, 62(5), 2450-2462, 2014. doi:10.1109/TAP.2014.2308534.

[L10] R. Solimene, M. A. Maisto, G. Romeo, and R. Pierri, "On the Singular Spectrum of the Radiation Operator for Multiple and Extended Observation Domains," *International Journal of Antennas and Propagation*, article 585238, 2013. doi:10.1155/2013/585238.

[L11] D. Calvetti and E. Somersalo, "Spotlight, priorsketching and Bayesian approximation error paradigms," arXiv:2604.26254v1, 2026. Preprint. https://arxiv.org/html/2604.26254v1

## Evidence Availability

All new scripts and raw results accompany this draft. See `A5_REPRODUCIBILITY.md` for execution commands, `A5_FROZEN_PROTOCOL.md` for the locally registered random protocol, and `A5_LITERATURE_LEDGER.md` for access levels and unresolved prior-art questions. Original A1-A4 files were not modified. No real measurement, hardware validation, external model-provider pipeline, or automatic repository push was performed.
