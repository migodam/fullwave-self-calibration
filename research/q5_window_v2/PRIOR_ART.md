# Closest-prior audit (2026-09-11)

Access labels concern what was actually read, not whether a publisher calls a paper open access. Several PDF screenshot attempts returned internal errors; parsed full text and equations were available for the entries labelled full text. Figures are not used as evidence. No priority claim is inferred from an inaccessible source. These are the references [R1]-[R9] used in the manuscript sections.

## [R1] Gilev et al., JQSRT 131 (2013), 202-214

Title: An optimization method for solving the inverse Mie problem based on adaptive algorithm for construction of interpolating database. DOI 10.1016/j.jqsrt.2012.08.001. **Full author-hosted PDF read**, notably Sections 2-4, printed pp. 204-206.

The single nonabsorbing sphere map is assumed continuously differentiable and one-to-one on its parameter domain. Equation (3) defines a coordinatewise tolerance neighborhood B(x,e); equation (5) defines d(x)=inf_{u outside B(x,e)}||f(u)-f(x)||, its database version d1, and a neighborhood radius d2. The paper explicitly distinguishes the continuous inverse requirement from its discrete database condition. It does not provide the present shared complex gain, receiver shift, reference drift and forward-discrepancy model. Nevertheless, this is direct prior art against claiming that a task-tolerance-dependent finite forward separation, or the warning that discrete distances do not certify the continuum, is newly invented here.

Source: https://scattering.ru/papers/Gilev%20et%20al.%20-%202013%20-%20An%20optimization%20method%20for%20solving%20the%20inverse%20Mie.pdf

## [R2] Dyatlov et al., Inverse Problems 28 (2012), 045012

Title: An optimization method with precomputed starting points for solving the inverse Mie problem. DOI 10.1088/0266-5611/28/4/045012. **Full author-hosted PDF read**, especially Sections 2.2-2.3 and the covering/start construction.

Equations (14)-(15) give the electric/magnetic Mie coefficients; (21)-(22) recover them from full-angle complex amplitude functions by orthogonality. With psi_n=x j_n and xi_n=x h_n, equation (23) gives

m^2 = [(a_n xi_n-psi_n)/(a_n xi_n'-psi_n')] [(b_n xi_n'-psi_n')/(b_n xi_n-psi_n)].

Thus electric/magnetic information for material identification is explicitly old. This identity uses absolutely normalized complex coefficients, not an arbitrary unknown common gain. The precomputed-start scheme covers the parameter domain; equation (29) uses sup_{x in C(z)}||f(z)-f(x)||. Generic multistart, mode extraction and material recovery from Mie coefficients are not contributions of V2. The narrowly checked gain-invariant ratio monotonicity/noise budget differs in hypotheses, but its priority is not established by that difference alone.

Source: https://scattering.ru/papers/Dyatlov%20et%20al.%20-%202012%20-%20An%20optimization%20method%20with%20precomputed%20starting%20p.pdf

## [R3] Romanov, Maltsev and Yurkin, Optics & Laser Technology 161 (2023), 109141

Title: Retrieving refractive index of single spheres using the phase spectrum of light-scattering pattern. DOI 10.1016/j.optlastec.2023.109141. **Full author-hosted PDF read**, Sections 2-3.

Equation (1) is an angular intensity profile over 10-65 degrees. Equation (2) Fourier-transforms that real profile. Equations (4)-(7) use a refractive-index-dependent coordinate and an approximate spectral phase shift, F(v) approximately exp(i v c0/sqrt(m)) F_RGD(v/sqrt(m)). This phase is the Fourier phase of an intensity pattern, not the carrier phase of the complex scattered electric field. Its restricted working prior and synthetic uniform-noise experiment differ from the Q5 proper-complex-noise nuisance model. Phase-related spherical material inversion is therefore not new merely because the present title contains 'phase'. No claim that its approximations or noise bounds solve class C is made.

Source: https://scattering.ru/papers/Romanov%20et%20al.%20-%202023%20-%20Retrieving%20refractive%20index%20of%20single%20spheres%20usin.pdf

## [R4] Yurkin and Mishchenko, Physical Review A 97 (2018), 043824

Title: Volume integral equation for electromagnetic scattering: Rigorous derivation and analysis for a set of multilayered particles with piecewise-smooth boundaries in a passive host medium. DOI 10.1103/PhysRevA.97.043824. **Full text obtained and formulation inspected**.

This supports the importance of a correctly defined continuous Maxwell volume integral operator for piecewise material domains. It is not a source for the new numerical value 43/500, nor for the claim that a mesh difference bounds the continuum residual. The distributional local contribution and interface assumptions must survive the transition to any implemented residual verifier. V2's explicit two-ball estimate is derived in its manuscript, not attributed to this paper.

Source: https://scattering.ru/papers/Yurkin%20and%20Mishchenko%20-%202018%20-%20Volume%20integral%20equation%20for%20electromagnetic%20scatt.pdf

## [R5] Costabel, Dauge and Nedaiasl, arXiv:2302.13159v3 (2023)

Title: Stability Analysis of a Simple Discretization Method for a Class of Strongly Singular Integral Equations. **Full text read**, Sections 1.2, 1.5 and 3.5.

Equations (1.5)-(1.8) separate the distributional Maxwell VIE and its principal-value/local-term form. Lemma 3.14 and equation (3.44) bound a quasi-static resolvent by the distance to the numerical-range interval [-1/d,1-1/d]. The paper also demonstrates that the corresponding DDA discretization can have a wider numerical range than the continuum operator. Thus coercivity/resolvent reasoning is established prior work, and small or convergent-looking DDA diagnostics must not be promoted to continuum guarantees. V2 adds a finite-frequency two-ball numerical constant under its own restrictive geometry/contrast assumptions, not a new general stability principle.

Source: https://arxiv.org/pdf/2302.13159

## [R6] Miller et al., Optics Express 24 (2016), 3329-3364

Title: Fundamental limits to optical response in absorptive systems. DOI 10.1364/OE.24.003329. **Full author manuscript read**, notably equations (26)-(32).

The induced-current energy bounds depend on material factors including |chi|^2/Im chi. They constrain optical response, not automatically the distinguishability of two materials after gain/geometry profiling. Passivity, response strength and inverse information are different claims. This work is prior art for absorption-based polarization norm bounds; V2 instead uses a static/dynamic split to obtain its stated class-C constant. No shape-independent information limit is claimed.

Source: https://arxiv.org/pdf/1503.03781

## [R7] Tsitsas, Journal of Computational Mathematics 31(5) (2013), 439-448

Title: A Low-Frequency Electromagnetic Near-Field Inverse Problem for a Spherical Scatterer. DOI 10.4208/jcm.1304-m4388. **Publisher abstract/metadata only; decisive full-text formulas unavailable in this environment.** The migrated publisher site displays an inconsistent 2018 date; the journal issue is 2013.

The abstract describes interior electric-dipole excitation, a low-frequency closed approximation of the secondary field at the dipole and complex permittivity recovery. This is a direct sphere/near-field/material-identification neighbor and cannot be dismissed from an abstract. Its source-normalization, admissible-material and stability assumptions still need full-text checking before a novelty claim. It remains an explicit originality blocker.

Source: https://www.global-sci.com/jcm/article/view/12126

## [R8] Doicu, Efremenko and Yurkin, JQSRT 362 (2026), 110034

Title: An optimization tool for inverse problems of multilayered spherical particles. DOI 10.1016/j.jqsrt.2026.110034. **DLR institutional record/abstract read; linked full-PDF fetch failed.** The record was deposited in June 2026; publisher volume date is October 2026. It is relevant published/online work, not silently excluded by the training cutoff.

The record describes T-matrix recurrences/derivatives, multistart/global-local inversion, optional data transformations and covariance/posterior summaries. Without full text we do not assert which gain/reference nuisance assumptions it excludes. This entry blocks any claim that a standard spherical T-matrix plus multistart optimizer is a new method.

Source: https://elib.dlr.de/225211/

## [R9] Related direct neighbors not fully audited

Ludlow and Everitt, JOSA A 17 (2000), 2229-2235, DOI 10.1364/JOSAA.17.002229: cited within [R2] for the analytic inverse-Mie coefficient identity; original full text not read.

Takenaka and Moriyama, Optics Letters 37 (2012), 3432, DOI 10.1364/OL.37.003432: field-equivalence inversion without prior incident-field knowledge; abstract-level access only. It must be compared using its measured field/equivalent-current assumptions, not equated to an arbitrary unknown detector gain.

Tsitsas, ICEAA 2011, DOI 10.1109/ICEAA.2011.6046512: low-frequency spherical inverse-medium work; metadata/abstract only.

## Originality verdict

Not claimed new: finite inverse moduli; GLS/profile equivalence; gain cancellation by ratios; electric/magnetic Mie information; multistart; ordinary sphere solvers; operator coercivity; generic minimum-cost set cover. Newly supplied in this research record: explicit rationally covered constants, the sharp conditional amplitude-class boundary and its finite window example, a finite class-D projective range certificate, and the declared class-C residual-transfer constant. Their literature priority is unproved. Full class-C finite recovery and hardware feasibility are also unproved. These facts preclude a submission-ready TAP novelty statement.
