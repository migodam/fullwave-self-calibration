# A4 closure addendum: geometry information is not material recovery

10 September 2026. Read alongside the unchanged A4 working manuscript. This
supplement does not upgrade it to submission readiness.

## Interpretation to retain in the main text

The exact radial/tangential spectrum concerns receiver geometry after profiling
radial material and a complex frame coefficient. Under independently free
frequency gains, the same acquisition cannot identify the material coefficient.
Consequently the fixed-relative-SNR optimum kR=1 is not a material-reconstruction
window. A finite-SNR material design must profile geometry and electronics rather
than reuse the geometry spectrum with a different interpretation.

For nonzero measured blocks y_b=g_b f_b(alpha), unrestricted independent complex
gains make two materials indistinguishable precisely when their block responses
are complex-collinear. A radial single-mode response is an exact full-wave
instance. Bounded-gain counterexamples additionally require admissible compensating
gains. Shared-frequency gains impose a common, rather than blockwise, multiplier.

For a differentiable weak-contrast expansion f(s)=s v1+s² v2+O(s³), with v1 nonzero,
the shared-gain-profiled log-scale derivative is s² Q_v1 v2+O(s³). Thus a measured
second-order component transverse to the first-order response is sufficient for
a local visible direction. A nonzero internal second-order term alone is not.
Neither statement establishes novelty or global stable recovery.

## New causal control

Six deterministic contrast scales (.01,.03,.1,.3,1,2) were tested at k=18 and
receiver radius .6 m in the previous two-sphere, four-illumination configuration.
The complete lmax=4 Treams cluster was compared with the coherent sum of exact
isolated-sphere fields, retaining within-sphere response but removing rescattering
between spheres. Two central log-derivative steps were checked. A physical
first-Born dyadic volume integral supplied a homogeneous negative control.

At every tested scale, the isolated-sphere sum retained a slightly larger
shared-gain-profiled scale sensitivity than the interacting cluster. At scales
.01 and2 the respective relative sensitivities were .001260 versus .001243,
and .072064 versus .071205. Therefore inter-sphere rescattering is not necessary
for this observed visibility and did not improve this metric in these cases.
No universal statement about multiple-scattering benefit follows. Internal
response and the relative responses of distinct spheres were not further isolated.

The Born scale direction vanished after gain profiling to1.40e-16. Its two-grid
quadrature difference was .188%; this is not a continuum error bound. The maximum
relative step change of the profiled full-wave derivatives was4.92e-7. Complete
supplement runtime was4.99s, excluding previous experiments and development.
These are mechanism diagnostics, not noisy inversion outcomes or multi-material
recovery. The original complete-cluster lmax3/4 check does not itself establish
isolated-sphere truncation convergence.

## Manuscript consequence

Retain the restricted modal geometry theorem and its exact limits. Withdraw any
interpretation that it proves a material window or that the new cluster examples
establish a benefit caused by inter-sphere rescattering. Keep the failed risk
selector and conditional coverage tests as supporting negative evidence. The next
substantive gate is independently unknown material regions, physically justified
gain sharing/reference, finite separation above an error floor, and actual joint
recovery. Priority against electromagnetic dipole localization remains open.

Protocol: ROUND2_PROTOCOL.md. Raw results: results/interaction_ablation.json.
Source and protocol hashes are stored with the results. Original A4 and round-one
files are unchanged. No fresh native Agentic-AI-Scientist model loop was executed.
