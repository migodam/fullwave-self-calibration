# Receiver phase sensitivity versus gain-profiled calibration information

Parent derivation, 2026-09-07. **Conditional theorem; no global originality claim.** This is a full-wave receiver-side result, not a Born approximation and not a SOM advantage theorem. Numerical verification and nearest-prior comparison are separate tasks.

## 1. Scope and expansion assumption

Fix a wavenumber k, a bounded scatterer, and world-fixed illuminations t. Observe a fixed Cartesian field component c at receiver r=R n, with |n|=1. Assume a smooth outgoing-field expansion

$$
H_{ct}(R,n)=\frac{e^{ikR}}R\left[a_{ct}(n)+\frac{b_{ct}(R,n)}R\right].
$$

On the finite sampled angular set and its neighborhoods, a and its first angular derivatives are bounded; b, its angular derivatives and R times its radial derivative are uniformly bounded as R tends to infinity. This explicit C1 remainder assumption is stronger than merely writing an O(R^-2) field expansion. The material and induced currents are fixed while differentiating receiver locations. Full multiple scattering inside the object may enter a and b without changing the argument.

Assume at least one sampled leading far-field component is nonzero, so the stacked field norm is asymptotic to a nonzero constant times R^-1. Gains may be an independent complex number per receiver component shared across illuminations. A common complex gain per receiver vector is also enough to absorb the leading propagation factor. Adding independent transmitter gains or material nuisance can only reduce the visible geometry norm further. Gains locked across different frequencies require a separate analysis: the k-dependent phase below need not belong to that more restricted nuisance space.

## 2. Exact differential decomposition

For a common translation h, define u=n·h and v=(I-n n^T)h. The receiver differential is

$$
D_h= u\,\partial_R+R^{-1}v\cdot\nabla_{\mathbb S^2}.
$$

Writing g(R)=e^{ikR}/R, differentiation and subtraction of a rowwise complex gain perturbation give

$$
D_h H-(ik-R^{-1})u H
=\frac{g(R)}R\left[
v\cdot\nabla_{\mathbb S^2}a
+u\left(\partial_R b-\frac bR\right)
+\frac1R v\cdot\nabla_{\mathbb S^2}b
\right].
$$

The subtracted term is an allowed infinitesimal gain change, the same for all illuminations at a receiver/component. Its O(k |H|) phase sensitivity is therefore **not unique geometry evidence** under the stated gain model.

Let P_G be the Euclidean orthogonal projector onto the correctly realified gain tangent. For receivers R_r=s_r R with fixed positive s_r and fixed angular directions, the bounded remainder assumptions imply

$$
\|(I-P_G)D H\|=O(R^{-2}),\qquad
\frac{\|(I-P_G)D H\|}{\|H\|}=O(R^{-1}).
$$

Proof: projection annihilates the subtracted gain vector and cannot increase norm. The angular term is O(R^-2), remaining terms O(R^-3), and the stacked field norm is Theta(R^-1). This proves an upper bound, not a nonzero lower bound: geometry may be even less observable. Uniformly conditioned whitening preserves the scaling. Under fixed per-measurement relative SNR, whitening is proportional to R, so the visible Jacobian is O(R^-1) and its Fisher information is O(R^-2). Under fixed absolute noise the corresponding upper scaling is R^-2 and R^-4. Do not confuse these noise experiments.

For a single receiver moving radially with its direction fixed, v=0 and the surviving raw derivative is O(R^-3); the relative derivative is O(R^-2). A common array translation is not radial at every receiver, so this stronger radial bound does not normally apply to the whole array.

## 3. What the result does and does not imply

1. Increasing propagation phase sensitivity does not automatically increase self-calibration information when independent complex channel gains are unknown.
2. Receiver-side angular diversity/parallax, electronic constraints, rotations and reference observations can matter more than the unprofiled k-scale phase derivative. The absolute constants depend on the full-wave far-field pattern and may worsen with k; this is a fixed-k, increasing-distance theorem, not a uniform high-frequency limit.
3. The conclusion is independent of numerical current basis and scattering strength under the explicit expansion. SOM cannot defeat the missing information. It may still provide computational benefits, which need separate tests.
4. A zero leading component, an angular radiation null, an unbounded derivative, source movement changing the induced current, range-dependent antenna orientation/pattern, or correlated whitening with unbounded condition number requires revisiting the assumptions.
5. A receiver gain shared across all its vector components still absorbs the leading scalar propagation factor, but its smaller nuisance space may preserve more polarization information in the remainder. A gain tied across frequencies does **not** in general absorb the k-dependent leading term. The existing three-frequency calibration study uses frequency-shared illumination gains, not the unrestricted single-frequency receiver-gain model here.

## 4. Algorithmic consequence to test

Use the nuisance-profiled target Jacobian, not |d phase/d position| or unprofiled Fisher information, to decide whether a proposed range/frequency acquisition is useful. A far receiver at the same raw SNR can have large apparent phase sensitivity yet weak identifiable translation. If projected information is inadequate, the admissible actions are additional angular diversity, a justified electronic reference, or a report of uncertainty; expanding a current basis cannot manufacture missing information.

The mechanism is plausibly related to longstanding near/far-field array-calibration ambiguity results. Closest-prior search has surfaced earlier array self-calibration and gain/phase identifiability work. This derivation is currently an independently obtained specialization, **not** a certified first discovery.
