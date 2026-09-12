# A2 bounded FFT VIE and Gaussian component ports

`Gaussian/code/a2/physics.py` implements the finite-grid 2-D scalar outgoing VIE
with a Toeplitz FFT internal Green action and its Euclidean discrete adjoint.
The default backend keeps the original point/off-diagonal plus equal-area-disk
self rule.  An explicitly opt-in `cell_integrated=True` backend uses 4-by-4
Gauss-Legendre source-cell quadrature and an adaptively split polar integral of
the square self cell.  It is an independent finite-grid quadrature check, not
a replacement accepted by theory or data fitting.  `ports.py` first defines every Gaussian
component's local response through the full discretized local equation, then
trains optional finite ports from incident and cross-component response
snapshots.  These ports are separate from the fixed HG patches in A1.

The N16 dense comparison tests FFT forward, adjoint, and a material directional
derivative.  `run_quadrature_ports.py` evaluates N32/N64/N128 same-physics
scattered-field comparisons, a scalar cylinder-series convention control, and
three separation/overlap cases across 3/6/12/24 response ports per component.
The block-Jacobi preconditioner uses finite-grid component ownership blocks;
Gaussian overlap and sampled boundary tails are audited rather than denied.
No result states that a finite Gaussian is compactly supported, that a chosen
port rank is sufficient generally, or that the network coordinates are
power-normalized scattering parameters.
