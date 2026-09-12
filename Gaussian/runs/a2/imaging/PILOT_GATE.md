# Independent-grid imaging pilot gate

The frozen configuration requests 30 development objects (three families ×
ten), five inverse methods, N64 inverse solves, and N128 cell-integrated data.
One preregistered `gaussian_close_00` pilot completed with the following wall
times: ordinary LM 16.38 s, fixed twofold 19.61 s, adaptive tangent 13.58 s,
random complement 19.70 s, and generalized B/A ratio selector 46.70 s.  Total
time was 115.98 s per object, projecting to about 58.0 minutes for 30 objects.

This failed the parent-specified 30-minute scheduling gate. The parent then
explicitly authorized completion of the same frozen 30-object development run;
the remaining 29 objects were launched without changing method, family,
iteration cap, seed, or data-generation configuration. The initial gate did
not pass and is retained here as the historical scheduling estimate.

The adaptive rank rule in this development configuration is a gradient-capture
heuristic. It is not a certificate: gradient inclusion in a subspace does not
in general imply inclusion of the Hessian-corrected step because of cross
coupling. A later safeguarded rule must use a residual/correction model-gap
check; no such guarantee is claimed for the recorded pilot.
