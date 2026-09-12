# Fixed-geometry Fresnel measured-cylinder diagnostic

`Gaussian/code/a2/measured.py` uses the local scientific-use Fresnel
`dielTM_dec8f.exp` copy under the established conjugation and back-arc
incident-source-factor convention.  The acquisition geometry is declared and
fixed: 2/4/6 GHz, 12 uniformly spaced published source views, all available
receiver records, and global even receiver channels for training with odd
channels held out.  No geometry calibration is performed.

The primary fits profile a complex gain independently per frequency from the
scattered training data.  It is an additional nuisance parameter, so the run
is a gain-profiled model-mismatch diagnostic and retains material-amplitude
ambiguity.  `fixed_gain_refits.py` starts from selected states and refits with
gain fixed to one.  The softened known-shape disk is an oracle diagnostic, not
a fair general inverse baseline.  The disk image discrepancy is a proxy only,
not an approximation-floor theorem.  This one 2-D cylinder does not justify
an NN conclusion, a diverse-target claim, a calibrated material estimate, or
3-D validation.
