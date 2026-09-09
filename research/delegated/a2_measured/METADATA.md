# METADATA: 2001 Fresnel TM dataset `dielTM_dec8f.exp`

> INTERPRETATION CORRECTION: the documented exp(+i omega t) convention is valid;
> its outgoing spatial phase is -k*d, not +k*d. Object-center direction is an
> unknown parameter, not automatically a gauge. See the parent notice at the end.

Evidence levels: `[fulltext]` = primary document read in full; `[metadata]` =
bibliographic/registry identity only; `[measured-file]` = fact read directly
from the `.exp` file; `[empirical]` = computed from the measured values in this
pass; `[unknown]` = not found in any read source; do not guess.

## Canonical sources

- `[fulltext]` K. Belkebir, M. Saillard, "Guest Editors' Introduction: Special
  section: Testing inversion algorithms against experimental data," Inverse
  Problems **17** (2001) 1565-1571, DOI 10.1088/0266-5611/17/6/301. Author copy
  read in full from
  `https://www.fresnel.fr/perso/belkebir/Articles/Ip01Introduction_Belkebir.pdf`
  (7 pages; local copy `sources/Ip01Introduction_Belkebir.pdf` + `.txt`).
- `[metadata]` IOP supplementary-data page for the same DOI (canonical domain,
  no challenge) was the download route used by the A2 literature worker; see
  `research/delegated/a2_literature/public_data_audit.md`.
- `[metadata]` Weir 1974 waveguide permittivity technique is cited by the
  intro for the dielectric measurement; not re-read here.

## File layout (intro Section 6; verified against the file)

Seven whitespace columns after a 10-line `#` header:

| col | meaning | details |
|---|---|---|
| 1 | source/view index | 1..36; angle = (col1-1)*10 deg; source at distance `de` from setup centre |
| 2 | receiver index | 1..72; angle = (col2-1)*5 deg; receiver at distance `dr` from centre |
| 3 | operating frequency | GHz (1..8 here); **not** a sequential index |
| 4 | Re(total electric field) | object present |
| 5 | Im(total electric field) | object present |
| 6 | Re(incident electric field) | measured on the circle **without** the target |
| 7 | Im(incident electric field) | same as col 6 |

So scattered = (col4 + i col5) - (col6 + i col7). All columns are electric
field; TM means E parallel to the cylinder axis. Time dependence `exp(+i w t)`
stated in the intro and confirmed `[empirical]`: the measured incident phase
advances as +k*d from the emitter (coherence 0.92 vs 0.21 for the opposite
sign, view 1 @ 1 GHz).

## Geometry `[fulltext]`

- emitter-centre distance `de = 720 mm +/- 3 mm`
- receiver-centre distance `dr = 760 mm +/- 3 mm`
- target rotation 0..350 deg in 10 deg steps (36 views)
- receiver sweep 60..300 deg in 5 deg steps (49 receivers per view)
- dielectric targets: one or two filled circular cylinders, radius `a = 15 mm`
- cylinder centre "about 30 mm from the azimuthal positioner axis"
- relative permittivity real part `eps_r = 3 +/- 0.3` (waveguide technique)

`[measured-file]` receiver windows: view 1 = indices 13..61 (60..300 deg),
each later view shifts +2 indices (wrap mod 72), consistent with the 60..300
deg sweep re-anchored per view.

## Explicit unknowns

- `[unknown]` direction of the 30 mm offset (only the magnitude is published;
  the direction is not stated). The +x anchor used for the nominal fit is an
  arbitrary gauge choice.
- `[unknown]` exact offset magnitude (only "about 30 mm").
- `[unknown]` loss tangent / imaginary permittivity of the dielectric.
- `[unknown]` material name of the filled cylinder.
- `[unknown]` field units and global calibration scale (the intro does not
  state units; values are ~0.1-1).
- `[unknown]` emitter horn aperture/phase-centre details beyond `de` and the
  per-view incident measurement; the intro states the incident field in the
  target area "has to be estimated, or optimized along the reconstruction
  process".
- `[unknown]` receiver antenna factors as a function of angle.

## Data integrity `[measured-file]`

- SHA256 `476cc9d1cfc98797545ab4adf69302dc5aeb45848a24cf8d7b3d222940cc79eb`
  (matches the literature worker's `SHA256SUMS.txt`)
- 14112 rows = 36 views x 49 receivers x 8 frequencies; 7 columns; 10 header
  lines; all finite; zero duplicate (view, receiver, freq) keys; full
  8-frequency coverage in every recorded cell.
- `[empirical]` incident field is rotationally consistent across views
  (median relative deviation < 1%), i.e. view v pattern equals view 1 rotated
  by the source step.

## Unit semantics

The fits use scattered = total - incident with a per-frequency complex gain,
so results are scale-free; no physical unit conversion is attempted.
# Parent interpretation warning, 2026-09-06

The primary descriptor's exp(+i omega t) statement is valid (p. 1570), but below
"phase advances as +k*d" has the outgoing sign reversed. Under this time
convention outgoing waves have exp(-i*k*d). An unpublished object-center
direction is not automatically a gauge with fixed sensor coordinates. See the
parent `a2_research/FINAL_REVIEW.md`; original text below is retained for audit.
The full parent warning above refers to the original metadata text, retained
unchanged for traceability.
