# INVALIDATED: Family-1 harness snapshot (pre-self-cell-correction)

Status: **INVALIDATED** for scientific use. Retained verbatim as evidence; do
not cite any numbers, figures, or source in this folder as an admitted Family-1
result.

Date archived (UTC): 2026-09-03T09:23:39Z
Experiment root:
`/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry`

## Why this snapshot is invalidated

A parent mathematical audit (`context/PARENT_CORRECTIONS.md`, and workshop
rule 16 in `context/workshop.md`) found that the implemented self-cell rule in
this snapshot is incomplete. For the declared Green function

    g(r, r') = (i/4) H_0^{(1)}(k_b |r - r'|)

with equal-area disk radius a = h/sqrt(pi), the correct disk integral is

    I_self = (i*pi*a/(2*k_b)) * H_1^{(1)}(k_b*a) - 1/k_b^2,

because lim_{r->0} r H_1^{(1)}(k_b r)/k_b = -2i/(pi*k_b^2). This snapshot's
`self_cell_green` omitted the lower-endpoint term `-1/k_b^2`, so its alleged
cell integral tends to `1/k_b^2` instead of zero as h -> 0 and the diagonal
`[G_D]_nn` does not tend to the corrected value.

Per rule 16 and PARENT_CORRECTIONS.md, a pilot produced with the incomplete
self-cell rule is not admissible evidence and must be rerun. The source-position
gradient sign in this snapshot (`grad_s g = -grad_z g`) had already been
corrected and was retained in the follow-up run; that part is not invalidated.

## Contents (byte-for-byte copies at archive time)

```text
src/helmholtz.py                       Helmholtz harness with incomplete self-cell rule
src/family1_pilot.py                   Family-1 pilot driver producing the invalidated JSON/figure
results/family1_pilot_results.json     Family-1 pilot results (invalidated)
figures/family1_fd_convergence.png     FD-convergence figure (invalidated)
logs/environment.md                    Prior environment/run log (describes invalidated run)
notes/family1_pilot_head_dump.txt      Prior run transcript
notes/helmholtz_tail_dump.txt          Prior source dump
notes/helmholtz_AB_tail_dump.txt       Prior source dump
```

## SHA-256 checksums

```text
00a1204f60ec52f9b89c08554eade1282547be435b59d4f36c26c2aa1024ec2a  src/helmholtz.py
d8def2e4ca1eb7fbfa6bc051f91ab806c77002067a952430a3b117dbfda046bd  src/family1_pilot.py
eacb59cbe24fc6f61c209f741bbb2b6c743bffd52bc0a2a9005d5ef18561dcdb  results/family1_pilot_results.json
05de02e68c6f24bac2af3c654a53b1fa355a2232a3687bbb9549a557e946ad8f  figures/family1_fd_convergence.png
abebf2afbe662ee39d7eee8933dd3a6fe19e9d8f5fb7c8999dc0fd0201b73314  logs/environment.md
9b096fca531565b77fc08e64d3ecca490d4e3ae664e750996f8e06e57910e316  notes/family1_pilot_head_dump.txt
3060b4399e0d8fafca21a25983f555bc68c92d649afbafb6f1b38f6143e0c731  notes/helmholtz_tail_dump.txt
7667de4558b5e774868b045db63665ce16222651a73e6336f954c1b87db13a8a  notes/helmholtz_AB_tail_dump.txt
```

Checksums were verified against the original live files immediately before
copying and again on the copies after copying. A machine-readable copy of the
same digests is in `checksums.txt`.

## Replacement evidence

The corrected Family-1 rerun and new self-cell validations live in the normal
live locations:

- `src/helmholtz.py` (corrected `self_cell_green` and documentation)
- `src/test_family1_corrections.py` (unit tests, including a sign test that
  fails for the old gradient sign)
- `src/validate_self_cell.py` and `results/self_cell_quadrature_validation.json`
- `src/family1_grid_refinement.py` and `results/family1_grid_refinement.json`
- `results/family1_pilot_results.json` (corrected rerun)
- `notes/family1_self_cell_correction_report.md`
