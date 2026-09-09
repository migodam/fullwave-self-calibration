# Consolidated round-3 results (geometry-lifted tri-space SOM)
## How to read this
All numbers below are read directly from the result JSONs listed in each table caption; medians and 25th/75th-percentile IQRs are computed over the recorded noise seeds in those files. No physics was recomputed. Fields that are absent in a file are shown as `—`. Scientific notation is used only for values < 1e-3 or > 1e3.

### Table R3-1: E10 setting sweep (source: `results_e10_settings_sweep.json`)

| setting | M | k | aperture | rank G_s | sv_max | sv_min | lift r4 resid | lift r6 resid | Q^H R r4 | Q^H R r6 | subspace r4 n_vis (hid) | subspace r6 n_vis (hid) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M8_k12_full | 8 | 12 | full | 8 | 0.00329 | 0.0025 | 0.823 | 0.557 | 1.817e-18 | 4.364e-18 | 3 (hid 0) | 1 (hid 2) |
| M12_k12_full | 12 | 12 | full | 12 | 0.00308 | 0.00224 | 0.474 | 0.452 | 8.125e-18 | 1.170e-17 | 3 (hid 0) | 3 (hid 0) |
| M16_k12_full | 16 | 12 | full | 16 | 0.00365 | 0.00108 | 0.517 | 0.4 | 4.156e-18 | 6.018e-18 | 3 (hid 0) | 3 (hid 0) |
| M12_k8_full | 12 | 8 | full | 12 | 0.00472 | 6.859e-04 | 0.491 | 0.354 | 3.287e-18 | 8.543e-18 | 3 (hid 0) | 3 (hid 0) |
| M12_k16_full | 12 | 16 | full | 12 | 0.00288 | 0.00222 | 0.854 | 0.519 | 2.711e-18 | 9.925e-18 | 3 (hid 0) | 3 (hid 0) |
| M12_k12_half | 12 | 12 | half | 12 | 0.0042 | 1.578e-04 | 0.487 | 0.357 | 5.851e-18 | 8.339e-18 | 3 (hid 0) | 3 (hid 0) |
| M16_k12_threequarter | 16 | 12 | threequarter | 16 | 0.00402 | 2.676e-04 | 0.864 | 0.686 | 5.652e-18 | 4.047e-18 | 3 (hid 0) | 3 (hid 0) |

Q^H R columns are `orthogonality_QH_R_norm_over_R_norm` (|Q^H R|/|R|). n_vis (hid) is written `n_vis(hidden_rank)` from each setting's subspace sweep.

**Table R3-1b: Pose error median over seeds** (source: `results_e10_settings_sweep.json`)

| setting | direct | reduced_r4 | reduced_r6 | wrongpose | known_pose | known_alpha |
|---|---|---|---|---|---|---|
| M8_k12_full | 0.0158 | 0.0158 | 0.0442 | 0.112 | 0 | 0.0106 |
| M12_k12_full | 0.00761 | 0.00761 | 0.00761 | 0.112 | 0 | 0.00897 |
| M16_k12_full | 0.0148 | 0.0148 | 0.0148 | 0.112 | 0 | 0.00803 |
| M12_k8_full | 0.00177 | 0.00177 | 0.00177 | 0.112 | 0 | 0.00165 |
| M12_k16_full | 0.00815 | 0.00815 | 0.00815 | 0.112 | 0 | 0.01 |
| M12_k12_half | 0.0104 | 0.0102 | 0.0271 | 0.112 | 0 | 0.86 |
| M16_k12_threequarter | 0.0102 | 0.0102 | 0.0102 | 0.112 | 0 | 0.00837 |

**Table R3-1c: Map error median over seeds** (source: `results_e10_settings_sweep.json`)

| setting | direct | reduced_r4 | reduced_r6 | wrongpose | known_pose | known_alpha |
|---|---|---|---|---|---|---|
| M8_k12_full | 0.0124 | 0.0124 | 0.0274 | 0.209 | 0.0119 | 0 |
| M12_k12_full | 0.0121 | 0.0121 | 0.0121 | 0.26 | 0.0114 | 0 |
| M16_k12_full | 0.015 | 0.015 | 0.015 | 0.24 | 0.0116 | 0 |
| M12_k8_full | 0.0118 | 0.0118 | 0.0118 | 0.242 | 0.0117 | 0 |
| M12_k16_full | 0.00707 | 0.00707 | 0.00707 | 0.245 | 0.00922 | 0 |
| M12_k12_half | 0.0284 | 0.0285 | 0.0312 | 1.33 | 0.00812 | 0 |
| M16_k12_threequarter | 0.00644 | 0.00644 | 0.00644 | 0.254 | 0.00641 | 0 |

### Table R3-2: E10b M12/M16 high-rank transition (source: `results_e10b_m12_m16_highrank.json`)

| M | k | r_transition | hidden_rank | n_vis | has_hidden_rank_3 | empty-lift trailing r |
|---|---|---|---|---|---|---|
| 12 | 12 | 10 | 2 | 1 | False | 11,12 |
| 16 | 12 | 14 | 2 | 1 | False | 15,16 |

**Table R3-2b: E10b median over seeds (direct and reduced at the selected transition/empty-lift ranks)** (source: `results_e10b_m12_m16_highrank.json`)

| M | method | r | regime | pose med | visible med | hidden med | map med |
|---|---|---|---|---|---|---|---|
| 12 | direct | — | direct reference | 0.00761 | 0.00221 | 0.00725 | 0.0121 |
| 12 | reduced_r10 | 10 | hidden | 0.0759 | 0.0239 | 0.072 | 0.0463 |
| 12 | reduced_r11 | 11 | empty_lift_diagnostic | 0.00761 | 0.00761 | 1.062e-18 | 0.0121 |
| 16 | direct | — | direct reference | 0.0148 | 0.00339 | 0.0144 | 0.015 |
| 16 | reduced_r14 | 14 | hidden | 0.111 | 0.0642 | 0.0905 | 0.0895 |
| 16 | reduced_r15 | 15 | empty_lift_diagnostic | 0.0148 | 0.0148 | 1.227e-18 | 0.015 |

### Table R3-3: E6b multi-transmitter robustness (source: `results_e6b_multitx_robust.json`)

Main groups requested for the paper: monopole L2 tx_only direct and wrongpose, monopole L3 tx_only direct, monopole L2 co_moving direct, dipole L2 tx_only direct. med = 50th percentile over seeds; IQR = 25th-75th percentile interval.

**Table R3-3a: medians over seeds**

| group | SNR dB | success | pose med | theta med | tx med | ty med | map med | final res med |
|---|---|---|---|---|---|---|---|---|
| dipole L2 tx_only direct_tx | 20 | 8/8 | 0.0084 | 0.008 | 0.00178 | 0.00101 | 0.0216 | 0.095 |
| dipole L2 tx_only direct_tx | 30 | 8/8 | 0.00266 | 0.0025 | 5.668e-04 | 3.194e-04 | 0.00683 | 0.0301 |
| dipole L2 tx_only direct_tx | 40 | 8/8 | 8.412e-04 | 7.921e-04 | 1.798e-04 | 1.014e-04 | 0.00216 | 0.00952 |
| monopole L2 tx_only direct_tx | 15 | 8/8 | 0.0152 | 0.0145 | 0.00318 | 0.00194 | 0.0421 | 0.168 |
| monopole L2 tx_only direct_tx | 20 | 8/8 | 0.00854 | 0.00808 | 0.00181 | 0.00106 | 0.0234 | 0.0952 |
| monopole L2 tx_only direct_tx | 25 | 8/8 | 0.00481 | 0.00455 | 0.00102 | 5.898e-04 | 0.0131 | 0.0536 |
| monopole L2 tx_only direct_tx | 30 | 8/8 | 0.00271 | 0.00256 | 5.764e-04 | 3.291e-04 | 0.00735 | 0.0301 |
| monopole L2 tx_only direct_tx | 35 | 8/8 | 0.00152 | 0.00144 | 3.250e-04 | 1.843e-04 | 0.00413 | 0.0169 |
| monopole L2 tx_only direct_tx | 40 | 8/8 | 8.566e-04 | 8.103e-04 | 1.830e-04 | 1.034e-04 | 0.00232 | 0.00952 |
| monopole L2 tx_only wrongpose | 15 | 8/8 | 0.112 | 0.05 | 0.08 | 0.06 | 0.299 | 0.714 |
| monopole L2 tx_only wrongpose | 20 | 8/8 | 0.112 | 0.05 | 0.08 | 0.06 | 0.292 | 0.708 |
| monopole L2 tx_only wrongpose | 25 | 8/8 | 0.112 | 0.05 | 0.08 | 0.06 | 0.289 | 0.706 |
| monopole L2 tx_only wrongpose | 30 | 8/8 | 0.112 | 0.05 | 0.08 | 0.06 | 0.288 | 0.705 |
| monopole L2 tx_only wrongpose | 35 | 8/8 | 0.112 | 0.05 | 0.08 | 0.06 | 0.288 | 0.705 |
| monopole L2 tx_only wrongpose | 40 | 8/8 | 0.112 | 0.05 | 0.08 | 0.06 | 0.288 | 0.705 |
| monopole L2 co_moving direct_comoving | 20 | 8/8 | 0.0222 | 0.0195 | 0.00261 | 0.00589 | 0.034 | 0.0897 |
| monopole L2 co_moving direct_comoving | 30 | 8/8 | 0.00692 | 0.00641 | 8.125e-04 | 0.00176 | 0.0106 | 0.0286 |
| monopole L2 co_moving direct_comoving | 40 | 8/8 | 0.00223 | 0.00207 | 2.560e-04 | 5.483e-04 | 0.00332 | 0.00904 |
| monopole L3 tx_only direct_tx | 20 | 8/8 | 0.00422 | 0.00368 | 0.00181 | 0.00107 | 0.0265 | 0.0945 |
| monopole L3 tx_only direct_tx | 30 | 8/8 | 0.00133 | 0.00117 | 5.451e-04 | 3.467e-04 | 0.00837 | 0.03 |
| monopole L3 tx_only direct_tx | 40 | 8/8 | 4.192e-04 | 3.711e-04 | 1.697e-04 | 1.104e-04 | 0.00265 | 0.0095 |

**Table R3-3b: IQR (q25-q75) over seeds**

| group | SNR dB | pose IQR | theta IQR | tx IQR | ty IQR | map IQR | final res IQR |
|---|---|---|---|---|---|---|---|
| dipole L2 tx_only direct_tx | 20 | 0.00524-0.0113 | 0.00471-0.0109 | 9.746e-04-0.00241 | 7.902e-04-0.00138 | 0.0195-0.0286 | 0.0933-0.096 |
| dipole L2 tx_only direct_tx | 30 | 0.00162-0.00364 | 0.00149-0.00355 | 3.048e-04-7.584e-04 | 2.389e-04-4.257e-04 | 0.00599-0.00911 | 0.0297-0.0305 |
| dipole L2 tx_only direct_tx | 40 | 5.085e-04-0.00115 | 4.666e-04-0.00112 | 9.608e-05-2.395e-04 | 7.447e-05-1.335e-04 | 0.00188-0.00289 | 0.00941-0.00966 |
| monopole L2 tx_only direct_tx | 15 | 0.00935-0.0195 | 0.00827-0.0189 | 0.00238-0.00462 | 0.0013-0.00265 | 0.0352-0.0517 | 0.162-0.169 |
| monopole L2 tx_only direct_tx | 20 | 0.00521-0.0112 | 0.00477-0.0109 | 0.00132-0.00259 | 7.101e-04-0.00144 | 0.0199-0.029 | 0.0928-0.0962 |
| monopole L2 tx_only direct_tx | 25 | 0.0029-0.00633 | 0.00267-0.00619 | 7.353e-04-0.00146 | 3.922e-04-7.975e-04 | 0.0112-0.0163 | 0.0526-0.0541 |
| monopole L2 tx_only direct_tx | 30 | 0.00162-0.00357 | 0.00151-0.00349 | 4.117e-04-8.188e-04 | 2.183e-04-4.431e-04 | 0.00633-0.00916 | 0.0297-0.0305 |
| monopole L2 tx_only direct_tx | 35 | 9.084e-04-0.00201 | 8.440e-04-0.00197 | 2.309e-04-4.601e-04 | 1.221e-04-2.477e-04 | 0.00356-0.00515 | 0.0167-0.0172 |
| monopole L2 tx_only direct_tx | 40 | 5.096e-04-0.00113 | 4.740e-04-0.00111 | 1.296e-04-2.586e-04 | 6.846e-05-1.388e-04 | 0.002-0.0029 | 0.00941-0.00965 |
| monopole L2 tx_only wrongpose | 15 | 0.112-0.112 | 0.05-0.05 | 0.08-0.08 | 0.06-0.06 | 0.288-0.308 | 0.709-0.726 |
| monopole L2 tx_only wrongpose | 20 | 0.112-0.112 | 0.05-0.05 | 0.08-0.08 | 0.06-0.06 | 0.287-0.298 | 0.705-0.715 |
| monopole L2 tx_only wrongpose | 25 | 0.112-0.112 | 0.05-0.05 | 0.08-0.08 | 0.06-0.06 | 0.287-0.293 | 0.704-0.71 |
| monopole L2 tx_only wrongpose | 30 | 0.112-0.112 | 0.05-0.05 | 0.08-0.08 | 0.06-0.06 | 0.287-0.29 | 0.704-0.707 |
| monopole L2 tx_only wrongpose | 35 | 0.112-0.112 | 0.05-0.05 | 0.08-0.08 | 0.06-0.06 | 0.287-0.289 | 0.704-0.706 |
| monopole L2 tx_only wrongpose | 40 | 0.112-0.112 | 0.05-0.05 | 0.08-0.08 | 0.06-0.06 | 0.287-0.288 | 0.704-0.705 |
| monopole L2 co_moving direct_comoving | 20 | 0.0204-0.0306 | 0.0182-0.0265 | 0.00188-0.00415 | 0.00319-0.00916 | 0.0239-0.0516 | 0.0853-0.092 |
| monopole L2 co_moving direct_comoving | 30 | 0.00604-0.00991 | 0.00523-0.00875 | 5.915e-04-0.00132 | 0.00103-0.00283 | 0.00744-0.0163 | 0.0268-0.0293 |
| monopole L2 co_moving direct_comoving | 40 | 0.00185-0.00314 | 0.00159-0.00278 | 1.868e-04-4.182e-04 | 3.275e-04-8.873e-04 | 0.00234-0.00515 | 0.00846-0.00928 |
| monopole L3 tx_only direct_tx | 20 | 0.00376-0.00617 | 0.00241-0.00551 | 0.00102-0.00244 | 8.217e-04-0.00153 | 0.024-0.0305 | 0.0936-0.0955 |
| monopole L3 tx_only direct_tx | 30 | 0.00118-0.00197 | 7.722e-04-0.00177 | 3.419e-04-7.732e-04 | 2.604e-04-4.761e-04 | 0.00763-0.00959 | 0.0299-0.0303 |
| monopole L3 tx_only direct_tx | 40 | 3.728e-04-6.259e-04 | 2.454e-04-5.618e-04 | 1.099e-04-2.445e-04 | 8.241e-05-1.497e-04 | 0.00242-0.00302 | 0.00945-0.00958 |

**Table R3-3c: dipole FD checks at eps=1e-4 (max column-wise relative error)**

| pose_mode | matrix | convention | eps | max_rel_err | per-column |
|---|---|---|---|---|---|
| tx_only | B_t | receiver_fixed | 1.000e-04 | 1.740e-07 | 1.276e-07, 1.740e-07, 2.404e-08 |
| co_moving | B_t | receiver_fixed | 1.000e-04 | 1.763e-07 | 1.254e-07, 1.763e-07, 2.536e-08 |
| co_moving | B_total | co_moving | 1.000e-04 | 2.169e-07 | 2.169e-07, 1.836e-07, 2.039e-08 |

**Table R3-3d: dipole B_t gauge-rank summary (from `gauge_partA_dipole`)**

| model | L | complex_rank | sv max | sv min | realified colnormed sv | theta resid rel |
|---|---|---|---|---|---|---|
| dipole | 1 | 3 | 0.00555 | 3.308e-06 | 1.41, 1.01, 0.0583 | 0.00798 |
| dipole | 2 | 3 | 0.00667 | 0.00123 | 1.05, 0.996, 0.948 | 0.992 |

### Table R3-4: E11 VP state-null test (source: `results_e11_vp_state_null.json`)


**Table R3-4a: median over seeds {0,1}, zero p_init (r=6 retained basis)**

| method | success | pose med | visible med | hidden med | map med | data res med | discarded cur med |
|---|---|---|---|---|---|---|---|
| direct | 2/2 | 0.0103 | 0.00249 | 0.00941 | 0.0178 | 0.025 | — |
| reduced_r6 | 2/2 | 0.0442 | 0.00258 | 0.0441 | 0.0346 | 0.0402 | — |
| vp_hard_r6 | 0/2 | 0.112 | 0.103 | 0.0441 | 0.336 | 0.438 | 0.861 |
| vp_soft_r6 | 2/2 | 0.0166 | 0.00316 | 0.016 | 0.838 | 0.0623 | 0.869 |

Note: the file also stores a diagnostic rerun vp_hard_r6_pert0 (init_label `pert_pos`, seed 0). It is kept in `runs_all` of the summary JSON but excluded from the seed medians above.

**Table R3-4b: linearized VP pose signature (B_red singular values after removing alpha columns, at alpha_init/p_init)**

| signature | visible pose dof above cut | B_red singular values | pose col norm 1 | pose col norm 2 | pose col norm 3 |
|---|---|---|---|---|---|
| hard | 3 | 8.74, 0.00258, 0.00151 | 0.00293 | 0.00275 | 8.89 |
| soft | 3 | 5.138e-04, 2.238e-04, 1.091e-04 | 4.963e-04 | 3.076e-04 | 1.265e-04 |
| full_physics_alpha_only | 3 | 9.514e-04, 4.533e-04, 2.110e-04 | 9.043e-04 | 6.316e-04 | 2.482e-04 |

Subspace anchor (r6 at alpha_init/p_init): n_vis = 1, hidden_rank = 2, hid_sv = 5.679e-04, 5.372e-19, 1.710e-19.

Interpretation (first stored sentence): Central question: can variable projection (hard state equality cur=P(p) J_phys, c eliminated) recover the two hidden pose directions that reduced_r6 freezes?

### Table R3-5: prior anchor points (sources: `results_e5_final.json`, `results_e9_grid_scenes.json`, `results_e6.json`, `results_e7b.json (listed as results_e7b_soft_state.json)`)


**Table R3-5a: E5 final, medians over seeds by p_init and method (lossless reduced_r4 vs direct; r6 hidden freeze)**

| p_init | method | r | pose med | visible med | hidden med | map med | success |
|---|---|---|---|---|---|---|---|
| pert_neg | direct | — | 0.00674 | 0.00198 | 0.00594 | 0.0148 | 2/2 |
| pert_neg | reduced_r4 | 4 | 0.00674 | 0.00674 | 3.086e-18 | 0.0148 | 2/2 |
| pert_neg | reduced_r6 | 6 | 0.0442 | 0.00271 | 0.0441 | 0.032 | 2/2 |
| pert_neg | wrongpose | — | — | — | — | — | 0/0 |
| pert_pos | direct | — | 0.00667 | 0.00195 | 0.00624 | 0.0148 | 2/2 |
| pert_pos | reduced_r4 | 4 | 0.00667 | 0.00667 | 2.059e-18 | 0.0148 | 2/2 |
| pert_pos | reduced_r6 | 6 | 0.15 | 0.0293 | 0.147 | 0.152 | 2/2 |
| pert_pos | wrongpose | — | — | — | — | — | 0/0 |
| zero | direct | — | 0.00674 | 0.00198 | 0.00594 | 0.0148 | 2/2 |
| zero | reduced_r4 | 4 | 0.00674 | 0.00674 | 3.086e-18 | 0.0148 | 2/2 |
| zero | reduced_r6 | 6 | 0.0442 | 0.00271 | 0.0441 | 0.032 | 2/2 |
| zero | wrongpose | — | 0.112 | 0.103 | 0.0441 | 0.204 | 2/2 |

**Table R3-5b: E9 grid/scene medians by scene and method (stored `summary_groups`)**

| scene | N | method | pose med | map med | final res med | success |
|---|---|---|---|---|---|---|
| alt | 16 | direct | 0.0121 | 0.0315 | 0.0256 | 3/3 |
| alt | 16 | reduced_r4 | 0.0124 | 0.0293 | 0.0256 | 3/3 |
| alt | 16 | reduced_r6 | 0.0445 | 0.036 | 0.0467 | 3/3 |
| alt | 16 | wrongpose | 0.112 | 0.617 | 0.259 | 3/3 |
| base | 24 | direct | 0.0158 | 0.0124 | 0.0278 | 3/3 |
| base | 24 | reduced_r4 | 0.0158 | 0.0124 | 0.0278 | 3/3 |
| base | 24 | reduced_r6 | 0.0463 | 0.0322 | 0.046 | 3/3 |
| base | 24 | wrongpose | 0.112 | 0.209 | 0.282 | 3/3 |
| base | 32 | direct | 0.0158 | 0.0124 | 0.0278 | 3/3 |
| base | 32 | reduced_r4 | 0.0158 | 0.0124 | 0.0278 | 3/3 |
| base | 32 | reduced_r6 | 0.0432 | 0.0265 | 0.0446 | 3/3 |
| base | 32 | wrongpose | 0.112 | 0.209 | 0.282 | 3/3 |

**Table R3-5c: E6 monopole B_t gauge-rank rows (`gauge_partA`, monopole model)**

| model | L | complex_rank | smallest complex sv | theta resid rel |
|---|---|---|---|---|
| monopole | 1 | 2 | 7.402e-19 | 2.534e-15 |
| monopole | 2 | 3 | 8.877e-04 | 0.989 |
| monopole | 3 | 3 | 0.00121 | 0.979 |
FD pass at eps=1e-4 (stored `fd.pass_at_eps_1e-4_max_rel`): 1.930e-07.

**Table R3-5d: E7b soft-state null - stored seed medians (direct/reduced_r6 baselines + state_c_soft lambda sweep)**

| method | lambda | pose med | visible med | hidden med | map med | full data res med | success |
|---|---|---|---|---|---|---|---|
| direct | — | 0.0103 | 0.00249 | 0.00941 | 0.0178 | 0.025 | 2/2 |
| reduced_r6 | — | 0.0442 | 0.00258 | 0.0441 | 0.0346 | 0.0402 | 2/2 |
| state_c_soft | 0 | 0.0847 | 0.0802 | 0.0273 | 0.445 | 0.462 | 2/2 |
| state_c_soft | 1.000e-04 | 0.464 | 0.0895 | 0.455 | 0.122 | 0.35 | 0/2 |
| state_c_soft | 0.001 | 0.47 | 0.0982 | 0.459 | 0.2 | 0.375 | 2/2 |
| state_c_soft | 0.01 | 0.814 | 0.141 | 0.802 | 1 | 1 | 0/2 |
| state_c_soft | 0.03 | 0.828 | 0.143 | 0.815 | 1 | 1 | 0/2 |
| state_c_soft | 0.1 | 0.831 | 0.144 | 0.819 | 1 | 1 | 0/2 |
| state_c_soft | 0.3 | 0.839 | 0.146 | 0.826 | 1 | 1 | 0/2 |
| state_c_soft | 1 | 0.82 | 0.138 | 0.808 | 1 | 1 | 0/2 |
| state_c_soft | 3 | 0.82 | 0.138 | 0.808 | 1 | 1 | 0/2 |
Every state_c_soft row has hidden-error median >= 0.0273, larger than the direct baseline's 0.00941; consistent with the stored E7b goal/notes that state consistency does not recover the hidden pose directions that reduced_r6 freezes.
