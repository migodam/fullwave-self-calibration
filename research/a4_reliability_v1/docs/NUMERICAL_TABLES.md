# Executed A4 numerical summary

## Nine matched baselines

| Method | Median geometry (mm) | Median material (%) | Mean scaled task loss | Median charged wall time (s) |
|---|---:|---:|---:|---:|
| coarse_to_fine | 0.6355 | 0.856 | 0.050674 | 0.061069 |
| conditional_eem | 0.5562 | 0.9181 | 0.18598 | 0.040222 |
| full_fine | 0.6355 | 0.856 | 0.050674 | 0.047384 |
| isotropic | 2.447 | 7.397 | 6.8081 | 0.039566 |
| low_only | 1.978 | 6.29 | 5.355 | 0.013245 |
| rank_one | 7.169 | 15.91 | 45.238 | 0.042294 |
| raw_all | 9.947 | 26.82 | 64.006 | 0.015863 |
| risk_controller | 0.6355 | 0.856 | 0.050674 | 0.084678 |
| sampled_eem | 0.5424 | 0.4274 | 0.071704 | 0.039828 |

## Coverage-aware policy by stratum

| Stratum | Scenes | Accepted | Wrong accepted | Correct rejected |
|---|---:|---:|---:|---:|
| antipodal | 4 | 4 | 0 | 0 |
| bank_missing | 4 | 4 | 0 | 0 |
| easy | 4 | 4 | 0 | 0 |
| low_snr | 4 | 0 | 0 | 4 |
| model_discrepancy | 4 | 0 | 0 | 4 |
| tangent_two_world | 4 | 4 | 4 | 0 |

## Paired controller/fine result

{
  "controller_actions": [
    "full_fine",
    "full_fine",
    "full_fine",
    "full_fine",
    "full_fine",
    "full_fine",
    "full_fine",
    "full_fine"
  ],
  "median_paired_cost_ratio": 1.8614462593778212,
  "paired_ratios": [
    1.7205674132984206,
    2.0386342596569964,
    2.1231034298939364,
    1.9177002996338306,
    1.7440324395606437,
    1.8051922191218117,
    1.7725688413369625,
    2.055452903936487
  ],
  "descriptive_bootstrap_95_interval": [
    1.7440324395606437,
    2.055452903936487
  ],
  "max_parameter_difference": 0.0,
  "note": "Only eight scenes. Timing includes the independent pilot and charged construction/screening; hardware mode-synthesis cost unavailable."
}

## Exact-modal strata

{
  "n": 16,
  "accepted": 12,
  "wrong_accepted": 0,
  "zero_events_one_sided_95_upper_if_iid_population": 0.1707497229824809,
  "interpretation": "Stratified finite demonstration, not an iid population assurance or a calibrated sub-percent empirical failure rate."
}

## Reference correction, fresh sample

{
  "fresh_scenes": 16,
  "median_material_percent": 0.45734111105057024,
  "max_material_visible_norm_without_reference": 4.580667833169738e-06,
  "min_material_visible_norm_with_reference": 18.261856351495958,
  "maximum_no_reference_profile_loss_range": 1.2079226507921703e-13,
  "scope": "geometry known externally; three independent noisy complex gain references"
}
