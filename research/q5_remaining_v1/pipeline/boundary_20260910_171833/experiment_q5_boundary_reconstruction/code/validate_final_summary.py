#!/usr/bin/env python3
"""Validate results/FINAL_SUMMARY.json against the Part I contract."""
import hashlib
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S = json.load(open(os.path.join(ROOT, "results", "FINAL_SUMMARY.json")))
OK = True


def check(label, cond, extra=""):
    global OK
    print(("PASS " if cond else "FAIL ") + label + (" " + str(extra) if extra else ""))
    OK = OK and bool(cond)


for k in ("hypothesis", "overall_verdict", "counts", "constant_table", "artifacts",
          "environment", "declared_limits", "failures_and_caveats"):
    check("top-level key %s" % k, k in S)

t = S["constant_table"]
allowed = {"REPRODUCED", "REPRODUCED_WEAKER", "NOT_REPRODUCED", "INCOMPARABLE"}
check("verdicts all in allowed set", all(r["verdict"] in allowed for r in t),
      sorted({r["verdict"] for r in t}))
check("row keys", all(set(("name", "supplied", "computed", "rel_diff", "verdict")) <= set(r)
                      for r in t))
c = S["counts"]
from collections import Counter
cnt = Counter(r["verdict"] for r in t)
check("counts match table", cnt["REPRODUCED"] == c["REPRODUCED"]
      and cnt["REPRODUCED_WEAKER"] == c["REPRODUCED_WEAKER"]
      and cnt["NOT_REPRODUCED"] == c["NOT_REPRODUCED"]
      and cnt["INCOMPARABLE"] == c["INCOMPARABLE"], dict(cnt))
check("constants_compared == rows", c["constants_compared"] == len(t), (c["constants_compared"], len(t)))
check("leaves == 866", c["leaves_in_inputs_modal_boundary_json"] == 866)
check("declared_limits has 6 verbatim items", len(S["declared_limits"]) == 6)
required_caveats = ("qprime_lower_amplitude_derivative_lower_not_bit_reproduced",
                    "section3_gain_instance_inadmissible",
                    "logdet_vs_gls_distinction",
                    "annulus_constrained_breakdown",
                    "geometry_nuisance_breaks_sigma_min_Avis_bound",
                    "ill_conditioned_LS_solver_artifact")
ids = [f["id"] for f in S["failures_and_caveats"]]
for r in required_caveats:
    check("caveat %s" % r, r in ids)

bad = []
for a in S["artifacts"]:
    p = os.path.join(ROOT, a["path"])
    if not os.path.exists(p):
        bad.append((a["path"], "missing"))
        continue
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    if h != a["sha256"]:
        bad.append((a["path"], "hash mismatch"))
check("all artifact paths exist with matching sha256", not bad, bad)

env = S["environment"]
check("environment keys", {"interpreter", "python_version", "numpy_version",
                           "mpmath_version", "mpmath_provenance"} <= set(env))
check("mpmath provenance sha256 list path", bool(env["mpmath_provenance"]["sha256_list_path"]))
check("mpmath provenance list hash present", bool(env["mpmath_provenance"]["sha256_list_sha256"]))
h = S["part_H_counterexample"]
check("lambda interval", h["admissible_lambda_interval"]["lambda_min_float"] < 1.25)
check("no-rescale counterexamples >= 2", len(h["no_rescale_counterexamples"]) >= 2)
check("all H residual <= 1e-15", h["all_delta_mu_zero_to_1e15"])
g = S["part_G_mie_fix"]
check("grid numerically unchanged", g["differing_leaves_in_numerical_sections"] == 0)
check("three-way grid unchanged flag", g["three_way_grid_unchanged"] is True)
print("OVERALL", "PASS" if OK else "FAIL")
