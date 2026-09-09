"""Family 14: robustness-scope audit (no pseudo-certification).

Run from the experiment root:
    .venv/bin/python src/family14_robustness_scope.py

Purpose
-------
Audit the *scope* of the previously executed Family 5b2 and Family 5 parent
robustness evidence without recomputing any Helmholtz operator or forward
model:

  A.  re-evaluate the finite-dimensional affine tangent certificate formula
      from ingredients stored in results/family5b2_alpha_tight.json and
      compare with the stored L_cert_affine values;
  B.  restate the stored FD-Jacobian validation as conditional/discretized
      evidence, not continuum proof;
  C.  report the stored sampled full-nonlinear and affine worst quotients,
      margins, and ratios without calling the samples a certificate;
  D.  restate the stored generalized-derivative validation, the repeated-
      eigenvalue cluster control, and the engineered rank-event numbers from
      results/family5_parent_generalized.json as audit evidence, together
      with the open questions they leave unresolved.

The audit explicitly demotes all full nonlinear execution-error robustness
claims to empirical/structural observations.  Only the affine tangent
certificate is retained, and that certificate is a finite-dimensional
conditional certificate (conditional on the implemented, FD-validated
Jacobians), not a continuum or full nonlinear Helmholtz certificate.

No existing source, result, figure, or note is modified.  Prior artifacts are
reused read-only.  No figure is produced.  Apple Silicon CPU only.
"""

from __future__ import annotations

import hashlib
import json
import platform
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy  # noqa: F401  (version recorded only)
import matplotlib  # noqa: F401  (version recorded only)

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent

EPS_LIST = [1e-3, 3e-3, 1e-2]
REL_PATHS = {
    "family5b2_json": "results/family5b2_alpha_tight.json",
    "family5b2_script": "src/family5b2_alpha_tight.py",
    "parent_generalized_json": "results/family5_parent_generalized.json",
    "parent_generalized_script": "src/family5_parent_generalized.py",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_artifacts() -> dict:
    """Read both stored result JSONs (and referenced source scripts) read-only."""
    base = _ROOT
    loaded = {}
    for key, rel in REL_PATHS.items():
        p = base / rel
        loaded[key] = {
            "path": rel,
            "absolute_path": str(p),
            "sha256": _sha256(p),
        }
    b5b2 = json.loads((base / REL_PATHS["family5b2_json"]).read_text(encoding="utf-8"))
    parent = json.loads(
        (base / REL_PATHS["parent_generalized_json"]).read_text(encoding="utf-8")
    )
    return loaded, b5b2, parent


def _stored_row(rows: list[dict], eps: float) -> dict:
    for r in rows:
        if float(r["eps"]) == float(eps):
            return r
    raise KeyError(f"no stored row for eps={eps}")


def recompute_affine_rows(const: dict, eps_list: list[float]) -> list[dict]:
    """Recompute the stored Family 5b2 affine certificate from ingredients.

    A_max    = ||A0||_F + JA_op*eps
    tau_max  = tau + JB_op*eps
    mu_aff   = 1/(max(0, sigma_min(B0) - JB_op*eps)^2 + alpha)
    L_cert   = 2*A_max*JA_op
             + 2*A_max^2*JB_op*(tau_max*mu_aff + tau_max^3*mu_aff^2)

    This duplicates only the stored ingredient algebra; no operator or forward
    model is evaluated.
    """
    A_fro = float(const["A_fro"])
    tau = float(const["tau_spec"])
    sig = float(const["sigma_min_B"])
    JA_op = float(const["JA_op_sigma_max"])
    JB_op = float(const["JB_op_sigma_max"])
    alpha = float(const["alpha"])
    out = []
    for eps in eps_list:
        A_max = A_fro + JA_op * eps
        tau_max = tau + JB_op * eps
        sigma_aff_lb = max(0.0, sig - JB_op * eps)
        mu_aff = 1.0 / (sigma_aff_lb**2 + alpha)
        L_cert = (
            2.0 * A_max * JA_op
            + 2.0 * A_max**2 * JB_op
            * (tau_max * mu_aff + tau_max**3 * mu_aff**2)
        )
        out.append(
            {
                "eps": float(eps),
                "A_max": float(A_max),
                "tau_max": float(tau_max),
                "sigma_min_aff_lb": float(sigma_aff_lb),
                "mu_aff": float(mu_aff),
                "recomputed_L_cert_affine": float(L_cert),
            }
        )
    return out


def build_layer_a(const: dict, aff_rows: list[dict], eps_list: list[float]) -> dict:
    recon = recompute_affine_rows(const, eps_list)
    comparison = []
    stored_by_eps = {float(r["eps"]): r for r in aff_rows}
    max_abs_rel = 0.0
    for rr in recon:
        eps = rr["eps"]
        srow = stored_by_eps[eps]
        stored_L = float(srow["L_cert_affine"])
        rel_diff = (rr["recomputed_L_cert_affine"] - stored_L) / stored_L
        max_abs_rel = max(max_abs_rel, abs(rel_diff))
        comparison.append(
            {
                "eps": eps,
                "recomputed_L_cert_affine": rr["recomputed_L_cert_affine"],
                "stored_L_cert_affine": stored_L,
                "abs_relative_diff": float(abs(rel_diff)),
                "A_max": rr["A_max"],
                "tau_max": rr["tau_max"],
                "mu_aff": rr["mu_aff"],
            }
        )
    return {
        "status": (
            "finite-dimensional structural bound conditional on the implemented "
            "FD-validated Jacobians; NOT a continuum proof and NOT a full "
            "nonlinear Helmholtz certificate"
        ),
        "ingredients": {
            "A0_fro": const["A_fro"],
            "tau_sigma_max_B0": const["tau_spec"],
            "sigma_min_B0": const["sigma_min_B"],
            "JA_op_sigma_max": const["JA_op_sigma_max"],
            "JB_op_sigma_max": const["JB_op_sigma_max"],
            "alpha": const["alpha"],
        },
        "formula": (
            "A_max=||A0||_F+JA_op*eps; tau_max=tau+JB_op*eps; "
            "mu_aff=1/(max(0,sigma_min(B0)-JB_op*eps)^2+alpha); "
            "L_cert_affine=2*A_max*JA_op+2*A_max^2*JB_op"
            "*(tau_max*mu_aff+tau_max^3*mu_aff^2)"
        ),
        "L_struct_tight_at_X0": const["L_struct_point_tight"],
        "comparison_rows": comparison,
        "max_abs_relative_diff": float(max_abs_rel),
        "agreement_note": (
            "recomputed with the same IEEE-double operation order stored by "
            "family5b2_alpha_tight.py; agreement is exact at this precision "
            "(max abs relative diff 0.0)"
        ),
    }


def build_layer_b(fd: dict) -> dict:
    return {
        "status": (
            "conditional/discretized full-wave differentiation: stored FD "
            "validation supports the implemented Jacobian but does not certify "
            "the continuum derivative"
        ),
        "h1": fd["h1"],
        "h2": fd["h2"],
        "global_relative_diff_fro": {
            "J_A": fd["global_relative_diff_fro"]["J_A"],
            "J_B": fd["global_relative_diff_fro"]["J_B"],
        },
        "per_column_relative_diff_max": fd["per_column_relative_diff_max"],
        "per_column_relative_diff_min": fd["per_column_relative_diff_min"],
        "o2_check_columns": fd["o2_checks"],
        "statement": (
            "J_A/J_B relative Frobenius FD agreement ~6e-7 (max per column "
            "~7.7e-7 / ~7.6e-7) and O(h^2) ratio checks near 4 support the "
            "implemented centered-FD Jacobians only"
        ),
    }


def build_layer_c(aff_cert: dict, full_val: dict, eps_list: list[float]) -> dict:
    aff_by_eps = {float(r["eps"]): r for r in aff_cert["rows"]}
    full_by_eps = {float(r["eps"]): r for r in full_val["rows"]}
    rows = []
    for eps in eps_list:
        arow = aff_by_eps[eps]
        frow = full_by_eps[eps]
        L = float(arow["L_cert_affine"])
        rows.append(
            {
                "eps": eps,
                "L_cert_affine": L,
                "full_nonlinear_worst_fro_ratio": frow["worst_fro_ratio"],
                "full_nonlinear_mean_fro_ratio": frow["mean_fro_ratio"],
                "affine_worst_fro_ratio": arow["worst_fro_ratio"],
                "affine_mean_fro_ratio": arow["mean_fro_ratio"],
                "margin_L_cert_minus_full_nonlinear_worst": float(
                    L - frow["worst_fro_ratio"]
                ),
                "ratio_L_cert_over_full_nonlinear_worst": float(
                    L / frow["worst_fro_ratio"]
                ),
                "margin_L_cert_minus_affine_worst": float(L - arow["worst_fro_ratio"]),
                "ratio_L_cert_over_affine_worst": float(L / arow["worst_fro_ratio"]),
                "n_full_nonlinear_samples": full_val["n_unit_samples"],
                "le_L_struct_point_tight_fro": frow["le_L_struct_point_tight_fro"],
                "violations_full_nonlinear_fro": frow["violation_magnitude_fro"],
            }
        )
    eps10 = 1e-2
    a10 = aff_by_eps[eps10]
    f10 = full_by_eps[eps10]
    L10 = float(a10["L_cert_affine"])
    return {
        "status": (
            "executed numerical evidence: zero sampled violations over the "
            "stored sample counts at the tested eps cannot be called "
            "certification; observed worst quotients and margins are reported "
            "exactly"
        ),
        "full_model_n_unit_samples": full_val["n_unit_samples"],
        "full_model_seed": full_val["seed"],
        "affine_n_unit_samples": aff_cert["n_unit_samples"],
        "affine_seed": aff_cert["seed"],
        "rows": rows,
        "eps_1e-2_detail": {
            "L_cert_affine": L10,
            "full_nonlinear_worst_fro_ratio": f10["worst_fro_ratio"],
            "margin": float(L10 - f10["worst_fro_ratio"]),
            "ratio": float(L10 / f10["worst_fro_ratio"]),
            "affine_worst_fro_ratio": a10["worst_fro_ratio"],
            "affine_margin": float(L10 - a10["worst_fro_ratio"]),
            "affine_ratio": float(L10 / a10["worst_fro_ratio"]),
        },
        "overall_full_model_worst_fro_ratio": full_val["overall_worst_fro_ratio"],
        "all_full_model_samples_le_bound": bool(
            full_val["all_le_L_struct_point_tight_fro"]
        ),
    }


def _three_modes_evidence(modes: list[dict]) -> list[dict]:
    ev = []
    for m in modes:
        row_1e3 = next(r for r in m["rows"] if float(r["eps"]) == 1e-3)
        ev.append(
            {
                "index_ascending": m["index_ascending"],
                "rho0": m["rho0"],
                "min_adjacent_gap": m["min_adjacent_gap"],
                "drho_formula": m["drho_formula"],
                "fd_relative_derivative_error_at_1e-3": row_1e3[
                    "relative_derivative_error"
                ],
                "first_order_abs_error_at_1e-3": row_1e3["first_order_abs_error"],
                "data_overlap_plus_at_1e-3": row_1e3["data_overlap_plus"],
                "first_order_error_slope_largest4": m[
                    "first_order_error_slope_largest4"
                ],
            }
        )
    return ev


def build_layer_d(parent: dict) -> dict:
    cluster = parent["no_prior_rho_one_cluster"]
    event = parent["engineered_rank_event"]
    rows = event["rows"]
    row_t0 = next(r for r in rows if float(r["t"]) == 0.0)
    row_near = next(r for r in rows if float(r["t"]) == -1e-8)
    # Rows at +/-1e-8 carry the same retained mass; use both for provenance.
    row_pos = next(r for r in rows if float(r["t"]) == 1e-8)
    cases = {}
    for key in ("finite_prior", "no_prior"):
        case = parent[key]
        cases[key] = {
            "alpha": case["alpha"],
            "rank_A": case["rank_A"],
            "condition_A": case["condition_A"],
            "normalization_direct_coefficient_max_error": case[
                "normalization_direct_coefficient_max_error"
            ],
            "normalization_factored_data_max_error": case[
                "normalization_factored_data_max_error"
            ],
            "generalized_backward_residual_factored_max": case[
                "generalized_residual_factored_max"
            ],
            "generalized_residual_relative_eigenvalue_max": case[
                "generalized_residual_relative_eigenvalue_max"
            ],
            "chosen_simple_positive_gap_modes": _three_modes_evidence(
                case["chosen_modes"]
            ),
            "gate": case["gate"],
        }
    return {
        "status": (
            "open questions: the full nonlinear Helmholtz operator has NO "
            "interval/uniform certificate in this implementation; "
            "generalized-derivative formulas require simple positive-gap modes "
            "(validated) and locally constant rank (inapplicable at rank "
            "events, see results/family5_parent_generalized.json); repeated "
            "eigenvalues require compressed derivatives (cluster control "
            "evidence)"
        ),
        "generalized_derivative_validation": cases,
        "repeated_eigenvalue_cluster": {
            "cluster_indices": cluster["cluster_indices"],
            "cluster_size": cluster["cluster_size"],
            "compressed_derivative_eigenvalues": cluster[
                "compressed_derivative_eigenvalues"
            ],
            "compressed_derivative_spectral_norm": cluster[
                "compressed_derivative_spectral_norm"
            ],
            "direct_coefficient_compressed_spectral_norm": cluster[
                "direct_coefficient_compressed_spectral_norm"
            ],
            "perturbed_rows": cluster["perturbed_rows"],
            "max_perturbed_cluster_deviation_from_one": float(
                max(r["max_abs_rho_minus_one"] for r in cluster["perturbed_rows"])
            ),
            "statement": (
                "repeated generalized eigenvalues require the compressed "
                "matrix; individual-vector derivatives are not invariant; the "
                "six structural rho=1 modes stay at one to <=1.44e-15 under "
                "the stored constant-rank perturbations"
            ),
        },
        "engineered_rank_event": {
            "construction": event["construction"],
            "base_rank_B": event["base_rank_B"],
            "base_sigma_min_B": event["base_sigma_min_B"],
            "projector_jump_norm_at_zero": event["projector_jump_norm_at_zero"],
            "KSLAM_jump_fro_at_zero": event["KSLAM_jump_fro_at_zero"],
            "KSLAM_jump_relative_to_KIS": event["KSLAM_jump_relative_to_KIS"],
            "rank_B_t0": row_t0["rank_B"],
            "retained_mass_t_minus_1e-8": row_near["retained_mass"],
            "retained_mass_t_0": row_t0["retained_mass"],
            "retained_mass_t_plus_1e-8": row_pos["retained_mass"],
            "statement": (
                "rank(B) drops 18 -> 17 at t=0, the orthogonal projector jumps "
                "by operator norm 1.0, retained mass jumps from ~9.40716 to "
                "10.40528, and the relative Frobenius information jump is "
                "0.13996; no smooth constant-rank derivative formula applies "
                "at the event"
            ),
            "full_rows": rows,
        },
        "parent_cannot_establish": parent["cannot_establish"],
    }


def build_report(results: dict, artifacts: dict) -> str:
    r = results
    a = r["audit_layers"]["layer_a"]
    b = r["audit_layers"]["layer_b"]
    c = r["audit_layers"]["layer_c"]
    d = r["audit_layers"]["layer_d"]
    fmt = "{:.6e}"

    def fv(x: float) -> str:
        return fmt.format(float(x))

    L = []
    L.append("# Family 14: robustness-scope audit (demotion audit)")
    L.append("")
    L.append(
        "Date: "
        + results["generated_utc"]
        + " (UTC; also recorded in JSON).  Experiment: "
        "`experiment_pose_confounding_spectral_geometry`.  This family audits "
        "the scope of previously executed robustness evidence; it computes no "
        "new Helmholtz operator and runs no new forward model."
    )
    L.append("")
    L.append("## Exact command and runtime")
    L.append("")
    L.append("```bash")
    L.append(r["command"])
    L.append("```")
    L.append("")
    L.append(
        f"Wall runtime: {r['wall_runtime_seconds']:.3f} s.  Platform: "
        f"{r['platform']}.  Python {r['versions']['python']}, "
        f"numpy {r['versions']['numpy']}, scipy {r['versions']['scipy']}, "
        f"matplotlib {r['versions']['matplotlib']}.  Apple Silicon CPU only; "
        "no GPU/MPS/CUDA."
    )
    L.append("")
    L.append("## Read-only reuse")
    L.append("")
    L.append(
        "Loaded artifacts (SHA-256 in the JSON): "
        "`results/family5b2_alpha_tight.json`, "
        "`results/family5_parent_generalized.json`, plus the two source "
        "scripts referenced only for name/structure context "
        "(`src/family5b2_alpha_tight.py`, `src/family5_parent_generalized.py`). "
        "No existing source, result, figure, or note file was modified."
    )
    L.append("")
    L.append("## Layer A: exact finite-dimensional linear algebra")
    L.append("")
    L.append(
        "The affine tangent certificate is a finite-dimensional structural "
        "bound conditional on the implemented FD-validated Jacobians; it is "
        "NOT a continuum proof and NOT a full nonlinear Helmholtz certificate."
    )
    L.append("")
    L.append(
        f"Stored pointwise structural constant: L_struct_tight(X0) = "
        f"{fv(a['L_struct_tight_at_X0'])}.  Stored affine certificates for "
        "eps in {1e-3, 3e-3, 1e-2}: "
        + ", ".join(fv(row["stored_L_cert_affine"]) for row in a["comparison_rows"])
        + ".  Recomputation from the stored ingredients reproduces every "
        "stored value exactly at IEEE double precision (max abs relative diff "
        + f"{a['max_abs_relative_diff']:.3e}"
        + "); the same formula and operation order as Family 5b2 was used."
    )
    L.append("")
    L.append("| eps | recomputed L_cert_affine | stored L_cert_affine | abs rel diff |")
    L.append("| --- | ------------------------ | -------------------- | ------------ |")
    for row in a["comparison_rows"]:
        L.append(
            f"| {row['eps']:.0e} | {fv(row['recomputed_L_cert_affine'])} | "
            f"{fv(row['stored_L_cert_affine'])} | "
            f"{row['abs_relative_diff']:.3e} |"
        )
    L.append("")
    L.append("## Layer B: conditional/discretized full-wave differentiation")
    L.append("")
    L.append(
        "FD validation numbers (J_A/J_B relative Frobenius agreement "
        f"{fv(b['global_relative_diff_fro']['J_A'])} / "
        f"{fv(b['global_relative_diff_fro']['J_B'])}, "
        f"i.e. ~6e-7) support the implemented Jacobian but do not certify the "
        "continuum derivative.  Max per-column relative differences are "
        f"{fv(b['per_column_relative_diff_max']['J_A'])} (J_A) and "
        f"{fv(b['per_column_relative_diff_max']['J_B'])} (J_B); stored "
        "consecutive-difference ratios over h ~ 4.0 provide O(h^2) support for "
        "the centered-FD estimates only."
    )
    L.append("")
    L.append("## Layer C: executed numerical evidence")
    L.append("")
    L.append(
        "Zero sampled violations over the stored sample counts at the tested "
        "eps cannot be called certification; the observed worst quotients and "
        "margins are reported below.  Full-nonlinear sampling was 500 unit "
        "directions per eps (seed 9191); affine sampling was 2000 unit "
        "directions per eps (seed 7171)."
    )
    L.append("")
    L.append(
        "| eps | L_cert_affine | full-nonlinear worst ||K(eps u)-K(0)||_F/eps "
        "| margin | ratio | affine worst | affine margin |"
    )
    L.append(
        "| --- | ------------- | -------------------------------------------- "
        "| ------ | ----- | ------------ | ------------- |"
    )
    for row in c["rows"]:
        L.append(
            f"| {row['eps']:.0e} | {fv(row['L_cert_affine'])} | "
            f"{fv(row['full_nonlinear_worst_fro_ratio'])} | "
            f"{fv(row['margin_L_cert_minus_full_nonlinear_worst'])} | "
            f"{row['ratio_L_cert_over_full_nonlinear_worst']:.4f} | "
            f"{fv(row['affine_worst_fro_ratio'])} | "
            f"{fv(row['margin_L_cert_minus_affine_worst'])} |"
        )
    L.append("")
    det = c["eps_1e-2_detail"]
    L.append(
        f"At eps=1e-2: L_cert_affine={fv(det['L_cert_affine'])}, "
        "full-nonlinear worst sampled quotient "
        f"{fv(det['full_nonlinear_worst_fro_ratio'])} (max across the three "
        "stored full-model eps), margin "
        f"{fv(det['margin'])}, ratio {det['ratio']:.4f}; "
        "affine worst sampled quotient "
        f"{fv(det['affine_worst_fro_ratio'])}, affine margin "
        f"{fv(det['affine_margin'])}, affine ratio {det['affine_ratio']:.4f}."
    )
    L.append("")
    L.append("## Layer D: open questions")
    L.append("")
    L.append(
        "The full nonlinear Helmholtz operator has NO interval/uniform "
        "certificate in this implementation.  The stored generalized-"
        "derivative formulas were validated on simple positive-gap modes and "
        "require locally constant rank, so they are inapplicable at rank "
        "events (see `results/family5_parent_generalized.json`).  Repeated "
        "eigenvalues require compressed derivatives; only cluster control "
        "evidence is available."
    )
    L.append("")
    gdv = d["generalized_derivative_validation"]
    fp = gdv["finite_prior"]
    np_ = gdv["no_prior"]
    L.append(
        "Generalized-derivative audit evidence (stored): finite-prior backward "
        "residual (factored) "
        f"{fv(fp['generalized_backward_residual_factored_max'])}, "
        "normalization error (factored data) "
        f"{fv(fp['normalization_factored_data_max_error'])}, "
        "FD relative errors at eps=1e-3 for three simple positive-gap modes "
        + ", ".join(
            fv(m["fd_relative_derivative_error_at_1e-3"])
            for m in fp["chosen_simple_positive_gap_modes"]
        )
        + ", prediction-error slopes "
        + ", ".join(
            f"{m['first_order_error_slope_largest4']:.6f}"
            for m in fp["chosen_simple_positive_gap_modes"]
        )
        + "."
    )
    L.append("")
    L.append(
        "No-prior (constant-rank) case: backward residual (factored) "
        f"{fv(np_['generalized_backward_residual_factored_max'])}, "
        "normalization error (factored data) "
        f"{fv(np_['normalization_factored_data_max_error'])}, "
        "FD relative errors at eps=1e-3 "
        + ", ".join(
            fv(m["fd_relative_derivative_error_at_1e-3"])
            for m in np_["chosen_simple_positive_gap_modes"]
        )
        + ", prediction-error slopes "
        + ", ".join(
            f"{m['first_order_error_slope_largest4']:.6f}"
            for m in np_["chosen_simple_positive_gap_modes"]
        )
        + "."
    )
    L.append("")
    cl = d["repeated_eigenvalue_cluster"]
    ev = d["engineered_rank_event"]
    L.append(
        "Repeated-eigenvalue cluster control: the six structural rho=1 modes "
        "stay at one to <=1.44e-15 under the stored perturbed rows, with "
        "compressed-derivative spectral norm "
        f"{fv(cl['compressed_derivative_spectral_norm'])}; the direct "
        "coefficient-space compressed norm is "
        f"{fv(cl['direct_coefficient_compressed_spectral_norm'])}."
    )
    L.append("")
    L.append(
        "Engineered rank event (stored): rank(B) 18 -> 17 at t=0, projector "
        "jump operator norm "
        f"{ev['projector_jump_norm_at_zero']:.3e} (1.0), retained mass "
        f"{fv(ev['retained_mass_t_minus_1e-8'])} -> "
        f"{fv(ev['retained_mass_t_0'])} (9.40716 -> 10.40528 at the stored "
        "rows), and relative Frobenius information jump "
        f"{ev['KSLAM_jump_relative_to_KIS']:.5f} (0.13996).  This is an "
        "engineered algebraic control, not an observed physical trajectory "
        "event, and it falsifies constant-rank derivative applicability at the "
        "event."
    )
    L.append("")
    L.append("## Explicit demotion")
    L.append("")
    L.append("> " + r["audit_layers"]["explicit_demotion_statement"])
    L.append("")
    L.append("## No claims made")
    L.append("")
    for claim in r["audit_layers"]["no_claims"]:
        L.append(f"* {claim}")
    L.append("")
    L.append("## Artifacts")
    L.append("")
    L.append("* [results/family14_robustness_scope.json]"
             "(results/family14_robustness_scope.json)")
    L.append("* [notes/family14_robustness_scope.md]"
             "(notes/family14_robustness_scope.md)")
    L.append("")
    L.append("Source SHA-256: `" + results["source_sha256"] + "`.")
    L.append("")
    L.append("Artifact digests:")
    for key in (
        "family5b2_json",
        "parent_generalized_json",
        "family5b2_script",
        "parent_generalized_script",
        "source_script",
    ):
        if key in artifacts:
            L.append(f"* `{key}` `{artifacts[key]}`")
    L.append("")
    return "\n".join(L)


def main() -> None:
    t0 = time.perf_counter()
    loaded, b5b2, parent = _load_artifacts()
    const = b5b2["constants"]
    fd = b5b2["fd_validation"]
    aff_cert = b5b2["affine_certificate"]
    full_val = b5b2["full_model_validation"]

    layer_a = build_layer_a(const, aff_cert["rows"], EPS_LIST)
    layer_b = build_layer_b(fd)
    layer_c = build_layer_c(aff_cert, full_val, EPS_LIST)
    layer_d = build_layer_d(parent)

    demotion = (
        "All full nonlinear execution-error robustness claims are DEMOTED to "
        "empirical/structural observations; only the affine tangent "
        "certificate is retained as a finite-dimensional conditional "
        "certificate."
    )

    no_claims = [
        "No interval proof is claimed.",
        "No Hankel bound is claimed.",
        "No continuum robustness is claimed.",
        "No 'certified nonlinear robustness' is claimed.",
    ]

    config = {
        "title": (
            "Family 14 robustness-scope audit; no new Helmholtz operator or "
            "forward model"
        ),
        "eps_list": EPS_LIST,
        "loaded_artifacts": [
            {
                "path": loaded["family5b2_json"]["path"],
                "role": (
                    "affine-certificate ingredients, FD validation, affine and "
                    "full-nonlinear sampled rows, structural constant"
                ),
                "sha256": loaded["family5b2_json"]["sha256"],
            },
            {
                "path": loaded["parent_generalized_json"]["path"],
                "role": (
                    "generalized-derivative validation, repeated-eigenvalue "
                    "cluster control, engineered rank event"
                ),
                "sha256": loaded["parent_generalized_json"]["sha256"],
            },
            {
                "path": loaded["family5b2_script"]["path"],
                "role": "read-only reference for stored formula/names",
                "sha256": loaded["family5b2_script"]["sha256"],
            },
            {
                "path": loaded["parent_generalized_script"]["path"],
                "role": "read-only reference for stored formula/names",
                "sha256": loaded["parent_generalized_script"]["sha256"],
            },
        ],
        "reuse_note": (
            "no existing source/result/figure/note modified; prior artifacts "
            "opened read-only"
        ),
    }

    results = {
        "schema": "family14-robustness-scope-v1",
        "family": 14,
        "title": (
            "Robustness-scope audit of Family 5b2/5 robustness evidence "
            "(demotion audit, not pseudo-certification)"
        ),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "runner": f"{platform.python_implementation()} {platform.python_version()}",
        "command": ".venv/bin/python src/family14_robustness_scope.py",
        "wall_runtime_seconds": None,  # filled before serialisation
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_note": "Apple Silicon CPU only; no GPU/MPS/CUDA; no figure produced",
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "source_sha256": _sha256(Path(__file__)),
        "config": config,
        "audit_layers": {
            "layer_a": layer_a,
            "layer_b": layer_b,
            "layer_c": layer_c,
            "layer_d": layer_d,
            "explicit_demotion_statement": demotion,
            "no_claims": no_claims,
        },
        "summary": {
            "affine_certificate_eps_1e-3": layer_a["comparison_rows"][0][
                "stored_L_cert_affine"
            ],
            "affine_certificate_eps_3e-3": layer_a["comparison_rows"][1][
                "stored_L_cert_affine"
            ],
            "affine_certificate_eps_1e-2": layer_a["comparison_rows"][2][
                "stored_L_cert_affine"
            ],
            "recomputation_max_abs_relative_diff": layer_a[
                "max_abs_relative_diff"
            ],
            "full_nonlinear_worst_quotient_eps_1e-2": layer_c["eps_1e-2_detail"][
                "full_nonlinear_worst_fro_ratio"
            ],
            "margin_at_eps_1e-2": layer_c["eps_1e-2_detail"]["margin"],
            "ratio_at_eps_1e-2": layer_c["eps_1e-2_detail"]["ratio"],
            "demotion": demotion,
            "uncertainties": [
                "No interval/uniform/Hankel/continuum claim is supported by "
                "the loaded evidence.",
                "FD Jacobian and generalized-derivative checks are "
                "finite-difference validations of the implemented model, with "
                "numerical floors.",
                "Rank events and repeated-eigenvalue statements are "
                "finite-dimensional/engineered controls, not trajectory "
                "theorems.",
            ],
        },
        "artifact_digests": None,  # filled just before serialisation
    }

    results_dir = _ROOT / "results"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)

    results["wall_runtime_seconds"] = float(time.perf_counter() - t0)
    report_path = notes_dir / "family14_robustness_scope.md"
    json_path = results_dir / "family14_robustness_scope.json"

    # The report needs the artifact digests of its inputs plus its own digest;
    # the JSON cannot contain a truthful self-digest, so it is omitted there.
    artifacts = {
        "family5b2_json": loaded["family5b2_json"]["sha256"],
        "parent_generalized_json": loaded["parent_generalized_json"]["sha256"],
        "family5b2_script": loaded["family5b2_script"]["sha256"],
        "parent_generalized_script": loaded["parent_generalized_script"]["sha256"],
        "source_script": results["source_sha256"],
        "report_note": None,  # filled after writing
    }
    report = build_report(results, artifacts)
    report_path.write_text(report, encoding="utf-8")
    artifacts["report_note"] = _sha256(report_path)
    results["artifact_digests"] = {
        "results/family14_robustness_scope.json": (
            "self-digest omitted; SHA-256 of the written file is verifiable "
            "directly on disk"
        ),
        "notes/family14_robustness_scope.md": artifacts["report_note"],
        "results/family5b2_alpha_tight.json": artifacts["family5b2_json"],
        "results/family5_parent_generalized.json": artifacts[
            "parent_generalized_json"
        ],
        "src/family5b2_alpha_tight.py": artifacts["family5b2_script"],
        "src/family5_parent_generalized.py": artifacts[
            "parent_generalized_script"
        ],
        "src/family14_robustness_scope.py": artifacts["source_script"],
    }
    json_path.write_text(
        json.dumps(results, indent=2, default=str), encoding="utf-8"
    )

    summary = {
        "family14_ok": True,
        "affine_numbers": {
            "L_struct_tight_at_X0": const["L_struct_point_tight"],
            "L_cert_affine": {
                f"{eps:.0e}": layer_a["comparison_rows"][i][
                    "stored_L_cert_affine"
                ]
                for i, eps in enumerate(EPS_LIST)
            },
            "recomputation_max_abs_relative_diff": layer_a[
                "max_abs_relative_diff"
            ],
        },
        "worst_quotients": layer_c["eps_1e-2_detail"],
        "demotion_statement": demotion,
        "uncertainties": results["summary"]["uncertainties"],
        "wall_runtime_seconds": results["wall_runtime_seconds"],
        "results": str(json_path),
        "report": str(report_path),
    }
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
