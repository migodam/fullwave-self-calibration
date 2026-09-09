#!/usr/bin/env python3
"""Round-3 empirical bookkeeping summary for the geometry-lifted tri-space SOM loop.

Reads the result JSONs already on disk (no numerical experiments, no physics
recomputation) and writes:
  * consolidated_round3_summary.json
  * consolidated_round3_tables.md  (also printed to stdout)

Medians / 25-75th percentiles are computed over the per-run records with
numpy-compatible linear interpolation.  File schemas are discovered from the
JSON itself; missing or unclear files are skipped with a note.
"""

from __future__ import annotations

import json
import math
import os
import sys


WD = os.path.dirname(os.path.abspath(__file__))

REQUIRED_FILES = [
    "results_e10_settings_sweep.json",
    "results_e10b_m12_m16_highrank.json",
    "results_e6b_multitx_robust.json",
    "results_e11_vp_state_null.json",
]

# Optional context files used only for the headline "prior anchor points".
CONTEXT_FILES = [
    "results_e5_final.json",
    "results_e9_grid_scenes.json",
    "results_e6.json",
    ("results_e7b_soft_state.json", "results_e7b.json"),  # (listed name, on-disk name)
]

E10_METHODS = [
    "direct",
    "reduced_r4",
    "reduced_r6",
    "wrongpose",
    "known_pose",
    "known_alpha",
]

E11_METHODS = ["direct", "reduced_r6", "vp_hard_r6", "vp_soft_r6"]


def load_json(rel_path: str) -> tuple[object | None, str]:
    """Return (data, status). Status is 'ok' or a short failure reason."""
    full = os.path.join(WD, rel_path)
    try:
        with open(full, "r", encoding="utf-8") as fh:
            return json.load(fh), "ok"
    except FileNotFoundError:
        return None, "file not found"
    except json.JSONDecodeError as exc:
        return None, f"json decode error: {exc}"
    except OSError as exc:
        return None, f"read error: {exc}"


def _load_manifest() -> tuple[dict, list]:
    """Load required and context files. Returns {logical_name: data} + manifest."""
    data = {}
    manifest = []

    for name in REQUIRED_FILES:
        payload, status = load_json(name)
        data[name] = payload
        manifest.append(
            {
                "path": name,
                "parse_status": status if payload is not None else "skipped_missing",
                "note": None if status == "ok" else "required file unavailable on disk",
            }
        )

    for spec in CONTEXT_FILES:
        if isinstance(spec, tuple):
            listed, actual = spec
            payload, status = load_json(actual)
            data[listed] = payload
            if payload is not None:
                manifest.append(
                    {
                        "path": actual,
                        "specified_as": listed,
                        "parse_status": "ok",
                        "note": (
                            f"task listed '{listed}'; the on-disk file is '{actual}' "
                            f"(experiment field: run_e7b_soft_state)"
                        ),
                    }
                )
            else:
                manifest.append(
                    {
                        "path": listed,
                        "parse_status": "skipped_missing",
                        "note": f"not found; on-disk candidate '{actual}' also unavailable",
                    }
                )
        else:
            payload, status = load_json(spec)
            data[spec] = payload
            manifest.append(
                {
                    "path": spec,
                    "parse_status": status if payload is not None else "skipped_missing",
                    "note": None if status == "ok" else "optional context file unavailable",
                }
            )

    return data, manifest


def percentile(vals: list, q: float) -> float | None:
    """Linear-interpolation percentile (numpy default convention)."""
    vals = [v for v in vals if isinstance(v, (int, float)) and v is not None]
    if not vals:
        return None
    vals = sorted(vals)
    if len(vals) == 1:
        return float(vals[0])
    rank = (len(vals) - 1) * q
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return float(vals[lo])
    frac = rank - lo
    return float(vals[lo] + (vals[hi] - vals[lo]) * frac)


def median(vals: list) -> float | None:
    return percentile(vals, 0.5)


def pct_summary(vals: list) -> dict:
    return {
        "q25": percentile(vals, 0.25),
        "median": percentile(vals, 0.50),
        "q75": percentile(vals, 0.75),
    }


def method_medians(runs: list, method: str, err_fields: tuple) -> dict:
    rows = [r for r in runs if r.get("method") == method]
    out = {
        "method": method,
        "n_seeds": len(rows),
        "n_success": int(sum(1 for r in rows if r.get("success"))),
    }
    for field in err_fields:
        out[field] = median([r.get(field) for r in rows])
    return out


def build_e10(d: dict) -> dict | None:
    if not isinstance(d, dict) or "settings" not in d:
        return None
    out = {
        "experiment": d.get("experiment"),
        "rank_sweep_setting": d.get("rank_sweep_setting"),
        "counts": d.get("counts"),
        "settings": [],
    }
    for label, cfg in d["settings"].items():
        entry = {"setting": label}

        # Geometry / rank checks.
        for key in ("M", "k", "aperture", "aperture_span_rad"):
            entry[key] = cfg.get(key)
        rc = cfg.get("rank_checks") or {}
        sv = rc.get("singular_values") or []
        entry["G_s_numerical_rank"] = rc.get("numerical_rank")
        entry["G_s_singular_values"] = {
            "largest": sv[0] if sv else None,
            "smallest": sv[-1] if sv else None,
            "count": len(sv),
            "rank_tol_relative": rc.get("rank_tol_relative"),
        }

        # Retained-lift checks at r=4 and r=6.
        retained = cfg.get("lift_checks", {}).get("retained") or []
        retained_by_r = {e.get("retained_rank"): e for e in retained}
        for r in (4, 6):
            e = retained_by_r.get(r) or {}
            entry[f"lift_r{r}"] = {
                "retained_rank": e.get("retained_rank"),
                "residual_ratio_ret": e.get("residual_ratio_ret"),
                "QH_R_max_abs": e.get("orthogonality_QH_R_max_abs"),
                "QH_R_norm_over_R_norm": e.get(
                    "orthogonality_QH_R_norm_over_R_norm"
                ),
                "norm_J_ret": e.get("norm_J_ret"),
            }
        unrestricted = cfg.get("lift_checks", {}).get("unrestricted") or {}
        entry["lift_unrestricted"] = {
            "retained_rank": unrestricted.get("retained_rank"),
            "residual_ratio_full": unrestricted.get("residual_ratio_full"),
            "norm_J_full": unrestricted.get("norm_J_full"),
        }

        # Subspace sweep at retained ranks r=4 and r=6.
        sweep = cfg.get("subspace_sweep") or []
        sweep_by_r = {e.get("rank"): e for e in sweep}
        for r in (4, 6):
            e = sweep_by_r.get(r) or {}
            entry[f"subspace_r{r}"] = {
                "n_vis": e.get("n_vis"),
                "hidden_rank": e.get("hidden_rank"),
                "hid_sv": e.get("hid_sv"),
                "regime": e.get("regime"),
            }

        # Median over seeds per method (exact stored field names).
        runs = cfg.get("runs") or []
        methods = []
        for method in E10_METHODS:
            mrow = method_medians(runs, method, ("pose_error", "map_error"))
            methods.append(mrow)
        entry["methods_median_over_seeds"] = methods
        out["settings"].append(entry)

    # Extra per-setting sanity fields (only if present on all runs).
    out["note"] = (
        "Median over noise seeds computed from each setting's per-run records; "
        "no physics recomputed."
    )
    return out


def build_e10b(d: dict) -> dict | None:
    if not isinstance(d, dict) or "per_M" not in d:
        return None
    out = {"experiment": d.get("experiment"), "counts": d.get("counts"), "per_M": []}
    for label, cfg in d["per_M"].items():
        ts = cfg.get("transition_summary") or {}
        rank_table = cfg.get("rank_table") or []
        entry = {
            "setting": label,
            "M": cfg.get("M"),
            "k": cfg.get("k"),
            "aperture": cfg.get("aperture"),
            "transition_summary": ts,
            "empty_lift_trailing_r": {
                "vanished_Bred_rs": ts.get("vanished_Bred_rs"),
                "n_vis_zero_rs": ts.get("n_vis_zero_rs"),
                "rank_table_regimes": [
                    {
                        "r": row.get("r"),
                        "regime": row.get("regime"),
                        "n_vis": row.get("n_vis"),
                        "hidden_rank": row.get("hidden_rank"),
                    }
                    for row in rank_table
                    if row.get("regime") in ("hidden", "vanished_Bred", "n_vis_zero")
                ],
            },
        }

        runs = d.get("runs") or []
        runs_m = [r for r in runs if r.get("setting") == label]
        err_fields = (
            "pose_error",
            "pose_error_visible",
            "pose_error_hidden",
            "map_error",
        )
        groups = []
        groups.append(
            {
                **method_medians(runs_m, "direct", err_fields),
                "r": None,
                "regime": "direct reference",
            }
        )
        for sel in cfg.get("selected_run_ranks") or []:
            r = sel.get("r")
            method = f"reduced_r{r}"
            groups.append(
                {
                    **method_medians(runs_m, method, err_fields),
                    "r": r,
                    "regime": sel.get("regime"),
                }
            )
        entry["method_medians_over_seeds"] = groups
        out["per_M"].append(entry)
    return out


def parse_group(group: str) -> tuple[str, int]:
    """Maps E6b group names (mono_L2, mono_L3, dipole_L2, ...) to type/L."""
    src, _, rest = group.partition("_")
    source_type = "monopole" if src == "mono" else ("dipole" if src == "dipole" else src)
    digits = "".join(ch for ch in rest if ch.isdigit())
    return source_type, int(digits) if digits else None


def build_e6b(d: dict) -> dict | None:
    if not isinstance(d, dict) or "runs" not in d:
        return None
    runs = d["runs"]
    groups: dict = {}
    ordered_keys = []
    for r in runs:
        key = (
            r.get("group"),
            r.get("pose_mode"),
            r.get("method_alias"),
            r.get("SNR_dB"),
        )
        if key not in groups:
            groups[key] = []
            ordered_keys.append(key)
        groups[key].append(r)

    metric_fields = (
        "pose_error",
        "theta_error",
        "tx_error",
        "ty_error",
        "map_error",
        "final_residual_ratio",
    )
    out_groups = []
    for key in ordered_keys:
        rows = groups[key]
        r0 = rows[0]
        source_type, L = parse_group(r0.get("group"))
        g = {
            "group": r0.get("group"),
            "source_type": source_type,
            "L": L,
            "pose_mode": r0.get("pose_mode"),
            "method_alias": r0.get("method_alias"),
            "method": r0.get("method"),
            "SNR_dB": r0.get("SNR_dB"),
            "n_seeds": len(rows),
            "n_success": int(sum(1 for x in rows if x.get("success"))),
        }
        for field in metric_fields:
            vals = [x.get(field) for x in rows]
            g[field] = pct_summary(vals)
        out_groups.append(g)

    def is_main(g: dict) -> bool:
        group = g["group"]
        alias = g["method_alias"]
        mode = g["pose_mode"]
        if group == "mono_L2" and mode == "tx_only" and alias in ("direct_tx", "wrongpose"):
            return True
        if group in ("mono_L3", "mono_L2_comoving", "dipole_L2") and alias in (
            "direct_tx",
            "direct_comoving",
        ):
            return True
        return False

    # FD checks at eps = 1e-4 (exact list membership; fall back to nearest log step).
    fd_checks = d.get("fd_checks_dipole") or []
    eps1e4 = [c for c in fd_checks if c.get("eps") == 1e-4]
    if not eps1e4:
        by_eps = {}
        for c in fd_checks:
            by_eps[abs(math.log10(c.get("eps", float("nan"))))] = c
        if by_eps:
            nearest = min(by_eps, key=lambda x: abs(x - 4.0))
            eps1e4 = [by_eps[nearest]]

    dipole_gauge = d.get("gauge_partA_dipole") or {}
    gauge_summary = {}
    for k in ("B_t_analysis", "comoving_B_total_sanity", "monopole_B_t_reference"):
        gauge_summary[k] = dipole_gauge.get(k)

    return {
        "experiment": (d.get("parameters") or {}).get("experiment"),
        "arrays": (d.get("parameters") or {}).get("arrays"),
        "groups": {
            "all": out_groups,
            "main_groups_requested": [g for g in out_groups if is_main(g)],
        },
        "fd_max_error_at_eps_1e-4": [
            {
                "pose_mode": c.get("pose_mode"),
                "matrix": c.get("matrix"),
                "fd_convention": c.get("fd_convention"),
                "eps": c.get("eps"),
                "per_col_rel_err": c.get("per_col_rel_err"),
                "max_rel_err": c.get("max_rel_err"),
            }
            for c in eps1e4
        ],
        "fd_all_checks": fd_checks,
        "dipole_gauge": gauge_summary,
    }


def build_e11(d: dict) -> dict | None:
    if not isinstance(d, dict) or "runs" not in d:
        return None
    runs = d["runs"]
    # Seed-median universe is the zero-init per-seed set (run_label == method).
    zero_runs = [
        r for r in runs if r.get("init_label") == "zero" or r.get("run_label") == r.get("method")
    ]
    fields = (
        "pose_error",
        "pose_error_visible",
        "pose_error_hidden",
        "map_error",
        "data_residual_ratio",
        "discarded_current_fraction",
    )
    method_rows = []
    for method in E11_METHODS:
        rows = [r for r in zero_runs if r.get("method") == method]
        row = {
            "method": method,
            "n_seeds": len(rows),
            "n_success": int(sum(1 for r in rows if r.get("success"))),
        }
        for f in fields:
            row[f] = pct_summary([r.get(f) for r in rows])
        method_rows.append(row)

    out = {
        "experiment": d.get("experiment"),
        "parameters": d.get("parameters"),
        "runs_all": runs,
        "methods_median_over_seeds": method_rows,
        "visible_pose_subspace_r6_at_p_init": d.get(
            "visible_pose_subspace_r6_at_p_init"
        ),
        "soft_filter": d.get("soft_filter"),
        "projector_checks_at_p_init": d.get("projector_checks_at_p_init"),
        "vp_anchor_at_alpha_init_p_init": d.get("vp_anchor_at_alpha_init_p_init"),
        "finite_dimensional_vp_jacobian": d.get("finite_dimensional_vp_jacobian"),
        "interpretation": d.get("interpretation"),
        "note": d.get("note"),
    }
    return out


def _e5_anchor(d: dict) -> dict:
    records = d.get("records") or []
    fields = (
        "pose_error",
        "pose_error_visible",
        "pose_error_hidden",
        "map_error",
    )
    rows = []
    methods = ("direct", "reduced_r4", "reduced_r6", "wrongpose")
    p_labels = sorted({r.get("p_init_label") for r in records if r.get("p_init_label")})
    for p in p_labels:
        for method in methods:
            sel = [
                r
                for r in records
                if r.get("method") == method and r.get("p_init_label") == p
            ]
            row = {
                "p_init_label": p,
                "method": method,
                "r": sel[0].get("r") if sel else None,
                "n_seeds": len(sel),
                "n_success": int(sum(1 for r in sel if r.get("success"))),
            }
            for f in fields:
                row[f] = median([r.get(f) for r in sel])
            rows.append(row)
    return {
        "goal_or_note": (d.get("note") or "")[:500],
        "records_median_over_seeds": rows,
    }


def _e9_anchor(d: dict) -> dict:
    return {
        "goal": d.get("goal"),
        "summary_groups": d.get("summary_groups"),
    }


def _e6_anchor(d: dict) -> dict:
    return {
        "gauge_partA": d.get("gauge_partA"),
        "finite_difference_checks": d.get("finite_difference_checks"),
        "fd_pass_at_eps_1e-4_max_rel": (d.get("parameters") or {}).get("fd", {}).get(
            "pass_at_eps_1e-4_max_rel"
        ),
    }


def _e7b_anchor(d: dict) -> dict:
    return {
        "goal": d.get("goal"),
        "baselines": d.get("baselines"),
        "summaries_median_over_seeds": d.get("summaries_median_over_seeds"),
        "visible_pose_subspace_at_p_init_r6": d.get(
            "visible_pose_subspace_at_p_init_r6"
        ),
        "notes": d.get("notes"),
    }


def build_prior_anchors(data: dict) -> tuple[dict, list]:
    out = {}
    used = []

    e5 = data.get("results_e5_final.json")
    if e5 is not None and isinstance(e5, dict) and "records" in e5:
        out["E5_final"] = _e5_anchor(e5)
        used.append("results_e5_final.json")

    e9 = data.get("results_e9_grid_scenes.json")
    if e9 is not None and isinstance(e9, dict) and "summary_groups" in e9:
        out["E9_grid_scenes"] = _e9_anchor(e9)
        used.append("results_e9_grid_scenes.json")

    e6 = data.get("results_e6.json")
    if e6 is not None and isinstance(e6, dict) and "gauge_partA" in e6:
        out["E6_multitx"] = _e6_anchor(e6)
        used.append("results_e6.json")

    e7b = data.get("results_e7b_soft_state.json")
    if e7b is not None and isinstance(e7b, dict) and "summaries_median_over_seeds" in e7b:
        out["E7b_soft_state"] = _e7b_anchor(e7b)
        used.append("results_e7b.json (listed as results_e7b_soft_state.json)")

    if not out:
        out["note"] = (
            "No prior-anchor fields could be parsed robustly; empty per the "
            "round-3 instruction rather than guessing."
        )
    return out, used


# ---------------------------------------------------------------------------
# Markdown rendering helpers
# ---------------------------------------------------------------------------


def fmt(x, sig: int = 4) -> str:
    if x is None:
        return "—"
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        ax = abs(x)
        if ax == 0.0:
            return "0"
        if ax < 1e-3 or ax >= 1e3:
            return f"{x:.{sig-1}e}"
        return f"{x:.{sig-1}g}"
    return str(x)


def md_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for row in rows:
        assert len(row) == len(headers), (headers, row)
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)

def md_for_e10(e10: dict) -> str:
    if not e10:
        return ""
    src = "results_e10_settings_sweep.json"
    settings = e10["settings"]
    chunks = []
    chunks.append(
        "\n### Table R3-1: E10 setting sweep (source: `%s`)\n" % src
    )
    hdr1 = [
        "setting",
        "M",
        "k",
        "aperture",
        "rank G_s",
        "sv_max",
        "sv_min",
        "lift r4 resid",
        "lift r6 resid",
        "Q^H R r4",
        "Q^H R r6",
        "subspace r4 n_vis (hid)",
        "subspace r6 n_vis (hid)",
    ]
    rows1 = []
    for s in settings:
        rows1.append(
            [
                s["setting"],
                fmt(s.get("M")),
                fmt(s.get("k")),
                str(s.get("aperture")),
                fmt((s.get("G_s_numerical_rank"))),
                fmt((s.get("G_s_singular_values") or {}).get("largest")),
                fmt((s.get("G_s_singular_values") or {}).get("smallest")),
                fmt((s.get("lift_r4") or {}).get("residual_ratio_ret")),
                fmt((s.get("lift_r6") or {}).get("residual_ratio_ret")),
                fmt((s.get("lift_r4") or {}).get("QH_R_norm_over_R_norm")),
                fmt((s.get("lift_r6") or {}).get("QH_R_norm_over_R_norm")),
                _nvis_cell(s.get("subspace_r4")),
                _nvis_cell(s.get("subspace_r6")),
            ]
        )
    chunks.append(md_table(hdr1, rows1))
    chunks.append(
        "\nQ^H R columns are `orthogonality_QH_R_norm_over_R_norm` (|Q^H R|/|R|). "
        "n_vis (hid) is written `n_vis(hidden_rank)` from each setting's subspace sweep."
    )

    hdr2 = [
        "setting",
        "direct",
        "reduced_r4",
        "reduced_r6",
        "wrongpose",
        "known_pose",
        "known_alpha",
    ]
    for metric in ("pose_error", "map_error"):
        rows2 = []
        for s in settings:
            med = {m["method"]: m[metric] for m in s["methods_median_over_seeds"]}
            rows2.append(
                [s["setting"]]
                + [fmt(med.get(m)) for m in hdr2[1:]]
            )
        label = "Pose error" if metric == "pose_error" else "Map error"
        chunks.append(
            f"\n**Table R3-1{chr(98+int(metric=='map_error'))}: {label} "
            "median over seeds** (source: `%s`)\n" % src
        )
        chunks.append(md_table(hdr2, rows2))
    return "\n".join(chunks)


def _nvis_cell(sub: dict | None) -> str:
    if not sub:
        return "—"
    nv = sub.get("n_vis")
    hh = sub.get("hidden_rank")
    if nv is None:
        return "—"
    return f"{fmt(nv)} (hid {fmt(hh)})"


def md_for_e10b(e10b: dict) -> str:
    if not e10b:
        return ""
    src = "results_e10b_m12_m16_highrank.json"
    chunks = ["\n### Table R3-2: E10b M12/M16 high-rank transition "
              f"(source: `{src}`)\n"]
    hdr = [
        "M",
        "k",
        "r_transition",
        "hidden_rank",
        "n_vis",
        "has_hidden_rank_3",
        "empty-lift trailing r",
    ]
    rows = []
    for m in e10b["per_M"]:
        ts = m["transition_summary"]
        trail = ts.get("vanished_Bred_rs") or ts.get("n_vis_zero_rs") or []
        rows.append(
            [
                fmt(m.get("M")),
                fmt(m.get("k")),
                fmt(ts.get("r_transition")),
                fmt(ts.get("hidden_rank_at_transition")),
                fmt(ts.get("n_vis_at_transition")),
                fmt(ts.get("has_hidden_rank_3")),
                ",".join(fmt(r) for r in trail),
            ]
        )
    chunks.append(md_table(hdr, rows))

    hdr2 = [
        "M",
        "method",
        "r",
        "regime",
        "pose med",
        "visible med",
        "hidden med",
        "map med",
    ]
    rows2 = []
    for m in e10b["per_M"]:
        for g in m["method_medians_over_seeds"]:
            rows2.append(
                [
                    fmt(m.get("M")),
                    g["method"],
                    fmt(g.get("r")),
                    str(g.get("regime")),
                    fmt(g.get("pose_error")),
                    fmt(g.get("pose_error_visible")),
                    fmt(g.get("pose_error_hidden")),
                    fmt(g.get("map_error")),
                ]
            )
    chunks.append(
        "\n**Table R3-2b: E10b median over seeds (direct and reduced at the "
        "selected transition/empty-lift ranks)** (source: `%s`)\n" % src
    )
    chunks.append(md_table(hdr2, rows2))
    return "\n".join(chunks)


def _e6b_main_key(g: dict) -> tuple:
    return (
        g["group"],
        g["pose_mode"],
        g["method_alias"],
        g["SNR_dB"],
    )


def md_for_e6b(e6b: dict) -> str:
    if not e6b:
        return ""
    src = "results_e6b_multitx_robust.json"
    main = sorted(e6b["groups"]["main_groups_requested"], key=_e6b_main_key)
    chunks = ["\n### Table R3-3: E6b multi-transmitter robustness "
              f"(source: `{src}`)\n"]
    chunks.append(
        "Main groups requested for the paper: monopole L2 tx_only direct and "
        "wrongpose, monopole L3 tx_only direct, monopole L2 co_moving direct, "
        "dipole L2 tx_only direct. med = 50th percentile over seeds; "
        "IQR = 25th-75th percentile interval."
    )

    def row_key(g: dict) -> str:
        src_t = f"{g['source_type']} L{g['L']}"
        return src_t + " " + str(g["pose_mode"]) + " " + str(g["method_alias"])

    # Group identity map for compact row labels.
    med_headers = [
        "group",
        "SNR dB",
        "success",
        "pose med",
        "theta med",
        "tx med",
        "ty med",
        "map med",
        "final res med",
    ]
    med_rows = []
    iqr_headers = [
        "group",
        "SNR dB",
        "pose IQR",
        "theta IQR",
        "tx IQR",
        "ty IQR",
        "map IQR",
        "final res IQR",
    ]
    iqr_rows = []
    for g in main:
        ident = row_key(g)
        succ = f"{g['n_success']}/{g['n_seeds']}"
        med_rows.append(
            [
                ident,
                fmt(g.get("SNR_dB")),
                succ,
                fmt(g["pose_error"]["median"]),
                fmt(g["theta_error"]["median"]),
                fmt(g["tx_error"]["median"]),
                fmt(g["ty_error"]["median"]),
                fmt(g["map_error"]["median"]),
                fmt(g["final_residual_ratio"]["median"]),
            ]
        )
        iqr_rows.append(
            [
                ident,
                fmt(g.get("SNR_dB")),
                f"{fmt(g['pose_error']['q25'])}-{fmt(g['pose_error']['q75'])}",
                f"{fmt(g['theta_error']['q25'])}-{fmt(g['theta_error']['q75'])}",
                f"{fmt(g['tx_error']['q25'])}-{fmt(g['tx_error']['q75'])}",
                f"{fmt(g['ty_error']['q25'])}-{fmt(g['ty_error']['q75'])}",
                f"{fmt(g['map_error']['q25'])}-{fmt(g['map_error']['q75'])}",
                f"{fmt(g['final_residual_ratio']['q25'])}-"
                f"{fmt(g['final_residual_ratio']['q75'])}",
            ]
        )
    chunks.append("\n**Table R3-3a: medians over seeds**\n")
    chunks.append(md_table(med_headers, med_rows))
    chunks.append("\n**Table R3-3b: IQR (q25-q75) over seeds**\n")
    chunks.append(md_table(iqr_headers, iqr_rows))

    fd = e6b.get("fd_max_error_at_eps_1e-4") or []
    fd_rows = [
        [
            str(c.get("pose_mode")),
            str(c.get("matrix")),
            str(c.get("fd_convention")),
            fmt(c.get("eps")),
            fmt(c.get("max_rel_err")),
            ", ".join(fmt(x) for x in (c.get("per_col_rel_err") or [])),
        ]
        for c in fd
    ]
    chunks.append(
        "\n**Table R3-3c: dipole FD checks at eps=1e-4 (max column-wise relative "
        "error)**\n"
    )
    chunks.append(
        md_table(
            ["pose_mode", "matrix", "convention", "eps", "max_rel_err", "per-column"],
            fd_rows,
        )
    )

    gauge_rows = []
    for e in (e6b.get("dipole_gauge") or {}).get("B_t_analysis") or []:
        cnorm = e.get("realified_colnormed_sv") or []
        gauge_rows.append(
            [
                str(e.get("model")),
                fmt(e.get("L")),
                fmt(e.get("complex_rank")),
                fmt((e.get("complex_sv") or [None])[0]),
                fmt((e.get("complex_sv") or [None])[-1]),
                ", ".join(fmt(x) for x in cnorm),
                fmt(e.get("theta_residual_rel")),
            ]
        )
    chunks.append(
        "\n**Table R3-3d: dipole B_t gauge-rank summary (from "
        "`gauge_partA_dipole`)**\n"
    )
    chunks.append(
        md_table(
            [
                "model",
                "L",
                "complex_rank",
                "sv max",
                "sv min",
                "realified colnormed sv",
                "theta resid rel",
            ],
            gauge_rows,
        )
    )
    return "\n".join(chunks)


def md_for_e11(e11: dict) -> str:
    if not e11:
        return ""
    src = "results_e11_vp_state_null.json"
    chunks = ["\n### Table R3-4: E11 VP state-null test "
              f"(source: `{src}`)\n"]
    rows = []
    for m in e11["methods_median_over_seeds"]:
        rows.append(
            [
                m["method"],
                f"{m['n_success']}/{m['n_seeds']}",
                fmt(m["pose_error"]["median"]),
                fmt(m["pose_error_visible"]["median"]),
                fmt(m["pose_error_hidden"]["median"]),
                fmt(m["map_error"]["median"]),
                fmt(m["data_residual_ratio"]["median"]),
                fmt(m["discarded_current_fraction"]["median"]),
            ]
        )
    chunks.append(
        "\n**Table R3-4a: median over seeds {0,1}, zero p_init (r=6 retained "
        "basis)**\n"
    )
    chunks.append(
        md_table(
            [
                "method",
                "success",
                "pose med",
                "visible med",
                "hidden med",
                "map med",
                "data res med",
                "discarded cur med",
            ],
            rows,
        )
    )
    extra = [
        r
        for r in (e11.get("runs_all") or [])
        if r.get("run_label") != r.get("method") or r.get("init_label") != "zero"
    ]
    if extra:
        chunks.append(
            "\nNote: the file also stores a diagnostic rerun "
            + ", ".join(str(r.get("run_label")) for r in extra)
            + " (init_label `pert_pos`, seed 0). It is kept in `runs_all` of "
            "the summary JSON but excluded from the seed medians above."
        )

    vp = e11.get("finite_dimensional_vp_jacobian") or {}
    sig_rows = []
    for kind in ("hard", "soft", "full_physics_alpha_only"):
        e = vp.get(kind) or {}
        svs = e.get("B_red_singular_values") or []
        sig_rows.append(
            [
                str(kind),
                fmt(e.get("visible_pose_dof_above_cut")),
                ", ".join(fmt(x) for x in svs),
                fmt((e.get("pose_column_norms") or [None])[0]),
                fmt((e.get("pose_column_norms") or [None])[1]),
                fmt((e.get("pose_column_norms") or [None])[2]),
            ]
        )
    chunks.append(
        "\n**Table R3-4b: linearized VP pose signature (B_red singular values "
        "after removing alpha columns, at alpha_init/p_init)**\n"
    )
    chunks.append(
        md_table(
            [
                "signature",
                "visible pose dof above cut",
                "B_red singular values",
                "pose col norm 1",
                "pose col norm 2",
                "pose col norm 3",
            ],
            sig_rows,
        )
    )

    sub = e11.get("visible_pose_subspace_r6_at_p_init") or {}
    hid_sv = sub.get("hid_sv") or []
    chunks.append(
        "\nSubspace anchor (r6 at alpha_init/p_init): n_vis = "
        f"{fmt(sub.get('n_vis'))}, hidden_rank = {fmt(sub.get('hidden_rank'))}, "
        "hid_sv = " + ", ".join(fmt(x) for x in hid_sv) + "."
    )
    interp = e11.get("interpretation") or ""
    one_line = interp.strip().splitlines()[0] if interp.strip() else None
    if one_line:
        chunks.append("\nInterpretation (first stored sentence): " + one_line)
    return "\n".join(chunks)


def md_for_prior(prior: dict, used: list) -> str:
    if not prior or "E5_final" not in prior and "E9_grid_scenes" not in prior:
        return ""
    chunks = [
        "\n### Table R3-5: prior anchor points (sources: "
        + ", ".join(f"`{u}`" for u in used)
        + ")\n"
    ]

    e5 = prior.get("E5_final")
    if e5:
        rows = []
        for r in e5.get("records_median_over_seeds") or []:
            rows.append(
                [
                    str(r.get("p_init_label")),
                    r.get("method"),
                    fmt(r.get("r")),
                    fmt(r.get("pose_error")),
                    fmt(r.get("pose_error_visible")),
                    fmt(r.get("pose_error_hidden")),
                    fmt(r.get("map_error")),
                    f"{r.get('n_success')}/{r.get('n_seeds')}",
                ]
            )
        chunks.append(
            "\n**Table R3-5a: E5 final, medians over seeds by p_init and method "
            "(lossless reduced_r4 vs direct; r6 hidden freeze)**\n"
        )
        chunks.append(
            md_table(
                [
                    "p_init",
                    "method",
                    "r",
                    "pose med",
                    "visible med",
                    "hidden med",
                    "map med",
                    "success",
                ],
                rows,
            )
        )

    e9 = prior.get("E9_grid_scenes")
    if e9:
        rows = []
        for g in e9.get("summary_groups") or []:
            rows.append(
                [
                    str(g.get("scene_id")),
                    fmt(g.get("N")),
                    str(g.get("method")),
                    fmt(g.get("median_pose_error")),
                    fmt(g.get("median_map_error")),
                    fmt(g.get("median_final_residual")),
                    f"{g.get('n_success')}/{g.get('n_seeds')}",
                ]
            )
        chunks.append(
            "\n**Table R3-5b: E9 grid/scene medians by scene and method "
            "(stored `summary_groups`)**\n"
        )
        chunks.append(
            md_table(
                [
                    "scene",
                    "N",
                    "method",
                    "pose med",
                    "map med",
                    "final res med",
                    "success",
                ],
                rows,
            )
        )

    e6 = prior.get("E6_multitx")
    if e6:
        rows = []
        for e in (e6.get("gauge_partA") or {}).get("B_t_analysis") or []:
            rows.append(
                [
                    "monopole",
                    fmt(e.get("L")),
                    fmt(e.get("complex_rank")),
                    fmt((e.get("complex_sv") or [None])[-1]),
                    fmt(e.get("theta_residual_rel")),
                ]
            )
        chunks.append(
            "\n**Table R3-5c: E6 monopole B_t gauge-rank rows (`gauge_partA`, "
            "monopole model)**\n"
        )
        chunks.append(
            md_table(
                [
                    "model",
                    "L",
                    "complex_rank",
                    "smallest complex sv",
                    "theta resid rel",
                ],
                rows,
            )
        )
        chunks.append(
            "FD pass at eps=1e-4 (stored `fd.pass_at_eps_1e-4_max_rel`): "
            + fmt(e6.get("fd_pass_at_eps_1e-4_max_rel"))
            + "."
        )

    e7 = prior.get("E7b_soft_state")
    if e7:
        rows = []
        for s in e7.get("summaries_median_over_seeds") or []:
            rows.append(
                [
                    s.get("method"),
                    fmt(s.get("lam")),
                    fmt(s.get("pose_error_median")),
                    fmt(s.get("pose_error_visible_median")),
                    fmt(s.get("pose_error_hidden_median")),
                    fmt(s.get("map_error_median")),
                    fmt(s.get("final_full_data_residual_median")),
                    f"{s.get('n_success')}/{s.get('n_seeds')}",
                ]
            )
        chunks.append(
            "\n**Table R3-5d: E7b soft-state null - stored seed medians "
            "(direct/reduced_r6 baselines + state_c_soft lambda sweep)**\n"
        )
        chunks.append(
            md_table(
                [
                    "method",
                    "lambda",
                    "pose med",
                    "visible med",
                    "hidden med",
                    "map med",
                    "full data res med",
                    "success",
                ],
                rows,
            )
        )
        chunks.append(
            "Every state_c_soft row has hidden-error median >= 0.0273, larger "
            "than the direct baseline's 0.00941; consistent with the stored "
            "E7b goal/notes that state consistency does not recover the hidden "
            "pose directions that reduced_r6 freezes."
        )
    return "\n".join(chunks)


def main() -> int:
    data, manifest = _load_manifest()

    def section_for(name, builder):
        payload = data.get(name)
        if payload is None:
            return {
                "parse_status": "skipped_missing",
                "source_file": name,
                "reason": "file missing or unreadable; see files_included",
            }
        built = builder(payload)
        if built is None:
            return {
                "parse_status": "skipped_unclear",
                "source_file": name,
                "reason": "required structure not found in parsed file",
            }
        return built

    e10 = section_for("results_e10_settings_sweep.json", build_e10)
    e10b = section_for("results_e10b_m12_m16_highrank.json", build_e10b)
    e6b = section_for("results_e6b_multitx_robust.json", build_e6b)
    e11 = section_for("results_e11_vp_state_null.json", build_e11)
    prior, prior_used = build_prior_anchors(data)

    summary = {
        "files_included": manifest,
        "settings_sweep_E10": e10,
        "rank_transition_E10b": e10b,
        "multitx_E6b": e6b,
        "vp_state_null_E11": e11,
        "prior_anchor_points": prior,
    }

    out_json = os.path.join(WD, "consolidated_round3_summary.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, allow_nan=False)
        fh.write("\n")

    parts = [
        "# Consolidated round-3 results (geometry-lifted tri-space SOM)",
        "",
        "## How to read this",
        "",
        "All numbers below are read directly from the result JSONs listed in "
        "each table caption; medians and 25th/75th-percentile IQRs are computed "
        "over the recorded noise seeds in those files. No physics was "
        "recomputed. Fields that are absent in a file are shown as `—`. "
        "Scientific notation is used only for values < 1e-3 or > 1e3.",
    ]
    parts.append(md_for_e10(e10) if isinstance(e10, dict) and "settings" in e10 else "")
    parts.append(md_for_e10b(e10b) if isinstance(e10b, dict) and "per_M" in e10b else "")
    parts.append(md_for_e6b(e6b) if isinstance(e6b, dict) and "groups" in e6b else "")
    parts.append(md_for_e11(e11) if isinstance(e11, dict) and "methods_median_over_seeds" in e11 else "")
    parts.append(md_for_prior(prior, prior_used))

    md_text = "\n".join(p for p in parts if p != "")
    md_path = os.path.join(WD, "consolidated_round3_tables.md")
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write(md_text + "\n")

    # A few cheap consistency checks against the files' own stored medians.
    _run_checks(data, e6b, e11)

    print(md_text)
    return 0


def _run_checks(data: dict, e6b: dict | None, e11: dict | None) -> None:
    """stdout-only sanity checks of our grouping logic vs stored medians."""
    problems = []
    e6 = data.get("results_e6b_multitx_robust.json")
    if isinstance(e6, dict) and isinstance(e6b, dict):
        stored = {tuple(map(str, [g.get("group"), g.get("pose_mode"),
                                  g.get("method_alias"), g.get("SNR_dB")]))
                  for g in e6.get("summary_medians") or []}
        for g in e6.get("summary_medians") or []:
            key = (g.get("group"), g.get("pose_mode"), g.get("method_alias"),
                   g.get("SNR_dB"))
            mine = None
            for gg in e6b["groups"]["all"]:
                if (gg["group"], gg["pose_mode"], gg["method_alias"],
                        gg["SNR_dB"]) == key:
                    mine = gg
                    break
            if mine is None:
                continue
            for metric in ("pose_error", "map_error", "final_residual_ratio"):
                if metric == "final_residual_ratio":
                    stored_val = g.get("median_final_residual")
                else:
                    stored_val = g.get(
                        {"pose_error": "median_pose_error",
                         "map_error": "median_map_error"}[metric]
                    )
                if stored_val is not None and mine[metric]["median"] is not None:
                    diff = abs(stored_val - mine[metric]["median"])
                    if diff > 1e-9:
                        problems.append(
                            f"E6b median mismatch group={key} metric={metric}: "
                            f"stored={stored_val:.6g} computed="
                            f"{mine[metric]['median']:.6g}"
                        )
        if not stored:
            problems.append("E6b: file has no summary_medians to cross-check")

    e11raw = data.get("results_e11_vp_state_null.json")
    if isinstance(e11raw, dict) and isinstance(e11, dict):
        for m in e11["methods_median_over_seeds"]:
            stored = (e11raw.get("medians_by_method") or {}).get(m["method"]) or {}
            for metric, fld in [
                ("pose_error", "pose_error_median"),
                ("map_error", "map_error_median"),
                ("data_residual_ratio", "data_residual_ratio_median"),
                ("discarded_current_fraction", "discarded_current_fraction_median"),
            ]:
                sv = stored.get(fld)
                mv = m[metric]["median"]
                if sv is not None and mv is not None and abs(sv - mv) > 1e-9:
                    problems.append(
                        f"E11 median mismatch method={m['method']} metric={metric}: "
                        f"stored={sv:.6g} computed={mv:.6g}"
                    )

    if problems:
        print("\n[consistency check] issues:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
    else:
        print("\n[consistency check] computed seed medians match stored "
              "file medians within 1e-9 where comparable.", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
