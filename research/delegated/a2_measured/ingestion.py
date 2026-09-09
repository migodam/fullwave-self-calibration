"""Ingestion and data-QA for the 2001 Fresnel TM dataset (dielTM_dec8f.exp).

Isolated worker (no child agents). Read-only with respect to the raw data in
research/delegated/a2_literature/data/2001_iop_17_6_301/; all outputs under
research/delegated/a2_measured/.

Resolved column semantics (from the 2001 special-section intro, Belkebir &
Saillard, Inverse Problems 17 (2001) 1565-1571, Section 6 "Layout of the
files"):
  col 1: source/view index, 1..36, angle (col1-1)*10 deg  (emitter distance de)
  col 2: receiver index, 1..72, angle (col2-1)*5 deg      (receiver distance dr)
  col 3: operating frequency in GHz (not a sequential index)
  col 4,5: real, imaginary part of the TOTAL electric field
  col 6,7: real, imaginary part of the INCIDENT electric field (target removed)
Time dependence exp(+i*omega*t); TM => E parallel to the cylinder axis.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA_FILE = (
    ROOT.parent
    / "a2_literature"
    / "data"
    / "2001_iop_17_6_301"
    / "dielTM_dec8f.exp"
)

COLUMNS = {
    "view": 0,
    "receiver": 1,
    "freq_ghz": 2,
    "total_re": 3,
    "total_im": 4,
    "incident_re": 5,
    "incident_im": 6,
}


def read_rows(path: Path | str = DATA_FILE) -> tuple[np.ndarray, list[str]]:
    """Parse a 2001 .exp file into a (n,7) float array plus header lines."""
    header: list[str] = []
    rows: list[list[float]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                header.append(stripped)
                continue
            rows.append([float(x) for x in stripped.split()])
    if not rows:
        raise ValueError(f"no data rows in {path}")
    return np.asarray(rows, dtype=float), header


def sha256(path: Path | str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def qa(arr: np.ndarray, header: list[str]) -> dict:
    """Structural QA: shape, counts, finiteness, duplicates, coverage, hash."""
    out: dict = {}
    out["n_rows"] = int(arr.shape[0])
    out["n_columns"] = int(arr.shape[1])
    out["header_lines"] = len(header)

    # Expected geometry from the file header + special-section intro
    expect_views, expect_freqs = 36, 8
    expect_total = expect_views * 49 * expect_freqs
    out["expected_rows_36x49x8"] = expect_total
    out["row_count_matches_expected"] = arr.shape[0] == expect_total

    out["finite_all"] = bool(np.all(np.isfinite(arr)))
    out["n_nonfinite"] = int(np.sum(~np.isfinite(arr)))

    views = arr[:, 0].astype(int)
    recs = arr[:, 1].astype(int)
    freqs = arr[:, 2].astype(int)

    out["view_min_max"] = [int(views.min()), int(views.max())]
    out["receiver_min_max"] = [int(recs.min()), int(recs.max())]
    out["freq_ghz_values"] = sorted({int(f) for f in np.unique(arr[:, 2])})

    out["views_present"] = sorted({int(v) for v in np.unique(views)})
    out["n_views"] = len(out["views_present"])
    out["views_are_1_36"] = out["views_present"] == list(range(1, 37))

    keys = [tuple(map(int, row[:3])) for row in arr]
    uniq = set(keys)
    out["n_unique_view_rec_freq_keys"] = len(uniq)
    out["n_duplicate_keys"] = len(keys) - len(uniq)

    # Per-view receiver window: 49 consecutive receiver indices out of 72
    per_view = {v: sorted({int(r) for r in recs[views == v]}) for v in range(1, 37)}
    out["receivers_per_view"] = {str(v): len(rs) for v, rs in per_view.items()}
    windows = []
    for v, rs in per_view.items():
        if len(rs) == 49:
            # contiguous mod 72 check
            wrap = [rs[0] + 72 if (rs[-1] - rs[0]) > 48 else rs[0]]
            windows.append((v, rs[0], rs[-1]))
    out["receiver_windows"] = [list(w) for w in windows]
    out["all_views_have_49_receivers"] = all(len(rs) == 49 for rs in per_view.values())

    # Frequency coverage per (view, receiver) cell
    cell_freqs = {}
    for k in uniq:
        cell_freqs.setdefault((k[0], k[1]), set()).add(k[2])
    missing = [k for k, fs in cell_freqs.items() if fs != set(range(1, 9))]
    out["cells_with_full_8freq_coverage"] = len(cell_freqs) - len(missing)
    out["cells_missing_freqs"] = len(missing)

    # Value statistics per channel (total, incident), per frequency
    stat = {}
    for fq in range(1, 9):
        sel = freqs == fq
        tot = arr[sel, 3] + 1j * arr[sel, 4]
        inc = arr[sel, 5] + 1j * arr[sel, 6]
        stat[f"f{fq}GHz"] = {
            "n_rows": int(sel.sum()),
            "total_abs_max": float(np.abs(tot).max()),
            "total_abs_mean": float(np.abs(tot).mean()),
            "incident_abs_max": float(np.abs(inc).max()),
            "incident_abs_mean": float(np.abs(inc).mean()),
            "scattered_abs_mean": float(np.abs(tot - inc).mean()),
            "scattered_abs_max": float(np.abs(tot - inc).max()),
        }
    out["per_frequency_stats"] = stat

    # Measured scattered field must be non-trivial (not pure noise)
    scattered = arr[:, 3] + 1j * arr[:, 4] - (arr[:, 5] + 1j * arr[:, 6])
    out["scattered_energy_ratio"] = float(
        np.sum(np.abs(scattered) ** 2) / np.sum(np.abs(arr[:, 5] + 1j * arr[:, 6]) ** 2)
    )
    out["hash_sha256"] = sha256(DATA_FILE)
    return out


def main() -> int:
    arr, header = read_rows()
    results = qa(arr, header)
    out_path = ROOT / "results_ingestion.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, sort_keys=True)
    print(json.dumps(results, indent=2))
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
