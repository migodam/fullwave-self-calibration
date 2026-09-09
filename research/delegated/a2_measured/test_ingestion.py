"""Executable checks for the dielTM_dec8f ingestion and model core.

Run with pytest or directly:  python test_ingestion.py
"""

from __future__ import annotations

import numpy as np

import ingestion
import model

EXPECTED_SHA256 = "476cc9d1cfc98797545ab4adf69302dc5aeb45848a24cf8d7b3d222940cc79eb"


def test_raw_structure():
    arr, header = ingestion.read_rows()
    assert arr.shape == (14112, 7), arr.shape
    assert len(header) == 10
    assert np.all(np.isfinite(arr))
    assert ingestion.sha256(ingestion.DATA_FILE) == EXPECTED_SHA256


def test_index_coverage():
    arr, _ = ingestion.read_rows()
    view = arr[:, 0].astype(int)
    rec = arr[:, 1].astype(int)
    freq = arr[:, 2].astype(int)
    assert set(view) == set(range(1, 37))
    assert rec.min() == 1 and rec.max() == 72
    assert set(freq) == set(range(1, 9))
    keys = [tuple(map(int, r[:3])) for r in arr]
    assert len(keys) == len(set(keys)), "duplicate (view, receiver, freq) keys"
    for v in range(1, 37):
        rs = sorted({int(r) for r in rec[view == v]})
        assert len(rs) == 49, (v, len(rs))
        start = 13 + 2 * (v - 1)
        expect = sorted({(start + i - 1) % 72 + 1 for i in range(49)})
        assert rs == expect, (v, rs, expect)


def test_semantics():
    """Column 3 is frequency in GHz; cols 4-7 are total/incident E fields."""
    arr, _ = ingestion.read_rows()
    view = arr[:, 0].astype(int)
    rec = arr[:, 1].astype(int)
    freq = arr[:, 2].astype(int)
    total = arr[:, 3] + 1j * arr[:, 4]
    incident = arr[:, 5] + 1j * arr[:, 6]
    # incident field rotational consistency: view 2 pattern == view 1 shifted
    idx = {}
    for i in range(arr.shape[0]):
        idx[(int(view[i]), int(rec[i]), int(freq[i]))] = incident[i]
    diffs = []
    for fq in range(1, 9):
        for i in range(15, 62):
            if (1, i - 2, fq) in idx and (2, i, fq) in idx:
                e1, e2 = idx[(1, i - 2, fq)], idx[(2, i, fq)]
                diffs.append(abs(e1 - e2) / max(abs(e1), abs(e2)))
    assert np.median(diffs) < 0.01, np.median(diffs)
    # time convention: incident phase coherence with exp(+i k d) dominates
    c0 = 299792458.0
    de, dr = 0.720, 0.760
    sel = (view == 1) & (freq == 1)
    rs = rec[sel]
    th = np.deg2rad((rs - 1) * 5.0)
    pos = np.stack([dr * np.cos(th), dr * np.sin(th)], axis=1)
    src = np.array([de, 0.0])
    d = np.linalg.norm(pos - src[None, :], axis=1)
    k0 = 2.0 * np.pi * 1e9 / c0
    e = incident[sel]
    # exp(+i w t) convention: phase advances as +k d from the emitter
    coh_plus = abs(np.sum(e * np.exp(1j * k0 * d))) / np.sum(abs(e))
    coh_minus = abs(np.sum(e * np.exp(-1j * k0 * d))) / np.sum(abs(e))
    assert coh_plus > coh_minus, (coh_plus, coh_minus)
    # scattered energy is non-trivial
    scat = total - incident
    assert np.sum(np.abs(scat) ** 2) / np.sum(np.abs(incident) ** 2) > 0.01


def test_model_self_checks():
    checks = model.self_checks()
    assert checks["jacobi_anger_max_err"] < 1e-10
    assert checks["bc_derivative_max_rel_err"] < 1e-10
    assert abs(checks["optical_theorem_ratio"] - 1.0) < 1e-10
    assert checks["pec_limit_max_err"] < 1e-3
    assert checks["series_tail_rel_last"] < 1e-10


def _main() -> int:
    for fn in [test_raw_structure, test_index_coverage, test_semantics, test_model_self_checks]:
        fn()
        print(f"PASS {fn.__name__}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
